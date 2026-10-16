#______________________________________________________
# ADSIM: IMPES solver
#
# IMplicit Pressure, Explicit Species. Implicit advective term: the SUMMED species
# equations are solved for the total concentration C_t^{n+1} first, then each
# species is advanced explicitly with the resulting gradient.
#
# Full derivation in src/IMPES_FORMULATION_NOTES.md (a reorganization of the
# equations of fully_explicit_solver.jl, no new governing equation). Summary:
#
# The ADSIM advective flux is a DIFFUSION operator on the total concentration,
# since the driving pressure depends on that total (P = R T C_t):
#
#     q^a_{α,i} = ∫ (R k C_α T / μ_α) ∇N_i · ∇C_t dΩ
#
# Pressure diffusivity D_P = R k C_α T / μ_α ≈ 2.6e-3 m²/s on the 2 mm carbonation
# mesh, three orders above the molecular 1.6e-6 m²/s. Forward Euler requires
# Δt ∝ h²/D_P: 7.85e-4 s, i.e. 1.7 million steps.
#
# With implicit ∇C_t, every species equation contains the unknown C_t^{n+1} but
# advances only its own C_α, therefore the closed equation for C_t^{n+1} is the sum
# of the species equations, (IMPES-P) below. The explicit remainder δq^a/δC_α is a
# first-order advection operator with an O(θ_g h/|v|) Courant limit.
#
# Own time loop; reuses the free functions of fully_explicit_solver.jl (mass
# assembly, element properties, kinetics closures, stab_tau, VTK output).
#
# Author: Luis Zambrano-Cruzatty
#______________________________________________________

using LinearAlgebra
using Base.Threads
using Printf
using Statistics

# Number of startup steps traced when ADSIM_IMPES_DEBUG=1.
const IMPES_DEBUG_STEPS = 12

# Terms that can set the adaptive step, in counting order. First five: stability and
# positivity bounds. "combined": harmonic sum of the diffusion-like terms, below all
# of them. "growth cap": step recovering towards the limit at dt_growth_max per step.
# "output snap": shortened step ending exactly at a snapshot time.
const DT_LIMIT_NAMES = ("Courant", "diffusion", "conduction", "reaction", "positivity",
                        "combined", "growth cap", "output snap")
const DT_LIMIT_COURANT  = 1
const DT_LIMIT_COMBINED = 6
const DT_LIMIT_GROWTH   = 7
const DT_LIMIT_SNAP     = 8


"""
    assemble_impes_operator!(A_elements, mesh, elem_props, μ_gas, C_g, T, R,
                             chunk_bounds, nchunks) -> A_elements

Assemble the per-element 4×4 operators of the implicit total-concentration equation,

    (A u)_i = ∫ Λ ∇N_i · ∇u dΩ ,    Λ = R k T Σ_α (C_α / μ_α)

`Λ`: total mobility of Eq. (IMPES-P), the species sum of the explicit advective flux
coefficient. The summed equation reproduces the species equations term by term.

`Λ` evaluated **at each Gauss point**, without element averaging. Required for exact
volume consistency: the species step contracts the same Gauss-point coefficients against
∇C_t^{n+1}; any mismatch makes Σ_α C_α^{n+1} drift from C_t^{n+1}. Element averaging
(reuse of cached `K_e`) breaks this.

Quadrature measure from `ShapeFunctions.get_weight` (includes `r_gp` for axisymmetry),
no special case.

Threaded over element chunks, no reduction: one `A_e` per element, no race,
thread-count-independent result.
"""
function assemble_impes_operator!(A_elements::Vector{Matrix{Float64}}, mesh, elem_props,
                                  μ_gas::Vector{Float64}, C_g::Matrix{Float64},
                                  T::Vector{Float64}, R::Float64,
                                  chunk_bounds::Vector{Int}, nchunks::Int)
    NG = length(μ_gas)

    @threads for c in 1:nchunks
        @inbounds for e in (chunk_bounds[c] + 1):chunk_bounds[c + 1]
            n1 = mesh.elements[e, 1]; n2 = mesh.elements[e, 2]
            n3 = mesh.elements[e, 3]; n4 = mesh.elements[e, 4]
            k_intrinsic = elem_props[e].permeability

            A_e = A_elements[e]
            fill!(A_e, 0.0)
            k_intrinsic > 0.0 || continue

            for p in 1:4
                N_p = ShapeFunctions.shape_funcs.N[p]

                T_gp = N_p[1]*T[n1] + N_p[2]*T[n2] + N_p[3]*T[n3] + N_p[4]*T[n4]

                # Σ_α C_α/μ_α at the Gauss point: C_t over the HARMONIC mole-fraction mean
                # viscosity. The reported Darcy velocity uses the arithmetic mean. Known
                # inconsistency of the formulation (~1% for CO2/air), kept for consistency
                # with the explicit solver. See IMPES_FORMULATION_NOTES.md §8.
                mobility = 0.0
                for g in 1:NG
                    C_gp = N_p[1]*C_g[n1, g] + N_p[2]*C_g[n2, g] +
                           N_p[3]*C_g[n3, g] + N_p[4]*C_g[n4, g]
                    mobility += C_gp / μ_gas[g]
                end

                Λ_p = R * k_intrinsic * T_gp * mobility
                Λ_p > 0.0 || continue

                cf = Λ_p * ShapeFunctions.get_weight(e, p) *
                     ShapeFunctions.shape_funcs.gauss_weights[p]
                dN_dx = ShapeFunctions.get_dN_dx(e, p)

                for i in 1:4
                    gi1 = dN_dx[i, 1]; gi2 = dN_dx[i, 2]
                    for j in 1:4
                        A_e[i, j] += cf * (gi1 * dN_dx[j, 1] + gi2 * dN_dx[j, 2])
                    end
                end
            end
        end
    end

    return A_elements
end


"""
    impes_system_diagonal!(diag, mesh, A_elements, M, dt) -> diag

Compute the diagonal of `S = M/Δt + A` for the Jacobi preconditioner. Mutates `diag`.

Serial: four scatter-adds per element is negligible next to the assembly; no per-chunk
buffers, thread-count-independent summation order (same as `apply_element_stiffness!`
in `time_step.jl`).
"""
function impes_system_diagonal!(diag::Vector{Float64}, mesh, A_elements::Vector{Matrix{Float64}},
                                M::Vector{Float64}, dt::Float64)
    @inbounds for i in eachindex(diag)
        diag[i] = M[i] / dt
    end

    @inbounds for e in 1:mesh.num_elements
        A_e = A_elements[e]
        for i in 1:4
            diag[mesh.elements[e, i]] += A_e[i, i]
        end
    end

    return diag
end


"""
    impes_matvec!(y, x, mesh, A_elements, M, dt) -> y

Compute the matrix-free product `y = (M/Δt) ⊙ x + Σ_e A_e x_e`, operator of Eq. (IMPES-P).
Mutates `y`.

`S` is SPD: `M` diagonal and strictly positive (checked at start-up), each `A_e` a sum
of `Λ_p ≥ 0` times outer products (positive semi-definite). Hence conjugate gradients.

Serial (see `impes_system_diagonal!`). Also faster: threading requires per-chunk nodal
buffers with an `nchunks·Nnodes` reduction per matvec, more than the `16·Nelements` flops
of the product.
"""
function impes_matvec!(y::Vector{Float64}, x::Vector{Float64}, mesh,
                       A_elements::Vector{Matrix{Float64}}, M::Vector{Float64}, dt::Float64)
    inv_dt = 1.0 / dt
    @inbounds for i in eachindex(y)
        y[i] = M[i] * inv_dt * x[i]
    end

    @inbounds for e in 1:mesh.num_elements
        A_e = A_elements[e]
        n1 = mesh.elements[e, 1]; n2 = mesh.elements[e, 2]
        n3 = mesh.elements[e, 3]; n4 = mesh.elements[e, 4]
        x1 = x[n1]; x2 = x[n2]; x3 = x[n3]; x4 = x[n4]

        (x1 == 0.0 && x2 == 0.0 && x3 == 0.0 && x4 == 0.0) && continue

        y[n1] += A_e[1,1]*x1 + A_e[1,2]*x2 + A_e[1,3]*x3 + A_e[1,4]*x4
        y[n2] += A_e[2,1]*x1 + A_e[2,2]*x2 + A_e[2,3]*x3 + A_e[2,4]*x4
        y[n3] += A_e[3,1]*x1 + A_e[3,2]*x2 + A_e[3,3]*x3 + A_e[3,4]*x4
        y[n4] += A_e[4,1]*x1 + A_e[4,2]*x2 + A_e[4,3]*x3 + A_e[4,4]*x4
    end

    return y
end


"""
    impes_pcg!(x, b, mesh, A_elements, M, dt, diag, is_dirichlet, workspace;
               tol, maxiter) -> (iterations, relative_residual, converged)

Solve `S x = b` with Jacobi-preconditioned conjugate gradients. Dirichlet nodes handled by
projection (the matrix is never formed). Mutates `x` and the workspace vectors.

Input `x`: full initial guess INCLUDING prescribed Dirichlet values. Output `x`: full
solution. Split `x = x_D + z`, `x_D` the prescribed field (zero off the Dirichlet set),
`z` the unknown (zero on the Dirichlet set):

    P S P z = P (b - S x_D)

`P` zeroes the Dirichlet entries. `P S P` is SPD on the free subspace, CG applies
unchanged, no row or column of `S` modified (matrix-free symmetric elimination).

Convergence on the relative residual over the free subspace. `maxiter` reached:
`converged = false`. The caller must report it, since a partially solved pressure field
gives a physically wrong velocity.
"""
function impes_pcg!(x::Vector{Float64}, b::Vector{Float64}, mesh,
                    A_elements::Vector{Matrix{Float64}}, M::Vector{Float64}, dt::Float64,
                    diag::Vector{Float64}, is_dirichlet::Vector{Bool},
                    r::Vector{Float64}, z::Vector{Float64}, p::Vector{Float64},
                    Ap::Vector{Float64}, x_dir::Vector{Float64}, dx::Vector{Float64};
                    tol::Float64 = 1.0e-10, maxiter::Int = 500)
    N = length(x)

    # Split off the prescribed part and form the effective right-hand side.
    @inbounds for i in 1:N
        x_dir[i] = is_dirichlet[i] ? x[i] : 0.0
    end
    impes_matvec!(Ap, x_dir, mesh, A_elements, M, dt)
    @inbounds for i in 1:N
        r[i] = is_dirichlet[i] ? 0.0 : (b[i] - Ap[i])
    end

    b_norm = sqrt(sum(abs2, r))
    # Empty free subspace, or prescribed field already satisfies every free equation
    if b_norm == 0.0
        @inbounds for i in 1:N
            x[i] = x_dir[i]
        end
        return 0, 0.0, true
    end

    # Warm start from the previous step solution. dx = free correction relative to the
    # incoming guess. Small per-step change, low iteration count.
    @inbounds for i in 1:N
        dx[i] = is_dirichlet[i] ? 0.0 : (x[i] - x_dir[i])
    end
    impes_matvec!(Ap, dx, mesh, A_elements, M, dt)
    @inbounds for i in 1:N
        is_dirichlet[i] ? (r[i] = 0.0) : (r[i] -= Ap[i])
    end

    rz_old = 0.0
    @inbounds for i in 1:N
        z[i] = is_dirichlet[i] ? 0.0 : r[i] / diag[i]
        p[i] = z[i]
        rz_old += r[i] * z[i]
    end

    iters = 0
    res = sqrt(sum(abs2, r)) / b_norm
    converged = res <= tol

    while !converged && iters < maxiter
        iters += 1

        impes_matvec!(Ap, p, mesh, A_elements, M, dt)
        pAp = 0.0
        @inbounds for i in 1:N
            is_dirichlet[i] && (Ap[i] = 0.0)
            pAp += p[i] * Ap[i]
        end
        # Non-positive curvature on an SPD operator: round-off limit reached. Stop
        pAp > 0.0 || break

        α = rz_old / pAp
        @inbounds for i in 1:N
            dx[i] += α * p[i]
            r[i] -= α * Ap[i]
        end

        res = sqrt(sum(abs2, r)) / b_norm
        if res <= tol
            converged = true
            break
        end

        rz_new = 0.0
        @inbounds for i in 1:N
            z[i] = is_dirichlet[i] ? 0.0 : r[i] / diag[i]
            rz_new += r[i] * z[i]
        end
        β = rz_new / rz_old
        rz_old = rz_new

        @inbounds for i in 1:N
            p[i] = z[i] + β * p[i]
        end
    end

    @inbounds for i in 1:N
        x[i] = x_dir[i] + dx[i]
    end

    return iters, res, converged
end


"""
    impes_courant_limit(mesh, elem_props, h_elem, v_gp_cache, NDim) -> Float64

Compute the explicit Courant limit of the species step.

    Δt <= min_e  θ_g,e · h_e / max_p |v_gp(e,p)|

`θ_g` in the numerator: the lumped mass contains `θ_g` and the advective flux does not;
the relevant spectrum is that of `M⁻¹ · (advection)`, `M ∝ θ_g`. Equivalently, `v` is the
Darcy flux per bulk area and `v/θ_g` the pore velocity of the front.

`Inf` for a quiescent field (no limit), e.g. the first step from a uniform state.
"""
function impes_courant_limit(mesh, elem_props, h_elem::Vector{Float64},
                             v_gp_cache::Array{Float64,3}, NDim::Int)
    dt_min = Inf

    @inbounds for e in 1:mesh.num_elements
        θ_g = elem_props[e].θ_g
        θ_g > 0.0 || continue
        h_e = h_elem[e]
        h_e > 0.0 || continue

        v_max2 = 0.0
        for p in 1:4
            s = 0.0
            for d in 1:NDim
                s += v_gp_cache[e, p, d]^2
            end
            v_max2 = max(v_max2, s)
        end
        v_max2 > 0.0 || continue

        dt_e = θ_g * h_e / sqrt(v_max2)
        dt_e < dt_min && (dt_min = dt_e)
    end

    return dt_min
end


"""
    mechanism_dt(report, name) -> Float64

Return the stability limit of one mechanism from the start-up `StabilityReport`.

Diffusion, conduction and reaction limits are recombined each step with the current
Courant limit. Reading them from `calculate_critical_time_step` avoids recomputation and
keeps the adaptive step consistent with the `--report-timestep` output.
"""
function mechanism_dt(report, name::AbstractString)
    report === nothing && return Inf
    for m in report.mechanisms
        if m.name == name
            return m.active ? m.dt : Inf
        end
    end
    return Inf
end


"""
    impes_solver(mesh, materials, calc_params, time_data, project_name, log_print,
                 initial_state=nothing, current_stage=1)

Integrate the ADSIM gas transport, reaction and energy equations with IMPES.

Selected by `[solver] time_integration = "IMPES"`. Same physics as
`fully_explicit_diffusion_solver` with implicit advection. Derivation in
`src/IMPES_FORMULATION_NOTES.md`; shared components in `FULLY_EXPLICIT_SOLVER_README.md`.

# Per step
1. refresh time-dependent boundary conditions
2. assemble the level-n explicit fluxes (diffusion, thermal pressure, gravity, SU) and
   the implicit operator `A^n` from the same Gauss-point coefficients
3. form `b^n`, apply Dirichlet rows, solve `(M/Δt + A^n) C_t^{n+1} = M C_t^n/Δt + b^n`
4. advective flux and Darcy velocities from `∇C_t^{n+1}`
5. Lagrange multipliers, species update, clamp
6. reaction, energy equation, boundary re-imposition, output
7. next Δt from the Courant limit of the current velocity field
"""
function impes_solver(mesh, materials, calc_params, time_data, project_name, log_print,
                      initial_state=nothing, current_stage::Int=1)
    log_print("\n[8/8] Starting IMPES solver (implicit advection)")
    log_print("   Using $(Threads.nthreads()) threads for parallel execution")

    # Access global variables
    global C_g, P, T, v, P_boundary, λ_bc
    global C_lime, C_caco3, C_lime_residual, binder_content, degree_of_carbonation, Caco3_max
    global dC_g_dt, dT_dt, dC_lime_dt
    global M_T, q_cond_T, q_adv_T, q_react_T, q_ext_T, q_heat
    global T_boundary, thermal_node_influences
    global boundary_node_influences, flow_node_influences
    global q_boundary

    negative_conc_warned = Dict{Int, Bool}()

    R = 8.314
    M_caco3 = 100.09 #g/mol
    ρ_caco3 = 2.71e6 #g/m³

    ΔH_r = materials.reactants.reaction_enthalpy
    k_o_reaction = materials.reactants.arrhenius_factor
    E_reaction = materials.reactants.activation_energy
    β_area_reaction = materials.reactants.interfacial_area_beta

    c_p_molar = [materials.gases[name].molar_heat_capacity
                 for name in materials.gas_dictionary]

    Nnodes = mesh.num_nodes
    Nelements = mesh.num_elements
    NGases = length(materials.gas_dictionary)

    # Species properties hoisted out of the hot loop
    μ_gas = [materials.gases[name].dynamic_viscosity for name in materials.gas_dictionary]
    D_gas = [materials.gases[name].diff_coefficient for name in materials.gas_dictionary]
    molar_mass_gas = [materials.gases[name].molar_mass for name in materials.gas_dictionary]

    solver_settings = calc_params["solver_settings"]
    calculate_diffusion = solver_settings["diffusion"] == 1
    calculate_advection = solver_settings["advection"] == 1
    calculate_gravity = solver_settings["gravity"] == 1
    calculate_reaction = solver_settings["reaction_kinetics"] == 1
    calculate_heat_conduction = solver_settings["heat_conduction"] == 1
    calculate_heat_advection = solver_settings["heat_advection"] == 1
    calculate_stabilization = get(solver_settings, "advection_stabilization", 0) == 1
    calculate_heat = calculate_heat_conduction || calculate_heat_advection

    pcg_tol = Float64(get(solver_settings, "impes_pcg_tolerance", 1.0e-10))
    pcg_maxiter = Int(get(solver_settings, "impes_pcg_maxiter", 500))

    # IMPES without advection: pressure system reduces to M/Δt·C_t = M/Δt·C_t + b,
    # i.e. the explicit update at higher cost
    if !calculate_advection
        log_print("   ⚠ WARNING: time_integration = \"IMPES\" with advection = 0. There is no")
        log_print("              advective term to make implicit, so this run is equivalent to")
        log_print("              the explicit solver at additional cost per step.")
    end

    if calculate_advection && !calculate_stabilization
        log_print("   ⚠ WARNING: advection_stabilization = 0. With the pressure feedback made")
        log_print("              implicit, the species step is advection dominated, where")
        log_print("              consistent Galerkin has no discrete maximum principle. Expect")
        log_print("              oscillation at the front and a high C_MIN clamp rate.")
    end

    log_print(@sprintf("   PCG: relative tolerance %.1e, max %d iterations", pcg_tol, pcg_maxiter))

    if calculate_reaction
        log_print(@sprintf("   Reaction enthalpy ΔH_r: %.4g J/mol CO2", ΔH_r))
        log_print(@sprintf("   Arrhenius factor k_o: %.4g m³/(mol·s), activation energy E: %.4g J/mol",
                           k_o_reaction, E_reaction))
        if k_o_reaction <= 0.0
            log_print("Warning: reaction kinetics is enabled but lime_arrhenius_factor = 0, " *
                      "so no carbonation will occur.")
        end
        if ΔH_r > 0.0
            log_print("Warning: reaction_enthalpy = $(ΔH_r) J/mol is positive, i.e. " *
                      "endothermic, so carbonation will cool the soil.")
        end
    end

    # Gravity vector, same convention as the explicit solver: physical acceleration
    # components, g enters v = -(k/μ)(∇p - ρ_g g) with a plus.
    gravity_params = calc_params["gravity"]
    g_vector = [gravity_params["x_component"], gravity_params["y_component"]] *
               gravity_params["magnitude"]

    # Time stepping. dt = STARTING step only, recomputed every step from the velocity field
    # (the governing Courant limit requires a velocity field)
    dt = time_data.actual_dt
    total_time = time_data.total_time
    load_step_time = time_data.time_per_step
    courant_number = time_data.courant_number
    data_saving_interval = calc_params["data_saving_interval"]

    # Solution-independent mechanism limits, from the start-up report (same values as
    # --report-timestep)
    dt_diffusion_fixed = mechanism_dt(time_data.stability, "Diffusive")
    dt_conduction_fixed = mechanism_dt(time_data.stability, "Thermal conduction")
    dt_reaction_fixed = mechanism_dt(time_data.stability, "Reactive")

    # Starting step = explicit-integrator step. Not a stability bound for IMPES, but the
    # correct starting value.
    #
    # Initial conditions may be sharp for the explicit species step: advection_1d starts
    # with a 40 mol/m³ jump across one element. The full diffusion-limited step makes the
    # minority species negative and the C_MIN clamp creates mass. Under IMPES the created
    # mass increases local pressure and velocity, creating more mass (runaway). The explicit
    # solver avoids this with a step three orders smaller.
    #
    # Geometric ramp from this start: a few dozen steps to resolve the startup transient.
    # dt_growth_max: maximum growth per step. 1.1 instead of a larger value: the species step
    # is first order in time, and fast step changes add consistency error. On advection_1d,
    # 1.2 gave 0.25% extra error at the first output compared to 1.1, which costs ~20 extra
    # steps out of a hundred.
    dt_growth_max = 1.1

    dt_advective_start = mechanism_dt(time_data.stability, "Advective [implicit]")
    dt_su_start = mechanism_dt(time_data.stability, "SU Stabilization [implicit]")
    let inv0 = 0.0
        for dtx in (dt_diffusion_fixed, dt_conduction_fixed, dt_advective_start, dt_su_start)
            isfinite(dtx) && dtx > 0.0 && (inv0 += 1.0 / dtx)
        end
        if inv0 > 0.0
            dt = min(dt, courant_number * min(1.0 / inv0, dt_reaction_fixed))
        end
    end
    log_print(@sprintf("   Starting step %.4g s (explicit-equivalent), growing by at most %.2gx per step",
                       dt, dt_growth_max))

    M = zeros(Float64, Nnodes)
    assemble_lumped_mass_vector!(M, mesh, materials)
    if any(M .<= 0.0)
        error("Lumped mass vector contains zero or negative values!")
    end

    K_elements = assemble_element_stiffness_matrices(mesh)
    elem_props = build_element_properties(mesh, materials, C_lime_residual)
    node_owner = build_node_owner_elements(mesh)

    # Element characteristic length, always required: SU parameter (if enabled) and Courant limit
    h_elem = [calculate_element_characteristic_length(mesh, e) for e in 1:Nelements]

    co2_gas_idx = findfirst(name -> name == "CO2", materials.gas_dictionary)

    total_concentration = vec(sum(C_g, dims=2))
    P = total_concentration .* R .* T

    if initial_state !== nothing
        current_time = initial_state.current_time
        output_counter = initial_state.output_counter
        next_output_time = initial_state.next_output_time
        log_print("      Continuing from time: $(current_time) $(calc_params["units"]["time_unit"])")
        log_print("      Next output at: $(next_output_time) $(calc_params["units"]["time_unit"])")
    else
        log_print("      Load step 0 (0.0%)")
        write_output_vtk(mesh, materials, 0, 0.0, project_name, total_concentration)
        current_time = 0.0
        next_output_time = load_step_time
        output_counter = 1
    end

    probe_writer = WriteProbes.init_node_probes(project_name, current_stage, mesh,
                                                materials.gas_dictionary,
                                                calc_params["probing_nodes"],
                                                Float64(data_saving_interval),
                                                calc_params["units"], current_time,
                                                log_print)

    WriteProbes.probe_sample!(probe_writer, current_time, C_g, dC_g_dt,
                              total_concentration, P, T, dT_dt, C_lime, C_caco3,
                              dC_lime_dt, degree_of_carbonation, binder_content, v;
                              force=true)

    # Flow vectors, function-local for type stability (as in the explicit solver).
    # Globals of the same name unused.
    q_diffusion = zeros(Float64, Nnodes, NGases)
    q_advection = zeros(Float64, Nnodes, NGases)
    q_gravitational = zeros(Float64, Nnodes, NGases)
    q_thermal_press = zeros(Float64, Nnodes, NGases)
    # SU stabilization in a separate array (merged into q_advection in the explicit solver).
    # Contracts ∇C_α instead of ∇C_t: per-species EXPLICIT term in b^n. Adding it to the
    # advective flux (assembled after the pressure solve) double counts it.
    # See IMPES_FORMULATION_NOTES.md §4.3 condition 2.
    q_stabilization = zeros(Float64, Nnodes, NGases)
    total_rate = zeros(Float64, Nnodes)
    q_source_sink = zeros(Float64, Nnodes)

    free_conc_sum = zeros(Float64, Nnodes)
    n_free_nodes = zeros(Int, Nnodes)

    nchunks = max(1, min(Nelements, 32))
    chunk_bounds = round.(Int, range(0, Nelements, length = nchunks + 1))

    qd_chunk = [zeros(Float64, Nnodes, NGases) for _ in 1:nchunks]
    qa_chunk = [zeros(Float64, Nnodes, NGases) for _ in 1:nchunks]
    qg_chunk = [zeros(Float64, Nnodes, NGases) for _ in 1:nchunks]
    qT_chunk = [zeros(Float64, Nnodes, NGases) for _ in 1:nchunks]
    qs_chunk = [zeros(Float64, Nnodes, NGases) for _ in 1:nchunks]
    v_chunk = [zeros(Float64, Nnodes, NDim) for _ in 1:nchunks]

    q_aux_chunk = [zeros(4) for _ in 1:nchunks]
    qT_aux_chunk = [zeros(4) for _ in 1:nchunks]
    qs_aux_chunk = [zeros(4) for _ in 1:nchunks]
    ρ_g_chunk = [zeros(4) for _ in 1:nchunks]
    a_su_chunk = [zeros(4) for _ in 1:nchunks]
    aT_su_chunk = [zeros(4) for _ in 1:nchunks]

    Q_dot_T = zeros(Float64, Nnodes)
    MT_chunk = [zeros(Float64, Nnodes) for _ in 1:nchunks]
    qd_T_chunk = [zeros(Float64, Nnodes) for _ in 1:nchunks]
    qr_T_chunk = [zeros(Float64, Nnodes) for _ in 1:nchunks]
    qa_T_chunk = [zeros(Float64, Nnodes) for _ in 1:nchunks]
    qh_chunk = [zeros(Float64, Nnodes, NDim) for _ in 1:nchunks]
    qhw_chunk = [zeros(Float64, Nnodes) for _ in 1:nchunks]
    heat_flux_weight = zeros(Float64, Nnodes)

    v_gp_cache = zeros(Float64, Nelements, 4, NDim)

    #______________________________________________________
    # IMPES pressure system: per-element operators and the CG workspace.
    #______________________________________________________
    A_elements = [zeros(Float64, 4, 4) for _ in 1:Nelements]
    S_diag = zeros(Float64, Nnodes)
    C_t_new = zeros(Float64, Nnodes)
    rhs = zeros(Float64, Nnodes)
    is_dirichlet = fill(false, Nnodes)
    pcg_r = zeros(Float64, Nnodes)
    pcg_z = zeros(Float64, Nnodes)
    pcg_p = zeros(Float64, Nnodes)
    pcg_Ap = zeros(Float64, Nnodes)
    pcg_xdir = zeros(Float64, Nnodes)
    pcg_dx = zeros(Float64, Nnodes)

    # Nodes whose total concentration is fixed by a Dirichlet row of (IMPES-P).
    #
    # Concentration or partial-pressure BC nodes have P_boundary zeroed for EVERY gas (both
    # apply_*_bc! loop over all species). Partially prescribed nodes (inconsistent summed
    # equation) cannot occur through those dictionaries, therefore the mask is exact.
    #
    # Pressure BC nodes added at solve time (value may follow a curve). Nodes with both:
    # pressure BC takes precedence (same end-of-step order as the explicit solver).
    fully_prescribed = [all(P_boundary[i, g] == 0 for g in 1:NGases) for i in 1:Nnodes]
    pressure_bc_nodes = collect(keys(mesh.absolute_pressure_bc))
    for i in 1:Nnodes
        is_dirichlet[i] = fully_prescribed[i]
    end
    for j in pressure_bc_nodes
        is_dirichlet[j] = true
    end
    log_print(@sprintf("   Pressure system: %d nodes, %d Dirichlet (%d pressure BC, %d prescribed composition)",
                       Nnodes, count(is_dirichlet), length(pressure_bc_nodes),
                       count(fully_prescribed)))

    #______________________________________________________
    # Time-dependent boundary conditions
    #______________________________________________________
    time_functions = calc_params["time_functions"]

    transient_flow = [(node_id, g, ids[g], mesh.uniform_flow_bc[node_id][g])
                      for (node_id, ids) in mesh.uniform_flow_bc_tf
                      for g in 1:min(length(ids), NGases)
                      if ids[g] != 0 && haskey(mesh.uniform_flow_bc, node_id)]

    transient_conc = [(node_id, g, ids[g], mesh.concentration_bc[node_id][g])
                      for (node_id, ids) in mesh.concentration_bc_tf
                      for g in 1:min(length(ids), NGases)
                      if ids[g] != 0 && haskey(mesh.concentration_bc, node_id)]

    transient_temp = [(node_id, id, mesh.temperature_bc[node_id])
                      for (node_id, id) in mesh.temperature_bc_tf
                      if id != 0 && haskey(mesh.temperature_bc, node_id)]

    has_transient_flow = !isempty(transient_flow)
    has_transient_conc = !isempty(transient_conc)
    has_transient_temp = !isempty(transient_temp)

    save_data = false
    negative_lime_warned = false

    # Clamp accounting, split into the two populations of the §5.1 positivity limit.
    #
    #   :roundoff - node at or below C_floor before the update: round-off skipped by the
    #               positivity bound (see C_floor guard). A species fully displaced by a
    #               front stays here for the rest of the run, one count per step per node.
    #   :material - node with real gas driven negative. Target of the "arresting an
    #               unstable mode" warning; created mass feeds the §5.1 pressure feedback.
    #
    # Both bins use the SAME C_floor as the positivity limit.
    negative_conc_count = zeros(Int, NGases)          # material
    negative_conc_count_roundoff = zeros(Int, NGases)
    # Moles created by resetting negative concentrations to zero, per gas (amount, measured
    # by the volume-consistency drift)
    clamped_mass = zeros(Float64, NGases)
    clamped_mass_roundoff = zeros(Float64, NGases)
    # TEMPORARY DIAGNOSTIC: remove after the clamp-location investigation is resolved.
    # Per-node, per-gas material-clamp count to classify top offenders against the BC
    # dictionaries.
    negative_conc_node_count = zeros(Int, Nnodes, NGases)
    negative_conc_node_roundoff = zeros(Int, Nnodes, NGases)
    # Reaction throttle, second known source of §4.3 volume-consistency error: the pressure
    # RHS uses the unthrottled rate, since the throttle depends on the post-solve advective
    # rate. Indexed by gas, one writer per entry in the threaded species loop.
    throttle_count = zeros(Int, NGases)
    throttled_mass = zeros(Float64, NGases)
    last_output_step = 0

    # Adaptive-step and solver-health diagnostics, reset at every output.
    pcg_iters_total = 0
    pcg_iters_max = 0
    pcg_res_max = 0.0
    pcg_solves = 0
    dt_min_seen = Inf
    dt_max_seen = 0.0
    consistency_max = 0.0
    dt_limit_idx = DT_LIMIT_GROWTH
    dt_limit_count = zeros(Int, length(DT_LIMIT_NAMES))

    # Reaction limit on the current field. dt_reaction_fixed (start-up report): lowest
    # mesh-file temperature and t = 0 lime inventory. Sizes the first step, then replaced
    # each step by the node-loop value.
    #
    # Both values typically agree for most of a run. λ is largest where the specimen is
    # coldest and least reacted, i.e. a region not yet reached by the front (start-up
    # assumption). The value changes once that region warms or reacts, and immediately
    # when the coldest point is a convective boundary (start-up uses initial and
    # prescribed temperatures only).
    rxn_lambda = zeros(Float64, Nnodes)
    dt_reaction = dt_reaction_fixed

    # Trial-step control. positivity_safety: margin of the positivity bound, applied to the
    # current step. Retry cap bounds the cost of a step with changing rates (exhaustion:
    # see accept/reject block).
    positivity_safety = 0.5
    positivity_max_retries = 8
    step_rejected = 0       # trials thrown away since the last output
    step_exhausted = 0      # steps that ran out of retries and fell back to the throttle

    # Assigned in the trial loop, read after it: scope-level locals
    t_next = current_time
    dt_positivity = Inf
    C_floor = 0.0
    dC_g_dt_rxn_store = zeros(Float64, Nnodes)

    # Absolute time tolerance for loop bound and output snapping, scaled with run length
    t_eps = 1.0e-12 * max(total_time, 1.0)

    # End of THIS stage. total_simulation_time = stage duration. The explicit solver uses a
    # step count (num_steps_per_load × num_load_steps) from the checkpoint time, e.g. resume
    # at t = 2000 with total_simulation_time = 2000 ends at t = 4000. The adaptive loop has
    # no step count, therefore the end time is explicit. Bounding by total_time would give
    # zero steps for a resumed stage.
    stage_start_time = current_time
    stage_end_time = current_time + total_time

    impes_debug = get(ENV, "ADSIM_IMPES_DEBUG", "0") != "0"

    step = 0
    while current_time < stage_end_time - t_eps
        step += 1

        #______________________________________________________
        # Refresh the time-dependent boundary conditions for this step.
        # A Neumann flux ENTERS the rate -> evaluate at t_n.
        # A Dirichlet value is IMPOSED on the new state -> evaluate at t_{n+1}.
        #______________________________________________________
        #______________________________________________________
        # Trial step, §5.1.
        #
        # Positivity bound unknown before the step: species rates -> advective flux ->
        # C_t^{n+1} from the implicit solve -> dt. Trial step, repeated at the reported
        # bound when too long.
        #
        # Whole body repeated, including the solve: dt enters the accumulation term of
        # (IMPES-P), therefore C_t^{n+1} is tied to its dt. A different species dt breaks
        # Σ_α C_α^{n+1} = C_t^{n+1}, the identity of the IMPES split.
        #______________________________________________________
        trial = 0
        while true
            trial += 1

            t_next = current_time + dt

            if has_transient_flow
                for (node_id, gas_idx, tf_id, q_const) in transient_flow
                    q_boundary[node_id, gas_idx] =
                        bc_value(time_functions, tf_id, q_const, current_time, dt) *
                        flow_node_influences[node_id]
                end
            end

            # Level-n total concentration, the accumulation term of (IMPES-P).
            for i in 1:Nnodes
                s = 0.0
                for g in 1:NGases
                    s += C_g[i, g]
                end
                total_concentration[i] = s
            end

            #reset flow vectors (q_boundary is prefilled and not reset here)
            q_diffusion .= 0.0
            q_advection .= 0.0
            q_gravitational .= 0.0
            q_thermal_press .= 0.0
            q_stabilization .= 0.0
            total_rate .= 0.0
            q_source_sink .= 0.0

            #______________________________________________________
            # Reaction kinetics, level n. Before the fluxes (q_source_sink enters b^n)
            #______________________________________________________
            if calculate_reaction && co2_gas_idx !== nothing
                @threads for node_id in 1:Nnodes
                    rxn_lambda[node_id] = 0.0
                    e = node_owner[node_id]
                    e == 0 && continue
                    props = elem_props[e]

                    # Total gas concentration for the area factor; C_aq uses the CO2 component
                    C_g_tot = 0.0
                    for g in 1:NGases
                        C_g_tot += C_g[node_id, g]
                    end

                    r = extent_of_reaction_rate(C_g[node_id, co2_gas_idx], C_g_tot,
                                                C_lime[node_id],
                                                props.residual_lime, props.θ_w, T[node_id],
                                                k_o_reaction, E_reaction, β_area_reaction)

                    dC_lime_dt[node_id] = -props.θ_w * r

                    if props.θ_g > 0.0
                        q_source_sink[node_id] = -M[node_id] * (props.θ_w / props.θ_g) * r
                    else
                        q_source_sink[node_id] = 0.0
                    end

                    # Reaction eigenvalue of the nodal CO2 equation at the current state. Two
                    # exponentials per node. Recomputed each step: lime consumption, heating,
                    # and cold faces at convective boundaries change every factor.
                    rxn_lambda[node_id] =
                        reaction_decay_constant(C_g_tot, C_lime[node_id], props.residual_lime,
                                                props.θ_w, props.θ_g, T[node_id],
                                                k_o_reaction, E_reaction, β_area_reaction)
                end

                # Factor 1/2 as in the start-up bound: forward Euler on decay rate λ is stable
                # to 2/λ and positive to 1/λ
                λ_rxn_max = maximum(rxn_lambda)
                dt_reaction = λ_rxn_max > 0.0 ? 1.0 / (2.0 * λ_rxn_max) : Inf
            end

            # Total mobility Λ frozen at level n (plain IMPES, linear pressure equation).
            # Lag justification and relinearization cost: IMPES_FORMULATION_NOTES.md §8.
            pcg_it = 0
            pcg_res = 0.0
            pcg_ok = true

            for c in 1:nchunks
                qd_chunk[c] .= 0.0
                qa_chunk[c] .= 0.0
                qg_chunk[c] .= 0.0
                qT_chunk[c] .= 0.0
                qs_chunk[c] .= 0.0
            end

            #______________________________________________________
            # Level-n explicit fluxes: diffusion, thermal pressure, gravity, SU.
            #
            # Advective flux NOT assembled here (requires C_t^{n+1}). All terms in this loop
            # enter b^n.
            #______________________________________________________
            @threads for c in 1:nchunks
                qd_c = qd_chunk[c]
                qg_c = qg_chunk[c]
                qT_c = qT_chunk[c]
                qs_c = qs_chunk[c]
                q_aux = q_aux_chunk[c]
                qT_aux = qT_aux_chunk[c]
                qs_aux = qs_aux_chunk[c]
                ρ_g_buf = ρ_g_chunk[c]
                a_su = a_su_chunk[c]
                aT_su = aT_su_chunk[c]

                for e in (chunk_bounds[c] + 1):chunk_bounds[c + 1]
                    nodes = mesh.elements[e, :]
                    props = elem_props[e]
                    θ_g = props.θ_g
                    τ = props.tortuosity
                    k_intrinsic = props.permeability

                    C_t = [total_concentration[nodes[i]] for i in 1:4]
                    T_e = [T[nodes[i]] for i in 1:4]

                    for gas_idx in 1:NGases
                        D_g = D_gas[gas_idx]
                        μ_g = μ_gas[gas_idx]
                        D_g_eff = calculate_diffusion ? (θ_g * D_g * τ) : 0.0

                        C_e = [C_g[nodes[i], gas_idx] for i in 1:4]

                        # Diffusion, Eq. (diffusive_flux)
                        if calculate_diffusion
                            mul!(q_aux, K_elements[e], C_e, θ_g * D_g * τ, 0.0)
                            for i in 1:4
                                qd_c[nodes[i], gas_idx] += q_aux[i]
                            end
                        end

                        # Thermal pressure and SU stabilization, both level n.
                        if calculate_advection
                            fill!(qT_aux, 0.0)
                            fill!(qs_aux, 0.0)

                            for p in 1:4
                                N_p = ShapeFunctions.shape_funcs.N[p]
                                C_gp = N_p' * C_e
                                T_gp = N_p' * T_e
                                C_tot_gp = N_p' * C_t
                                dV = ShapeFunctions.get_weight(e, p)
                                Wp = ShapeFunctions.shape_funcs.gauss_weights[p]
                                dN_dx = ShapeFunctions.get_dN_dx(e, p)

                                # Thermal contribution to the pressure gradient,
                                # Eq. (thermal_pressure_flux). Explicit, since ∇T is advanced
                                # explicitly by the energy equation.
                                qT_aux .+= (R * k_intrinsic * C_gp * C_tot_gp * dV * Wp / μ_g) .*
                                           (dN_dx * (dN_dx' * T_e))

                                # Streamline-upwind stabilization. τ* velocities at ∇C_t^n
                                # (term in b^n, formed before the pressure solve). Lag in the
                                # stabilization PARAMETER only; keeps the summed equation exactly
                                # consistent with the species equations.
                                # See IMPES_FORMULATION_NOTES.md §8.
                                if calculate_stabilization
                                    v_gp = (-k_intrinsic * R * T_gp / μ_g) .* (dN_dx' * C_t)
                                    v_gpT = (-k_intrinsic * R * C_tot_gp / μ_g) .* (dN_dx' * T_e)

                                    τ_g = stab_tau(v_gp, h_elem[e], D_g_eff)
                                    τ_gT = stab_tau(v_gpT, h_elem[e], D_g_eff)

                                    mul!(a_su, dN_dx, v_gp)
                                    qs_aux .+= (τ_g * dV * Wp) .* (a_su .* (a_su' * C_e))

                                    #Stabilizes the species concentration along the thermally driven
                                    #streamline, therefore C_e contracted with v_gpT. Contracting T_e adds no
                                    #diagonal weight to the C_g^i system (no discrete maximum principle) and
                                    #has units of K m^3/s instead of a molar flux.
                                    mul!(aT_su, dN_dx, v_gpT)
                                    qs_aux .+= (τ_gT * dV * Wp) .* (aT_su .* (aT_su' * C_e))
                                end
                            end

                            for i in 1:4
                                qT_c[nodes[i], gas_idx] += qT_aux[i]
                                qs_c[nodes[i], gas_idx] += qs_aux[i]
                            end
                        end

                        # Gravity, Eq. (gravitational_flux)
                        if calculate_gravity
                            fill!(q_aux, 0.0)
                            for p in 1:4
                                N_p = ShapeFunctions.shape_funcs.N[p]
                                C_gp = N_p' * C_e
                                dV = ShapeFunctions.get_weight(e, p)
                                Wp = ShapeFunctions.shape_funcs.gauss_weights[p]
                                dN_dx = ShapeFunctions.get_dN_dx(e, p)

                                fill!(ρ_g_buf, 0.0)
                                for i in 1:4
                                    for g in 1:NGases
                                        ρ_g_buf[i] += C_g[nodes[i], g] * molar_mass_gas[g]
                                    end
                                end
                                ρ_g_gp = N_p' * ρ_g_buf

                                q_aux .-= (k_intrinsic * C_gp * dV * Wp * ρ_g_gp / μ_g) .*
                                          (dN_dx * g_vector)
                            end
                            for i in 1:4
                                qg_c[nodes[i], gas_idx] += q_aux[i]
                            end
                        end
                    end
                end
            end

            for c in 1:nchunks
                q_diffusion .+= qd_chunk[c]
                q_gravitational .+= qg_chunk[c]
                q_thermal_press .+= qT_chunk[c]
                q_stabilization .+= qs_chunk[c]
            end

            #______________________________________________________
            # (IMPES-P): assemble A^n and solve for C_t^{n+1}.
            #______________________________________________________
            if calculate_advection
                assemble_impes_operator!(A_elements, mesh, elem_props, μ_gas, C_g, T, R,
                                         chunk_bounds, nchunks)
            else
                for e in 1:Nelements
                    fill!(A_elements[e], 0.0)
                end
            end
            impes_system_diagonal!(S_diag, mesh, A_elements, M, dt)

            # b^n = Σ_α [q^b - q^d - q^T - q^g - q^su]_α + q^r + M C_t^n/Δt
            #
            # P_boundary gates each species as in the explicit rate. Fully prescribed node:
            # zero contribution (Dirichlet row)
            for i in 1:Nnodes
                acc = 0.0
                for g in 1:NGases
                    acc += (q_boundary[i, g] - q_diffusion[i, g] - q_thermal_press[i, g] -
                            q_gravitational[i, g] - q_stabilization[i, g]) * P_boundary[i, g]
                end
                if calculate_reaction && co2_gas_idx !== nothing
                    # Unthrottled: the throttle depends on the post-solve advective rate. Where
                    # the throttle acts, §4.3 volume consistency error = throttled amount
                    # (measured by the consistency diagnostic below).
                    acc += q_source_sink[i] * P_boundary[i, co2_gas_idx]
                end
                rhs[i] = acc + M[i] * total_concentration[i] / dt
            end

            # Dirichlet values. Composition-prescribed nodes: given total. Pressure BC nodes:
            # t_{n+1} value (same target as the Lagrange multipliers).
            for i in 1:Nnodes
                C_t_new[i] = total_concentration[i]
            end
            for j in pressure_bc_nodes
                tf_id = get(mesh.absolute_pressure_tf, j, 0)
                P_bc = bc_value(time_functions, tf_id, mesh.absolute_pressure_bc[j], t_next)
                C_t_new[j] = P_bc / (R * T[j])
            end

            pcg_it, pcg_res, pcg_ok =
                impes_pcg!(C_t_new, rhs, mesh, A_elements, M, dt, S_diag, is_dirichlet,
                           pcg_r, pcg_z, pcg_p, pcg_Ap, pcg_xdir, pcg_dx;
                           tol = pcg_tol, maxiter = pcg_maxiter)

            pcg_solves += 1
            pcg_iters_total += pcg_it
            pcg_iters_max = max(pcg_iters_max, pcg_it)
            pcg_res_max = max(pcg_res_max, pcg_res)

            if !pcg_ok
                log_print(@sprintf("Warning: IMPES pressure solve did not converge at step %d (%d iterations, relative residual %.3e > %.1e).",
                                   step, pcg_it, pcg_res, pcg_tol))
                log_print("         The velocity field this step is not trustworthy. Lower the " *
                          "Courant number or raise impes_pcg_maxiter.")
            end

            if any(isnan, C_t_new)
                error("Simulation failed: NaN in the IMPES pressure solution at step $step")
            end

            #______________________________________________________
            # Advective flux from ∇C_t^{n+1}, and the Darcy velocities.
            #
            # One loop, same element gradient. C_gp and T_gp are the same interpolations used
            # for Λ in assemble_impes_operator!, therefore Σ_α C_α^{n+1} = C_t^{n+1} exactly.
            #______________________________________________________
            for c in 1:nchunks
                qa_chunk[c] .= 0.0
                v_chunk[c] .= 0.0
            end
            q_advection .= 0.0
            v .= 0.0

            # Nodal pressure at level n+1: velocity from a single ∇P^{n+1} (no product-rule split).
            # T^n used (energy equation not yet advanced); difference O(dt dT/dt), below the
            # scheme error.
            for i in 1:Nnodes
                P[i] = C_t_new[i] * R * T[i]
            end

            @threads for c in 1:nchunks
                qa_c = qa_chunk[c]
                v_c = v_chunk[c]
                q_aux = q_aux_chunk[c]
                ρ_g_vel = ρ_g_chunk[c]

                for e in (chunk_bounds[c] + 1):chunk_bounds[c + 1]
                    nodes = mesh.elements[e, :]
                    props = elem_props[e]
                    k_intrinsic = props.permeability

                    C_t_e = [C_t_new[nodes[i]] for i in 1:4]
                    T_e = [T[nodes[i]] for i in 1:4]
                    P_e = [P[nodes[i]] for i in 1:4]

                    # Advective flux with the implicit gradient, Eq. (IMPES-S)
                    if calculate_advection
                        for gas_idx in 1:NGases
                            μ_g = μ_gas[gas_idx]
                            C_e = [C_g[nodes[i], gas_idx] for i in 1:4]
                            fill!(q_aux, 0.0)

                            for p in 1:4
                                N_p = ShapeFunctions.shape_funcs.N[p]
                                C_gp = N_p' * C_e
                                T_gp = N_p' * T_e
                                dV = ShapeFunctions.get_weight(e, p)
                                Wp = ShapeFunctions.shape_funcs.gauss_weights[p]
                                dN_dx = ShapeFunctions.get_dN_dx(e, p)

                                q_aux .+= (R * k_intrinsic * C_gp * T_gp * dV * Wp / μ_g) .*
                                          (dN_dx * (dN_dx' * C_t_e))
                            end

                            for i in 1:4
                                qa_c[nodes[i], gas_idx] += q_aux[i]
                            end
                        end
                    end

                    # Darcy velocities, Eq. (Darcy_law). Same convention as the explicit
                    # solver (arithmetic composition-weighted viscosity)
                    for p in 1:4
                        N_p = ShapeFunctions.shape_funcs.N[p]
                        dV = ShapeFunctions.get_weight(e, p)
                        w = ShapeFunctions.shape_funcs.gauss_weights[p]
                        dN_dx = ShapeFunctions.get_dN_dx(e, p)

                        grad_P = dN_dx' * P_e
                        C_total_gp = N_p' * C_t_e

                        C_TOL = 1e-12
                        μ_g_weighted = 0.0
                        if C_total_gp > C_TOL
                            for g in 1:NGases
                                C_g_gp = 0.0
                                for i in 1:4
                                    C_g_gp += N_p[i] * C_g[nodes[i], g]
                                end
                                μ_g_weighted += (C_g_gp / C_total_gp) * μ_gas[g]
                            end
                        else
                            μ_g_weighted = mean(μ_gas)
                        end

                        v_gp = -(k_intrinsic / μ_g_weighted) * grad_P

                        mass_weight = props.θ_g * w * dV
                        for i in 1:4
                            v_c[nodes[i], :] += v_gp * N_p[i] * mass_weight
                        end

                        if calculate_gravity
                            ρ_g = ρ_g_vel
                            fill!(ρ_g, 0.0)
                            for i in 1:4
                                for g in 1:NGases
                                    ρ_g[i] += C_g[nodes[i], g] * molar_mass_gas[g]
                                end
                            end
                            ρ_g_gp = N_p' * ρ_g

                            v_g_gp = (k_intrinsic / μ_g_weighted) * ρ_g_gp * g_vector
                            for i in 1:4
                                v_c[nodes[i], :] += v_g_gp * N_p[i] * mass_weight
                            end
                            v_gp = v_gp + v_g_gp
                        end

                        for d in 1:NDim
                            v_gp_cache[e, p, d] = v_gp[d]
                        end
                    end
                end
            end

            for c in 1:nchunks
                q_advection .+= qa_chunk[c]
                v .+= v_chunk[c]
            end

            for i in 1:Nnodes
                if M[i] > 0.0
                    v[i, :] ./= M[i]
                end
            end

            #______________________________________________________
            # Species rates, Eq. (IMPES-S).
            #
            # UNTHROTTLED. The positivity test requires the physical rate: a throttled rate is
            # exactly -C/dt and returns the current step, giving no information to the
            # controller. Reaction part stored separately for the throttle on the accepted step.
            #______________________________________________________
            @threads for gas_idx in 1:NGases
                is_co2 = calculate_reaction && gas_idx == co2_gas_idx

                for i in 1:Nnodes
                    dC_g_dt[i, gas_idx] = ((q_boundary[i, gas_idx] - q_diffusion[i, gas_idx] -
                                            q_advection[i, gas_idx] - q_thermal_press[i, gas_idx] -
                                            q_gravitational[i, gas_idx] - q_stabilization[i, gas_idx]) *
                                           P_boundary[i, gas_idx]) / M[i]

                    if is_co2
                        rxn = (q_source_sink[i] * P_boundary[i, gas_idx]) / M[i]
                        dC_g_dt_rxn_store[i] = rxn
                        dC_g_dt[i, gas_idx] += rxn
                    end
                end
            end

            #______________________________________________________
            # Positivity limit, §5.1, and accept/reject test of the trial.
            #
            #     Δt <= min over (i,α) with dC_α/dt < 0 of  C_α / (-dC_α/dt)
            #
            # Nonlinear bound, independent of the Courant limit. Under IMPES the C_MIN clamp
            # creates mass when resetting a negative species to zero, and the pressure
            # equation responds immediately: higher total, pressure, velocity, overshoot.
            # Positive feedback, therefore the step must keep the clamp inactive.
            #
            # Nearly empty nodes skipped: a small negative rate there is round-off and would
            # collapse the step. Clamp still active for those nodes.
            #
            # C_floor outside the loop: same threshold for the clamp bins below. Previously
            # the limit skipped a node as round-off while the counter scored it as overshoot,
            # and runs with mostly displaced-species jitter were reported as unstable.
            #______________________________________________________
            C_floor = 1.0e-6 * max(maximum(C_g), 1.0)
            dt_positivity = Inf
            for gas_idx in 1:NGases
                for i in 1:Nnodes
                    rate = dC_g_dt[i, gas_idx]
                    rate < 0.0 || continue
                    C_g[i, gas_idx] > C_floor || continue
                    d = C_g[i, gas_idx] / (-rate)
                    d < dt_positivity && (dt_positivity = d)
                end
            end

            # Accept, or reject the trial and repeat with a shorter step.
            #
            # Applying the bound to the NEXT step cannot prevent the measured violation (step
            # already applied, mass already created). Rejection costs one extra pressure solve
            # and keeps species non-negative by construction.
            #
            # Replacement step = the bound itself instead of halving: overshoot by 2x or 50x
            # both resolved in one retry. positivity_safety = 0.5: margin for rate changes
            # over the step.
            #
            # Repeated failure: rates change with the shorter step. Retry count capped; on
            # exhaustion the step proceeds with throttle and clamp, and is counted for the log.
            if isfinite(dt_positivity) && dt > positivity_safety * dt_positivity
                if trial <= positivity_max_retries
                    dt = positivity_safety * dt_positivity
                    step_rejected += 1
                    continue
                end
                step_exhausted += 1
            end
            break
        end  # trial step

        #______________________________________________________
        # Reaction throttle, backstop of §5.1.
        #
        # Expected inactive with the trial loop; the output-block count verifies it. The
        # withheld amount equals the mismatch between Σ_α C_α^{n+1} and C_t^{n+1} (pressure
        # RHS uses the unthrottled rate), giving the volume-consistency warning a measured value.
        #______________________________________________________
        if calculate_reaction && co2_gas_idx !== nothing
            for i in 1:Nnodes
                rxn = dC_g_dt_rxn_store[i]
                rxn < 0.0 || continue
                C_g[i, co2_gas_idx] + dt * dC_g_dt[i, co2_gas_idx] < 0.0 || continue

                transport = dC_g_dt[i, co2_gas_idx] - rxn
                allowed = min(-C_g[i, co2_gas_idx] / dt - transport, 0.0)
                scale = clamp(allowed / rxn, 0.0, 1.0)

                throttle_count[co2_gas_idx] += 1
                throttled_mass[co2_gas_idx] += -rxn * (1.0 - scale) * dt * M[i]
                dC_g_dt[i, co2_gas_idx] = transport + rxn * scale
                dC_lime_dt[i] *= scale
            end
        end

        for i in 1:Nnodes
            for gas_idx in 1:NGases
                total_rate[i] += dC_g_dt[i, gas_idx]
            end
        end

        #______________________________________________________
        # Lagrangian multipliers at the pressure boundary, Eq. (§4.5).
        #
        # Same as the explicit solver. At a pressure BC node the (IMPES-P) row is REPLACED by
        # the Dirichlet value, therefore exact consistency does not apply and the species are
        # redistributed to match the prescribed total. Residual M(C_t^n + dt·total_rate - target)/dt
        # = that mismatch. Composition-weighted split prevents negative minority species at
        # the outlet.
        #______________________________________________________
        @threads for j in pressure_bc_nodes
            if haskey(boundary_node_influences, j)
                C_total_j = sum(C_g[j, :])

                n_free = 0
                free_sum = 0.0
                for g in 1:NGases
                    n_free += P_boundary[j, g]
                    free_sum += P_boundary[j, g] * C_g[j, g]
                end
                free_conc_sum[j] = free_sum
                n_free_nodes[j] = n_free

                if n_free > 0
                    λ_bc[j] = M[j] * (C_total_j + dt * total_rate[j] - C_t_new[j]) / dt
                else
                    λ_bc[j] = 0.0
                end
            else
                λ_bc[j] = 0.0
                free_conc_sum[j] = 0.0
                n_free_nodes[j] = 0
            end
        end

        #______________________________________________________
        # Species update and clamp.
        #______________________________________________________
        @threads for gas_idx in 1:NGases
            gas_name = materials.gas_dictionary[gas_idx]

            for i in 1:Nnodes
                lagrangian_correction = 0.0
                if haskey(boundary_node_influences, i) && haskey(mesh.absolute_pressure_bc, i)
                    C_TOL = 1e-12
                    denom = free_conc_sum[i]
                    weight = denom > C_TOL ? (C_g[i, gas_idx] / denom) :
                             (n_free_nodes[i] > 0 ? 1.0 / n_free_nodes[i] : 0.0)
                    lagrangian_correction = (λ_bc[i] / M[i]) * weight * P_boundary[i, gas_idx]
                end

                # Pre-update value, to classify a clamp as material or round-off
                C_before = C_g[i, gas_idx]

                C_g[i, gas_idx] += dt * (dC_g_dt[i, gas_idx] - lagrangian_correction)

                C_MIN = 1e-12
                if C_g[i, gas_idx] < C_MIN
                    if C_g[i, gas_idx] < 0.0
                        # Moles created by the reset. M[i] = lumped nodal mass, comparable to
                        # the domain inventory below.
                        created = -C_g[i, gas_idx] * M[i]

                        if C_before > C_floor
                            # Real gas driven negative (§5.1 population)
                            negative_conc_count[gas_idx] += 1
                            clamped_mass[gas_idx] += created
                            negative_conc_node_count[i, gas_idx] += 1
                            if !get(negative_conc_warned, gas_idx, false)
                                log_print("Warning: Negative concentration detected for gas $gas_name at step $step. Setting to zero.")
                                negative_conc_warned[gas_idx] = true
                            end
                        else
                            # Node already empty within the positivity-limit tolerance
                            negative_conc_count_roundoff[gas_idx] += 1
                            clamped_mass_roundoff[gas_idx] += created
                            negative_conc_node_roundoff[i, gas_idx] += 1
                        end
                    end
                    C_g[i, gas_idx] = 0.0
                end
            end

            if any(isnan.(C_g[:, gas_idx]))
                nan_nodes = findall(isnan.(C_g[:, gas_idx]))
                log_print("ERROR: NaN detected in gas $gas_name at step $step, nodes: $nan_nodes")
                error("Simulation failed due to NaN values")
            end
        end

        #______________________________________________________
        # Volume-consistency diagnostic, Eq. (§4.3).
        #
        # Sum of species updates must equal the implicit total at every non-Dirichlet node.
        # Nonzero value = defect (Step 1 and Step 2 Gauss-point coefficients differ), EXCEPT
        # where the reaction throttle or the C_MIN clamp acted. Relative to the local total.
        #______________________________________________________
        for i in 1:Nnodes
            is_dirichlet[i] && continue
            s = 0.0
            for g in 1:NGases
                s += C_g[i, g]
            end
            scale = max(abs(C_t_new[i]), 1.0e-12)
            consistency_max = max(consistency_max, abs(s - C_t_new[i]) / scale)
        end

        # Update reaction kinetic terms for lime concentration
        if calculate_reaction
            for i in 1:Nnodes
                ΔC_lime = dt * dC_lime_dt[i]
                if C_lime[i] + ΔC_lime < 0.0
                    if !negative_lime_warned
                        log_print("Warning: Negative lime concentration detected at step $step. Clamping to zero.")
                        negative_lime_warned = true
                    end
                    ΔC_lime = -C_lime[i]
                end
                C_lime[i] += ΔC_lime
                C_caco3[i] += -ΔC_lime
                binder_content[i] = C_caco3[i] * M_caco3 / ρ_caco3
                degree_of_carbonation[i] = Caco3_max[i] > 0.0 ? C_caco3[i] / Caco3_max[i] : 0.0
            end
        end

        #______________________________________________________
        # Energy equation, Eq. (FEM_compact_energy):
        #
        #     M^L_T · dT/dt = q̃_ext - q̃_d - q̃_a + q̃_r
        #
        # Same as the explicit solver. v_gp_cache contains the level-(n+1) velocity under
        # IMPES, i.e. the velocity used to transport the species.
        #______________________________________________________
        if calculate_heat
            ρ_w = materials.liquid.density
            c_w = materials.liquid.specific_heat

            if calculate_reaction
                @threads for i in 1:Nnodes
                    Q_dot_T[i] = ΔH_r * dC_lime_dt[i]
                end
            else
                Q_dot_T .= 0.0
            end

            M_T .= 0.0
            q_cond_T .= 0.0
            q_adv_T .= 0.0
            q_react_T .= 0.0
            q_heat .= 0.0
            heat_flux_weight .= 0.0
            for c in 1:nchunks
                MT_chunk[c] .= 0.0
                qd_T_chunk[c] .= 0.0
                qa_T_chunk[c] .= 0.0
                qr_T_chunk[c] .= 0.0
                qh_chunk[c] .= 0.0
                qhw_chunk[c] .= 0.0
            end

            @threads for c in 1:nchunks
                MT_c = MT_chunk[c]
                qd_c = qd_T_chunk[c]
                qa_c_T = qa_T_chunk[c]
                qr_c = qr_T_chunk[c]
                qh_c = qh_chunk[c]
                qhw_c = qhw_chunk[c]

                for e in (chunk_bounds[c] + 1):chunk_bounds[c + 1]
                    nodes = mesh.elements[e, :]
                    props = elem_props[e]
                    ρ_s = props.specific_gravity * ρ_w
                    solid_water = (1.0 - props.porosity) * ρ_s * props.specific_heat_solids +
                                  props.θ_w * ρ_w * c_w

                    T_e = [T[nodes[i]] for i in 1:4]

                    for p in 1:4
                        N_p = ShapeFunctions.shape_funcs.N[p]
                        Wp = ShapeFunctions.shape_funcs.gauss_weights[p]
                        dV = ShapeFunctions.get_weight(e, p) * Wp

                        ρc_g_gp = 0.0
                        for g in 1:NGases
                            C_gp = 0.0
                            for i in 1:4
                                C_gp += N_p[i] * C_g[nodes[i], g]
                            end
                            ρc_g_gp += C_gp * c_p_molar[g]
                        end
                        C_mix_gp = solid_water + props.θ_g * ρc_g_gp

                        for i in 1:4
                            MT_c[nodes[i]] += C_mix_gp * N_p[i] * dV
                        end

                        if calculate_heat_conduction && props.λ_e > 0.0
                            dN_dx = ShapeFunctions.get_dN_dx(e, p)
                            gradT = dN_dx' * T_e
                            qd_gp = dN_dx * gradT
                            for i in 1:4
                                qd_c[nodes[i]] += props.λ_e * qd_gp[i] * dV
                            end
                        end

                        # Advective flux, Eq. (thermal_advective_flux), NON-CONSERVATIVE form
                        # (ρc)_g v·∇T. Same as the explicit solver (see there for why the
                        # divergence form ∇·((ρc)_g v T) is inconsistent with this left-hand side).
                        if calculate_heat_advection
                            ρc_gp = ρc_g_gp
                            dN_dx_a = ShapeFunctions.get_dN_dx(e, p)
                            gradT_a = dN_dx_a' * T_e          # [NDim]
                            vg = @view v_gp_cache[e, p, :]
                            v_gradT = 0.0
                            for d in 1:NDim
                                v_gradT += vg[d] * gradT_a[d]
                            end
                            for i in 1:4
                                qa_c_T[nodes[i]] += ρc_gp * v_gradT * N_p[i] * dV
                            end
                        end

                        dN_dx_h = ShapeFunctions.get_dN_dx(e, p)
                        gradT_h = dN_dx_h' * T_e
                        T_gp_h = 0.0
                        for i in 1:4
                            T_gp_h += N_p[i] * T_e[i]
                        end
                        for i in 1:4
                            nd = nodes[i]
                            wN = N_p[i] * dV
                            qhw_c[nd] += wN
                            for d in 1:NDim
                                q_flux = -props.λ_e * gradT_h[d]
                                if calculate_heat_advection
                                    q_flux += ρc_g_gp * T_gp_h * v_gp_cache[e, p, d]
                                end
                                qh_c[nd, d] += q_flux * wN
                            end
                        end

                        if calculate_reaction
                            Q_gp = 0.0
                            for i in 1:4
                                Q_gp += N_p[i] * Q_dot_T[nodes[i]]
                            end
                            for i in 1:4
                                qr_c[nodes[i]] += Q_gp * N_p[i] * dV
                            end
                        end
                    end
                end
            end

            for c in 1:nchunks
                M_T .+= MT_chunk[c]
                q_cond_T .+= qd_T_chunk[c]
                q_adv_T .+= qa_T_chunk[c]
                q_react_T .+= qr_T_chunk[c]
                q_heat .+= qh_chunk[c]
                heat_flux_weight .+= qhw_chunk[c]
            end

            for i in 1:Nnodes
                if heat_flux_weight[i] > 0.0
                    for d in 1:NDim
                        q_heat[i, d] /= heat_flux_weight[i]
                    end
                end
            end

            q_ext_T .= 0.0
            for (node_id, params) in mesh.convective_heat_bc
                h_c = params[1]
                T_inf = params[2]
                le = get(thermal_node_influences, node_id, 0.0)
                q_ext_T[node_id] -= h_c * (T[node_id] - T_inf) * le
            end

            @threads for i in 1:Nnodes
                if M_T[i] > 0.0
                    dT_dt[i] = (q_ext_T[i] - q_cond_T[i] - q_adv_T[i] + q_react_T[i]) / M_T[i]
                else
                    dT_dt[i] = 0.0
                end
                dT_dt[i] *= T_boundary[i]
            end

            for i in 1:Nnodes
                T[i] += dt * dT_dt[i]
                if T[i] < 0.0
                    log_print("Warning: Temperature below absolute zero at node $i. Setting to 0.0 K.")
                    T[i] = 0.0
                end
            end
        end

        # Prescribed temperatures that follow a curve.
        if has_transient_temp
            for (node_id, tf_id, T_const) in transient_temp
                T[node_id] = bc_value(time_functions, tf_id, T_const, t_next)
            end
        end

        # Partial pressure boundaries, re-imposed against the updated temperature.
        for (node_id, partial_pressures) in mesh.partial_pressure_bc
            tf_ids = get(mesh.partial_pressure_bc_tf, node_id, Int[])
            for gas_idx in 1:NGases
                tf_id = gas_idx <= length(tf_ids) ? tf_ids[gas_idx] : 0
                p_partial = bc_value(time_functions, tf_id, partial_pressures[gas_idx], t_next)
                C_g[node_id, gas_idx] = p_partial / (R * T[node_id])
            end
        end

        # Prescribed concentrations that follow a curve.
        if has_transient_conc
            for (node_id, gas_idx, tf_id, C_const) in transient_conc
                C_g[node_id, gas_idx] = bc_value(time_functions, tf_id, C_const, t_next)
            end
        end

        total_concentration = vec(sum(C_g, dims=2))

        # Re-impose the total at pressure BC nodes with the updated temperature: reported
        # pressure = boundary value (instead of the T^n solve value)
        for node_id in pressure_bc_nodes
            tf_id = get(mesh.absolute_pressure_tf, node_id, 0)
            P_bc = bc_value(time_functions, tf_id, mesh.absolute_pressure_bc[node_id], t_next)
            total_concentration[node_id] = P_bc / (R * T[node_id])
        end

        P = total_concentration .* R .* T

        current_time += dt

        WriteProbes.probe_sample!(probe_writer, current_time, C_g, dC_g_dt,
                                  total_concentration, P, T, dT_dt, C_lime, C_caco3,
                                  dC_lime_dt, degree_of_carbonation, binder_content, v)

        #______________________________________________________
        # Next time step, Eq. (§5).
        #
        # Courant limit governs IMPES runs, from the current velocity field. Diffusion,
        # conduction and reaction limits: solution independent, from the start-up report.
        #
        # SU stabilization adds a term of the same size as the Courant one (τ*|v|² ≤ h|v|/2
        # gives a diffusive limit θ_g h/|v|), counted separately in the harmonic sum (same
        # as the explicit path).
        #______________________________________________________
        dt_courant = calculate_advection ?
            impes_courant_limit(mesh, elem_props, h_elem, v_gp_cache, NDim) : Inf

        inv_sum = 0.0
        for dtx in (dt_diffusion_fixed, dt_conduction_fixed, dt_courant,
                    calculate_stabilization ? dt_courant : Inf)
            isfinite(dtx) && dtx > 0.0 && (inv_sum += 1.0 / dtx)
        end
        dt_explicit = inv_sum > 0.0 ? 1.0 / inv_sum : Inf
        # 0.5 instead of the Courant number: consume at most half the available amount
        # (margin for rate changes over the step)
        dt_target = min(courant_number * min(dt_explicit, dt_reaction),
                        0.5 * dt_positivity)

        # Term that set the step. Candidates include their factors, therefore the minimum
        # is the term that produced dt_target.
        #
        # Recorded per step, reported as a share of the output interval. A single overwritten
        # name reports "growth cap" after every reduction and hides the term causing it.
        dt_limit_idx = DT_LIMIT_COURANT
        let cands = (courant_number * dt_courant,
                     courant_number * dt_diffusion_fixed,
                     courant_number * dt_conduction_fixed,
                     courant_number * dt_reaction,
                     0.5 * dt_positivity)
            best = Inf
            for k in eachindex(cands)
                if cands[k] < best
                    best = cands[k]
                    dt_limit_idx = k
                end
            end
            # Harmonic sum may be below every individual term
            if courant_number * dt_explicit < 0.95 * best
                best = courant_number * dt_explicit
                dt_limit_idx = DT_LIMIT_COMBINED
            end
            # Growth cap active during the opening ramp and after any reduction. Its share per
            # interval = fraction of the run spent recovering instead of at the limit.
            if dt_growth_max * dt < 0.95 * best
                dt_limit_idx = DT_LIMIT_GROWTH
            end
        end

        # No bound (e.g. pure Darcy run at rest): fall back to the output interval
        isfinite(dt_target) || (dt_target = Float64(load_step_time))

        # Growth limiter: prevents step oscillation after transient velocity spikes (accuracy
        # loss at the front)
        dt_target = min(dt_target, dt_growth_max * dt)

        dt_min_seen = min(dt_min_seen, dt)
        dt_max_seen = max(dt_max_seen, dt)

        # Per-step trace of the startup transient, enabled with ADSIM_IMPES_DEBUG=1.
        #
        # First IMPES_DEBUG_STEPS steps only: every failure observed so far started there
        # (ramping step, discontinuous initial condition; divergence visible by step four).
        # Output block covers the rest of the run.
        if impes_debug && step <= IMPES_DEBUG_STEPS
            vmax = 0.0
            for e in 1:Nelements, p in 1:4
                s = 0.0
                for d in 1:NDim; s += v_gp_cache[e, p, d]^2; end
                vmax = max(vmax, sqrt(s))
            end
            dCmax = 0.0; gmax = 0
            for g in 1:NGases, i in 1:Nnodes
                if abs(dC_g_dt[i, g]) * dt > dCmax
                    dCmax = abs(dC_g_dt[i, g]) * dt; gmax = g
                end
            end
            log_print(@sprintf("   [dbg] step %3d t=%.4g dt=%.4g (%s) |v|max=%.4g dt_courant=%.4g maxΔC=%.4g (gas %d) max|ΣC-C_t|=%.4g",
                               step, current_time, dt, dt_limit_name, vmax, dt_courant,
                               dCmax, gmax, maximum(abs.(C_t_new .- total_concentration))))
        end

        # Check if we need to save output
        if save_data || current_time >= stage_end_time - t_eps

            progress = 100.0 * (current_time - stage_start_time) / total_time
            log_print(@sprintf("      Load Step %d (%.1f%%), Time = %.4e %s",
                              output_counter, progress, current_time,
                              calc_params["units"]["time_unit"]))

            steps_since = max(1, step - last_output_step)
            # Rejected trials (one pressure solve each) shown next to the step count, only if any
            retry_note = step_rejected > 0 ?
                @sprintf(" (+%d rejected%s)", step_rejected,
                         step_exhausted > 0 ? @sprintf(", %d unresolved", step_exhausted) : "") : ""
            log_print(@sprintf("        Steps %d%s, dt %.4e to %.4e %s",
                               steps_since, retry_note, dt_min_seen, dt_max_seen,
                               calc_params["units"]["time_unit"]))

            # Share of each step-setting term over the interval, terms below 1% dropped
            let tot = sum(dt_limit_count)
                if tot > 0
                    parts = String[]
                    for k in sortperm(dt_limit_count, rev = true)
                        share = 100.0 * dt_limit_count[k] / tot
                        share < 1.0 && break
                        push!(parts, @sprintf("%s %.0f%%", DT_LIMIT_NAMES[k], share))
                    end
                    isempty(parts) || log_print("        Step set by: " * join(parts, ", "))
                end
            end
            fill!(dt_limit_count, 0)
            step_rejected = 0
            step_exhausted = 0
            if pcg_solves > 0
                log_print(@sprintf("        PCG: %.1f iterations mean, %d max, residual %.2e max",
                                   pcg_iters_total / pcg_solves, pcg_iters_max, pcg_res_max))
            end
            log_print(@sprintf("        Volume consistency |ΣC_α − C_t|/C_t: %.2e max", consistency_max))
            # Floor set by the pressure solve: relative residual pcg_tol gives a solution error
            # ~pcg_tol times the condition number (1e-8 on a stiff mesh is expected). Warning
            # only well above that floor.
            # Above the floor, report both known causes with their measured size. Both ~0 with
            # PCG at tolerance: operator and species flux Gauss-point coefficients differ
            # (solver defect).
            if consistency_max > 1.0e-6
                log_print(@sprintf("        ⚠ above PCG floor: throttle %d (%.2e mol), clamp %.2e mol",
                                   sum(throttle_count), sum(throttled_mass),
                                   sum(clamped_mass) + sum(clamped_mass_roundoff)))
            end
            pcg_iters_total = 0; pcg_iters_max = 0; pcg_res_max = 0.0; pcg_solves = 0
            dt_min_seen = Inf; dt_max_seen = 0.0; consistency_max = 0.0

            # Clamp report. "real": node with gas driven negative, the only bin indicating an
            # unstable mode. "empty": below C_floor (round-off, skipped by the step limit);
            # fully displaced species remain there for the rest of the run.
            total_material = sum(negative_conc_count)
            total_roundoff = sum(negative_conc_count_roundoff)
            if total_material + total_roundoff > 0
                rate = total_material / (steps_since * Nnodes * NGases)
                # Per-gas split only with more than one active gas
                active = [g for g in 1:NGases
                          if negative_conc_count[g] + negative_conc_count_roundoff[g] > 0]
                gas_split = length(active) > 1 ?
                    " (" * join((@sprintf("%s %d/%d", materials.gas_dictionary[g],
                                          negative_conc_count[g],
                                          negative_conc_count_roundoff[g])
                                 for g in active), ", ") * ")" :
                    (isempty(active) ? "" : " ($(materials.gas_dictionary[active[1]]))")
                log_print(@sprintf("        Clamp real/empty: %d/%d%s, created %.2e mol",
                                   total_material, total_roundoff, gas_split,
                                   sum(clamped_mass) + sum(clamped_mass_roundoff)))

                if rate > 1e-3
                    log_print(@sprintf("        ⚠ real clamp %.2f%% of nodal updates - unstable mode, not round-off.", 100 * rate))
                    log_print("          Lower the Courant number and check advection_stabilization.")
                end

                # TEMPORARY DIAGNOSTIC: remove after the clamp-location investigation is
                # resolved. Hits on pressure/composition-BC nodes indicate the Lagrangian
                # boundary correction; hits along the front do not.
                # Ranked on the real bin; empty bin as fallback to locate the front.
                node_material = vec(sum(negative_conc_node_count, dims = 2))
                node_roundoff = vec(sum(negative_conc_node_roundoff, dims = 2))
                rank_by = total_material > 0 ? node_material : node_roundoff
                tops = String[]
                for nid in sortperm(rank_by, rev = true)[1:min(6, Nnodes)]
                    rank_by[nid] == 0 && break
                    f = String[]
                    haskey(mesh.absolute_pressure_bc, nid) && push!(f, "P")
                    haskey(mesh.concentration_bc, nid) && push!(f, "C")
                    haskey(mesh.partial_pressure_bc, nid) && push!(f, "pp")
                    haskey(mesh.convective_heat_bc, nid) && push!(f, "heat")
                    haskey(mesh.uniform_flow_bc, nid) && push!(f, "flow")
                    isempty(f) && push!(f, "int")
                    push!(tops, @sprintf("%d:%d/%d[%s]", nid, node_material[nid],
                                         node_roundoff[nid], join(f, ",")))
                end
                isempty(tops) ||
                    log_print("          top " * (total_material > 0 ? "real" : "empty") *
                              ": " * join(tops, " "))
                fill!(negative_conc_node_count, 0)
                fill!(negative_conc_node_roundoff, 0)
            end
            fill!(negative_conc_count, 0)
            fill!(negative_conc_count_roundoff, 0)
            fill!(clamped_mass, 0.0)
            fill!(clamped_mass_roundoff, 0.0)
            fill!(throttle_count, 0)
            fill!(throttled_mass, 0.0)
            last_output_step = step

            write_output_vtk(mesh, materials, output_counter, current_time, project_name, total_concentration)

            next_output_time += load_step_time
            output_counter += 1
            save_data = false
            negative_lime_warned = false
        end

        # Stage done. Exit before sizing a zero-length step (would trigger the collapse guard)
        current_time >= stage_end_time - t_eps && break

        # Snap to the next output time or stage end. Stage-end cap replaces the fixed step
        # count of the explicit solver.
        t_snap = min(next_output_time, stage_end_time)
        if current_time + dt_target > t_snap - t_eps
            dt = t_snap - current_time
            # Stage end writes output through the loop condition, therefore only an output
            # boundary sets save_data.
            save_data = t_snap == next_output_time
            # Snap overrides the physics limit; counted separately such that a short
            # remainder is not reported as a stability event.
            dt_limit_idx = DT_LIMIT_SNAP
        else
            dt = dt_target
            save_data = false
        end
        dt_limit_count[dt_limit_idx] += 1

        # Collapsed step: error instead of an infinite loop
        if dt <= 0.0 || !isfinite(dt)
            error("IMPES time step collapsed to $(dt) at step $step, t = $(current_time). " *
                  "This usually means the velocity field diverged; check the PCG " *
                  "convergence warnings above.")
        end
    end

    log_print("   ✓ Time integration completed")
    log_print(@sprintf("   ✓ Final time: %.4e %s", current_time, calc_params["units"]["time_unit"]))
    log_print(@sprintf("   ✓ Total steps taken: %d", step))

    WriteProbes.close_probes!(probe_writer, log_print)

    return (current_time=current_time, output_counter=output_counter, next_output_time=next_output_time)
end
