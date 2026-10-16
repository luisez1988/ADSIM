#------------------------------------------------------------------------------
# ADSIM Time Step Module
# Critical time step calculation and time stepping for ADSIM FEM calculations
#------------------------------------------------------------------------------

"""
MechanismStability

Per-mechanism stability diagnostics, one entry per assembled spatial operator.

# Fields
- `name::String`: Mechanism name in the log ("Diffusive", "Advective", ...)
- `active::Bool`: Whether the solver assembles this term
- `λ_measured::Float64`: Largest eigenvalue of M̂⁻¹K from power iteration [1/s]
- `λ_analytic::Float64`: Largest eigenvalue from the closed form, `4·coef/h_min²` [1/s]
- `dt::Float64`: Stability limit of this mechanism only, `2/λ` [s]

`λ_measured` and `λ_analytic` are `NaN` for mechanisms without a spatial operator (reaction,
a nodal ODE) and for inactive mechanisms.
"""
struct MechanismStability
    name::String
    active::Bool
    λ_measured::Float64
    λ_analytic::Float64
    dt::Float64
end


"""
StabilityReport

Diagnostics of the critical time step calculation, reported by `--report-timestep` and
the run log.

`measured`: whether power iteration supplied the step. False only with an uninitialized
shape function cache: closed forms used, axisymmetric correction NOT applied (see
`measure_operator_eigenvalues`).
"""
struct StabilityReport
    mechanisms::Vector{MechanismStability}
    measured::Bool
    axisymmetric::Bool
    h_min::Float64
    C_max::Float64
    critical_dt::Float64
    limiting_scale::String
end


"""
TimeStepData

Time stepping information of the simulation.

# Fields
- `critical_dt::Float64`: Critical time step calculated from stability criteria [s]
- `actual_dt::Float64`: Actual time step used (critical_dt × Courant number) [s]
- `total_time::Float64`: Total simulation time requested by user [s]
- `time_per_step::Float64`: Time per load step (when data is saved) [s]
- `num_steps::Int`: Number of time steps required for the simulation
- `num_steps_per_load::Int`: Number of time steps per load step
- `num_load_steps::Int`: Number of load steps in the simulation
- `courant_number::Float64`: Courant number (0 < C_N ≤ 1)
- `h_min::Float64`: Minimum characteristic element size [m]
- `stability::Union{StabilityReport, Nothing}`: Diagnostics behind `critical_dt`
"""
mutable struct TimeStepData
    critical_dt::Float64
    actual_dt::Float64
    total_time::Float64
    time_per_step::Float64
    num_steps::Int
    num_steps_per_load::Int
    num_load_steps::Int
    courant_number::Float64
    h_min::Float64
    stability::Union{StabilityReport, Nothing}

    function TimeStepData()
        new(0.0, 0.0, 0.0, 0.0, 0, 0, 0, 0.0, 0.0, nothing)
    end
end


"""
    calculate_element_characteristic_length(mesh::MeshData, elem_id::Int) -> Float64

Compute the characteristic length h [m] of a quadrilateral element (1-based `elem_id`).

# Formula
Minimum of:
- Distance between midpoints of sides 1-2 and 3-4
- Distance between midpoints of sides 2-3 and 4-1
"""
function calculate_element_characteristic_length(mesh, elem_id::Int)
    # Get the four nodes of the element
    nodes = mesh.elements[elem_id, :]
    
    # Get coordinates of all four nodes
    x1, y1 = mesh.coordinates[nodes[1], 1], mesh.coordinates[nodes[1], 2]
    x2, y2 = mesh.coordinates[nodes[2], 1], mesh.coordinates[nodes[2], 2]
    x3, y3 = mesh.coordinates[nodes[3], 1], mesh.coordinates[nodes[3], 2]
    x4, y4 = mesh.coordinates[nodes[4], 1], mesh.coordinates[nodes[4], 2]
    
    # Calculate midpoints of opposite sides
    # Side 1-2 midpoint
    mx12 = (x1 + x2) / 2.0
    my12 = (y1 + y2) / 2.0
    
    # Side 3-4 midpoint
    mx34 = (x3 + x4) / 2.0
    my34 = (y3 + y4) / 2.0
    
    # Side 2-3 midpoint
    mx23 = (x2 + x3) / 2.0
    my23 = (y2 + y3) / 2.0
    
    # Side 4-1 midpoint
    mx41 = (x4 + x1) / 2.0
    my41 = (y4 + y1) / 2.0
    
    # Distance between opposite side midpoints
    d1 = sqrt((mx34 - mx12)^2 + (my34 - my12)^2)
    d2 = sqrt((mx41 - mx23)^2 + (my41 - my23)^2)
    
    # Return minimum distance as characteristic length
    return min(d1, d2)
end


"""
    find_minimum_characteristic_length(mesh::MeshData) -> Float64

Return the minimum element characteristic length h_min [m] of the mesh.
"""
function find_minimum_characteristic_length(mesh)
    h_min = Inf
    
    for elem_id in 1:mesh.num_elements
        h = calculate_element_characteristic_length(mesh, elem_id)
        h_min = min(h_min, h)
    end
    
    return h_min
end


"""
    get_maximum_diffusion_coefficient(materials) -> Float64

Return the maximum gas diffusion coefficient D_max [m²/s].
"""
function get_maximum_diffusion_coefficient(materials)
    D_max = 0.0
    
    for gas_name in materials.gas_dictionary
        gas = materials.gases[gas_name]
        D_max = max(D_max, gas.diff_coefficient)
    end
    
    return D_max
end


"""
    get_minimum_gas_viscosity(materials) -> Float64

Return the minimum gas dynamic viscosity μ_g [Pa·s].
"""
function get_minimum_gas_viscosity(materials)
    μ_min = Inf
    
    for gas_name in materials.gas_dictionary
        gas = materials.gases[gas_name]
        μ_min = min(μ_min, gas.dynamic_viscosity)
    end
    
    return μ_min
end


"""
    get_maximum_total_concentration(mesh::MeshData, NGases::Int) -> Float64

Return the largest TOTAL gas concentration [mol/m^3] over the initial condition and the
prescribed boundaries.

The advective flux, Eq. (advective_flux), contracts against the total concentration,
therefore the total scales the operator. Identical to the per-species maximum for one
gas; for more gases the per-species maximum makes the advective bound too permissive by
the ratio of total to largest species.

Boundary values included: a prescribed inlet can exceed the initial field (e.g. inlet at
40 mol/m^3, interior maximum 33).
"""
function get_maximum_total_concentration(mesh, NGases::Int;
                                         time_functions = nothing, t_end::Float64 = 0.0)
    C_max = 0.0

    for (_, concentrations) in mesh.initial_concentrations
        total = 0.0
        for gas_idx in 1:NGases
            total += concentrations[gas_idx]
        end
        C_max = max(C_max, total)
    end

    for (_, concentrations) in mesh.concentration_bc
        total = 0.0
        for gas_idx in 1:NGases
            total += concentrations[gas_idx]
        end
        C_max = max(C_max, total)
    end

    # Pressure-specified boundaries also prescribe a concentration (ideal gas law in
    # apply_partial_pressure_bc! and at end-of-step re-imposition). Carbonation benchmark:
    # interior 41 mol/m^3 of air, injection face CO2 at 113 kPa = 46.4 mol/m^3. Without
    # them the advective bound is 13% too permissive where the gradients are steepest.
    #
    # LOWEST temperature in the problem: C = P/(R T) is largest there (upper bound).
    R = 8.314
    T_conv = get_minimum_reference_temperature(mesh)
    if T_conv > 0.0
        for (_, partial_pressures) in mesh.partial_pressure_bc
            total = 0.0
            for gas_idx in 1:min(NGases, length(partial_pressures))
                total += partial_pressures[gas_idx]
            end
            C_max = max(C_max, total / (R * T_conv))
        end

        for (node_id, pressure) in mesh.absolute_pressure_bc
            # Time-driven pressure: bounded over the whole stage instead of sampled at t = 0.
            # A ramp from atmospheric to injection pressure would otherwise give a step
            # stable only at the start of the run.
            p_max = pressure
            if time_functions !== nothing && t_end > 0.0
                tf_id = get(mesh.absolute_pressure_tf, node_id, 0)
                if tf_id != 0
                    p_max = maximum_over_stage(time_functions, tf_id, pressure, t_end)
                end
            end
            C_max = max(C_max, p_max / (R * T_conv))
        end
    end

    return C_max
end


"""
    get_minimum_reference_temperature(mesh) -> Float64

Return the lowest mesh-file temperature over the initial field and prescribed temperature
boundaries; 293.15 K if none.

Converts prescribed pressures to concentrations for stability bounds: lowest temperature
gives the largest concentration and the smallest step.
"""
function get_minimum_reference_temperature(mesh)
    T_min = Inf

    for (_, T) in mesh.initial_temperature
        T > 0.0 && (T_min = min(T_min, T))
    end
    for (_, T) in mesh.temperature_bc
        T > 0.0 && (T_min = min(T_min, T))
    end

    return isfinite(T_min) ? T_min : 293.15
end


"""
    maximum_over_stage(time_functions, tf_id, base_value, t_end; nsamples) -> Float64

Return the largest value of a time function over `[0, t_end]` by dense sampling.

Sampling works for every `TimeFunction` kind (ramp, table, CSV) without duplicating the
`evaluate` dispatch. Spikes narrower than the sample spacing (a thousandth of the stage by
default) may be missed; `warn_if_aliased` flags tables sampled finer than the solver step.
"""
function maximum_over_stage(time_functions, tf_id::Int, base_value::Float64,
                            t_end::Float64; nsamples::Int = 1001)
    v_max = bc_value(time_functions, tf_id, base_value, 0.0)

    for k in 1:nsamples
        t = t_end * k / nsamples
        v_max = max(v_max, bc_value(time_functions, tf_id, base_value, t))
    end

    return v_max
end

"""
    get_maximum_initial_concentration(mesh::MeshData, NGases::Int) -> Float64

Return the maximum initial gas concentration C_g^i [mol/m³] over all elements and gases.
"""
function get_maximum_initial_concentration(mesh, NGases::Int)
    C_max = 0.0
    
    for (elem_id, concentrations) in mesh.initial_concentrations
        for gas_idx in 1:NGases
            C_max = max(C_max, concentrations[gas_idx])
        end
    end
    
    return C_max
end


"""
    get_maximum_co2_concentration(mesh::MeshData, materials) -> Float64

Return the maximum initial CO2 concentration C_CO2_max [mol/m³]; 0 if CO2 is not in the
gas dictionary.
"""
function get_maximum_co2_concentration(mesh, materials)
    # Find CO2 index in gas dictionary
    co2_idx = 0
    for (idx, gas_name) in enumerate(materials.gas_dictionary)
        if uppercase(gas_name) == "CO2" || uppercase(gas_name) == "CO₂"
            co2_idx = idx
            break
        end
    end
    
    # If CO2 not found, return 0
    if co2_idx == 0
        return 0.0
    end
    
    # Find maximum CO2 concentration
    C_co2_max = 0.0
    for (elem_id, concentrations) in mesh.initial_concentrations
        C_co2_max = max(C_co2_max, concentrations[co2_idx])
    end
    
    return C_co2_max
end


"""
    get_minimum_permeability_ratio(mesh::MeshData, materials, T_ref::Float64) -> Float64

Return the minimum of μ_g / (C_g^i × K × T × R) over all elements [s/m²], for the advective
time scale. R = 8.314 J/(mol·K).

# Arguments
- `T_ref`: Reference temperature [K]
- `C_max`: Largest total concentration; computed if not given
"""
function get_minimum_permeability_ratio(mesh, materials, T_ref::Float64,
                                        C_max::Float64 = NaN)
    R = 8.314  # Universal gas constant [J/(mol·K)]
    ratio_min = Inf

    # Get minimum gas viscosity
    μ_min = get_minimum_gas_viscosity(materials)

    # Largest total concentration. Passed by the caller such that closed form and measured
    # operator use the same scale (valid cross-check); recomputed for standalone calls.
    if !isfinite(C_max)
        C_max = get_maximum_total_concentration(mesh, length(materials.gas_dictionary))
    end
    
    # Loop through all elements to find minimum ratio
    for elem_id in 1:mesh.num_elements
        # Get material for this element
        material_idx = get_element_material(mesh, elem_id)
        
        if material_idx !== nothing
            # Get soil properties
            soil_name = materials.soil_dictionary[material_idx]
            soil = materials.soils[soil_name]
            
            # Calculate permeability ratio for this element
            # Using maximum concentration and minimum viscosity to get minimum ratio
            K = soil.intrinsic_permeability
            
            if K > 0.0 && C_max > 0.0
                # theta_g in the numerator: the advective flux, Eq. (advective_flux), has
                # no theta_g (unlike gas diffusion), therefore the lumped-mass theta_g remains:
                #     dC/dt = -(R K T C /(mu theta_g)) M_hat^-1 K_e C_tot
                # dt <= 2/lambda_max, lambda_max = 4 R K T C/(mu theta_g h^2).
                # Previous form (no theta_g, 4 instead of 2): 1/(2 theta_g) times the true
                # limit, 1.55x too permissive at theta_g = 0.308 (stable only for C_N <= 0.65).
                # Against 2/lambda_max of the assembled operator: 0.95 of the true limit.
                # Slightly conservative (nonlinear operator, element-maximum concentration
                # instead of Gauss-point values).
                θ_g_e = soil.porosity * (1.0 - soil.saturation)
                ratio = μ_min * θ_g_e / (2 * C_max * K * T_ref * R)
                ratio_min = min(ratio_min, ratio)
            end
        end
    end
    
    return ratio_min
end


"""
    get_maximum_reaction_parameters(mesh::MeshData, materials, C_co2_max::Float64) -> Float64

Return the maximum over the mesh of the CO2 reaction decay constant [1/s], for the
reactive time scale.

    (θ_w/θ_g) × k_T × a × K_H × R × T × (A_s - A_r)

Lime inventory at t = 0, area factor a at the largest total gas concentration.

# Arguments
- `T_ref`: Reference temperature [K]
"""
function get_maximum_reaction_parameters(mesh, materials, T_ref::Float64)
    param_max = 0.0
    M_lime = 74.093                 # Molar mass of Ca(OH)2 [g/mol]

    # Reaction kinetics, same for every material
    k_o = materials.reactants.arrhenius_factor
    E_a = materials.reactants.activation_energy
    β_area = materials.reactants.interfacial_area_beta

    # Area factor scales the rate and the decay constant, therefore enters the time step.
    # Increases with TOTAL gas pressure: largest total concentration gives the shortest step.
    # Species summed (argument of the factor). Same helper as the advective bound, including
    # prescribed inlets (an injection face above the initial field would give a step too long).
    #
    # T_ref for the pressure conversion: bound formed before the run from initial data. A
    # self-heating specimen at constant volume increases its pressure (and a), changes T and
    # consumes lime. The IMPES solver recomputes this bound each step; the returned value is
    # the starting estimate printed in the start-up report.
    C_g_max = get_maximum_total_concentration(mesh, length(materials.gas_dictionary))

    # Loop through all elements
    for elem_id in 1:mesh.num_elements
        # Get material for this element
        material_idx = get_element_material(mesh, elem_id)

        if material_idx !== nothing
            # Soil properties and reaction inputs of the same material
            soil_name = materials.soil_dictionary[material_idx]
            soil = materials.soils[soil_name]
            reaction = materials.reactions[soil_name]

            n = soil.porosity
            θ_w = n * soil.saturation        # Volumetric water content [-]
            θ_g = n * (1.0 - soil.saturation)  # Volumetric gas content [-]

            if θ_w > 0.0 && θ_g > 0.0
                # Lime above the residual at t = 0, per unit total volume. Monotonically
                # consumed, therefore t = 0 is the worst case for the whole run.
                C_lime_0 = (reaction.lime_content * soil.specific_gravity *
                            (1.0 - n) * 1e6) / M_lime          # [mol/m³ total]
                A_react = C_lime_0 * (1.0 - reaction.residual_lime)

                # Decay constant of the first-order equation. Same helper as the solver
                # adaptive step, consistent start-up and running values.
                #
                # No solubility cap: A_react = full lime above the residual instead of
                # A_aq,sat. Two orders of magnitude larger for typical mixtures, reactive
                # time scale correspondingly shorter.
                #
                # C_g_max passed as the area-factor argument (equivalent to a pre-multiplied
                # a_max at the same total concentration and T_ref).
                param = reaction_decay_constant(C_g_max, C_lime_0, C_lime_0 - A_react,
                                                θ_w, θ_g, T_ref, k_o, E_a, β_area)
                param_max = max(param_max, param)
            end
        end
    end

    return param_max
end


"""
    get_minimum_gas_volume_fraction(materials) -> Float64

Return the minimum gas volume fraction θ_g = n - θ_w [-] over all soils.
"""
function get_minimum_gas_volume_fraction(materials)
    θ_g_min = Inf
    
    for soil_name in materials.soil_dictionary
        soil = materials.soils[soil_name]
        n = soil.porosity
        θ_w = n * soil.saturation
        θ_g = n - θ_w
        
        θ_g_min = min(θ_g_min, θ_g)
    end

    return θ_g_min
end


#------------------------------------------------------------------------------
# Measured operator spectrum
#
# Forward Euler stability limit: Δt <= 2/λ_max, λ_max = largest eigenvalue of the
# assembled semi-discrete operator M̂⁻¹K. The closed forms below use 4·coef/h_min²,
# valid only for a uniform PLANE quad mesh.
#
# Axisymmetric formulation: incorrect at the axis. M̂ and K both include r_gp |det J|
# (`ShapeFunctions.get_weight`); away from the axis the radius cancels (plane value).
# In the element touching r = 0, r varies from 0 to h and does not cancel. Exact
# integrals for 1D radial linear elements:
#
#     interior node r_j:   M_j = h·r_j     K_jj = 2r_j/h    λ = 2K/M = 4/h²
#     axis node r = 0:     M_0 = h²/6      K_00 = 1/2       λ = 2K/M = 6/h²
#
# λ_max at the axis = 1.5x the plane value, closed-form step 1.5x too large. At the
# default Courant number 0.9: effective Courant 1.35, unstable on the axis nodes, bounded
# by the solver clamp into a standing sawtooth instead of visible divergence.
#
# λ_max measured directly from the assembled operator instead of special cases for
# axisymmetry, graded meshes, distorted elements or new element types. The assembly
# routines already include the correct measure.
#------------------------------------------------------------------------------

"""
    apply_element_stiffness!(y, x, mesh, K_elements, coef) -> y

Compute the matrix-free product `y = Σ_e coef[e] · K_e · x_e`, scattered to nodes. Mutates `y`.

`K_elements`: geometric stiffness matrices of `assemble_element_stiffness_matrices`
(measure includes `r_gp` for axisymmetry), valid for both formulations.

Serial on purpose: a few hundred calls at startup (milliseconds). No per-chunk buffers,
thread-count-independent summation order.
"""
function apply_element_stiffness!(y::Vector{Float64}, x::Vector{Float64}, mesh,
                                  K_elements::Vector{Matrix{Float64}}, coef::Vector{Float64})
    fill!(y, 0.0)

    @inbounds for e in 1:mesh.num_elements
        c = coef[e]
        (c == 0.0 || !isfinite(c)) && continue

        K_e = K_elements[e]
        n1 = mesh.elements[e, 1]; n2 = mesh.elements[e, 2]
        n3 = mesh.elements[e, 3]; n4 = mesh.elements[e, 4]
        x1 = x[n1]; x2 = x[n2]; x3 = x[n3]; x4 = x[n4]

        y[n1] += c * (K_e[1,1]*x1 + K_e[1,2]*x2 + K_e[1,3]*x3 + K_e[1,4]*x4)
        y[n2] += c * (K_e[2,1]*x1 + K_e[2,2]*x2 + K_e[2,3]*x3 + K_e[2,4]*x4)
        y[n3] += c * (K_e[3,1]*x1 + K_e[3,2]*x2 + K_e[3,3]*x3 + K_e[3,4]*x4)
        y[n4] += c * (K_e[4,1]*x1 + K_e[4,2]*x2 + K_e[4,3]*x3 + K_e[4,4]*x4)
    end

    return y
end


"""
    measure_max_eigenvalue(mesh, M, K_elements, coef; tol, maxiter) -> Float64

Compute the largest eigenvalue of `M̂⁻¹ Σ_e coef[e] K_e` by power iteration, i.e. the
largest `λ` of the generalised symmetric problem `K u = λ M u`.

`K` positive semi-definite, `M` positive diagonal: real non-negative eigenvalues,
convergence to `λ_max` without shift. Estimate: Rayleigh quotient `(xᵀKx)/(xᵀMx)`.

The Rayleigh quotient approaches `λ_max` FROM BELOW, therefore under-convergence gives a
step too LONG (unsafe). The top of a Laplacian-like spectrum is densely clustered and the
quotient plateaus early: on the 2 mm carbonation mesh, 50 iterations are 1.2% low and 200
are 0.12% low. 2000 iterations: within 0.006% of a dense eigensolve, ~50 ms once at
startup (default).

The caller never lets the measurement lengthen the step relative to the closed form (see
`calculate_critical_time_step`).

Deterministic starting vector, reproducible step. Smooth plus checkerboard components:
the top eigenvector is the highest-frequency mesh mode, slow to develop from a smooth start.

Returns 0.0 for an empty operator (all coefficients zero, or degenerate mass): no limit.
"""
function measure_max_eigenvalue(mesh, M::Vector{Float64}, K_elements::Vector{Matrix{Float64}},
                                coef::Vector{Float64}; tol::Float64 = 1e-6, maxiter::Int = 2000)
    Nnodes = length(M)
    Nnodes > 0 || return 0.0
    any(c -> c != 0.0 && isfinite(c), coef) || return 0.0

    x = Vector{Float64}(undef, Nnodes)
    y = Vector{Float64}(undef, Nnodes)

    # Deterministic seed. Golden angle in radians: no resonance with regular node
    # numbering. Alternating term: high-frequency content of the top eigenvector.
    @inbounds for i in 1:Nnodes
        x[i] = sin(2.399963229728653 * i) + (isodd(i) ? 0.5 : -0.5)
    end

    nrm = sqrt(sum(abs2, x))
    nrm > 0.0 || return 0.0
    x ./= nrm

    λ = 0.0
    for it in 1:maxiter
        apply_element_stiffness!(y, x, mesh, K_elements, coef)

        num = 0.0
        den = 0.0
        @inbounds for i in 1:Nnodes
            num += x[i] * y[i]
            den += M[i] * x[i] * x[i]
        end
        den > 0.0 || return 0.0
        λ_new = num / den

        # x <- M⁻¹Kx, renormalised. Massless nodes set to zero
        @inbounds for i in 1:Nnodes
            x[i] = M[i] > 0.0 ? y[i] / M[i] : 0.0
        end
        nrm = sqrt(sum(abs2, x))
        nrm > 0.0 || return max(λ_new, 0.0)
        x ./= nrm

        if it > 1 && abs(λ_new - λ) <= tol * max(abs(λ_new), eps())
            return max(λ_new, 0.0)
        end
        λ = λ_new
    end

    return max(λ, 0.0)
end


"""
    build_element_coefficients(mesh, materials, mechanism, D_max, C_max, T_ref) -> Vector{Float64}

Return the per-element scalar of the geometric stiffness for one mechanism, such that
`M̂⁻¹ Σ_e coef[e] K_e` reproduces the operator of `fully_explicit_diffusion_solver`.

- `:diffusion`  `θ_g · D_max · τ`, same scaling as the diffusive flux `θ_g * D_g * τ`.
  `θ_g` cancels with `M`, as in the closed form.
- `:advection`  `R · K · C_max · T_ref / μ_min`, same as the advective flux
  `R k C_g T / μ_g`. Largest total concentration and smallest viscosity: tightest bound
  of the frozen-coefficient linearisation.
- `:conduction` `λ_e` (effective conductivity), paired with the lumped thermal capacity
  instead of `M`.
- `:su_stabilization` `R · K · C_max · T_ref / μ_min`, identical to `:advection`. For any
  nodal vector `x`, Cauchy-Schwarz gives
  `xᵀK^SU_e x = Σ_gp τ*[v·∇C_h]² w ≤ (max τ*|v|²) xᵀK_e x`, therefore the SU Rayleigh
  quotient is bounded by an extra isotropic coefficient on the SAME geometric `K_e`.
  `τ* ≤ h_e/(2|v|)` unconditionally (Brooks & Hughes, independent of Pe), so
  `τ*|v|² ≤ h_e|v|/2`. With `v = -(k/μ_g)RT∇C_tot` and `|∇C_tot| ≤ C_max/h_e`, `h_e`
  cancels: `coef ≤ Rk C_max T_ref/(2μ_min) = coef_advection/2`. The thermal-pressure
  companion (`q_thermal_press`) gives the same bound (`T_ref` in place of `C_max`),
  therefore both SU terms together are bounded by `coef_advection`. Full derivation in
  `src/ADVECTION_STABILIZATION_NOTES.md`.
"""
function build_element_coefficients(mesh, materials, mechanism::Symbol,
                                    D_max::Float64, C_max::Float64, T_ref::Float64)
    R = 8.314
    Nelements = mesh.num_elements
    coef = zeros(Float64, Nelements)

    μ_min = get_minimum_gas_viscosity(materials)

    # Gas contribution to the effective conductivity, :conduction only (same as
    # get_minimum_conduction_time_scale)
    λ_g = 0.0
    NGases = length(materials.gas_dictionary)
    if mechanism === :conduction && NGases > 0
        for gas_name in materials.gas_dictionary
            λ_g += materials.gases[gas_name].thermal_conductivity
        end
        λ_g /= NGases
    end

    for e in 1:Nelements
        material_idx = get_element_material(mesh, e)
        material_idx === nothing && continue

        soil = materials.soils[materials.soil_dictionary[material_idx]]
        n = soil.porosity
        θ_w = n * soil.saturation
        θ_g = n * (1.0 - soil.saturation)

        if mechanism === :diffusion
            coef[e] = θ_g * D_max * soil.granular_tortuosity
        elseif mechanism === :advection || mechanism === :su_stabilization
            # su_stabilization: same formula (coef_su <= coef_advection, see docstring)
            K = soil.intrinsic_permeability
            coef[e] = (K > 0.0 && μ_min > 0.0) ? R * K * C_max * T_ref / μ_min : 0.0
        elseif mechanism === :conduction
            coef[e] = (1.0 - n) * soil.thermal_conductivity_solids +
                      θ_w * materials.liquid.thermal_conductivity + θ_g * λ_g
        else
            error("Unknown mechanism $(mechanism) in build_element_coefficients")
        end
    end

    return coef
end


"""
    assemble_thermal_capacity_vector!(C_T, mesh, materials) -> C_T

Assemble the lumped volumetric heat capacity per node, `C_a = Σ_gp C_mix N_a w_gp W_gp`.
Mutates `C_T`.

Energy-equation counterpart of `assemble_lumped_mass_vector!`: conduction operator
`C_T⁻¹ Σ_e λ_e K_e`. `C_mix` formed as in `get_minimum_conduction_time_scale`.
"""
function assemble_thermal_capacity_vector!(C_T::Vector{Float64}, mesh, materials)
    fill!(C_T, 0.0)

    NGases = length(materials.gas_dictionary)
    C_g_max = get_maximum_initial_concentration(mesh, NGases)

    ρc_g = 0.0
    for gas_name in materials.gas_dictionary
        ρc_g += C_g_max * materials.gases[gas_name].molar_heat_capacity
    end

    for e in 1:mesh.num_elements
        material_idx = get_element_material(mesh, e)
        material_idx === nothing && continue

        soil = materials.soils[materials.soil_dictionary[material_idx]]
        n = soil.porosity
        θ_w = n * soil.saturation
        θ_g = n * (1.0 - soil.saturation)
        ρ_s = soil.specific_gravity * materials.liquid.density

        C_mix = (1.0 - n) * ρ_s * soil.specific_heat_solids +
                θ_w * materials.liquid.density * materials.liquid.specific_heat +
                θ_g * ρc_g

        for p in 1:4
            N_p = ShapeFunctions.get_N(p)
            dV = ShapeFunctions.get_weight(e, p) * ShapeFunctions.shape_funcs.gauss_weights[p]
            for i in 1:4
                C_T[mesh.elements[e, i]] += C_mix * N_p[i] * dV
            end
        end
    end

    return C_T
end


"""
    measure_operator_eigenvalues(mesh, materials, D_max, C_max, T_ref, active) -> Dict{Symbol,Float64}

Return the measured `λ_max` per active mechanism; empty `Dict` if the measurement is not
possible.

Not possible only with an uninitialized shape function cache (tests and tooling without
`initialize_shape_functions!`). The caller then uses the closed forms WITHOUT the
axisymmetric correction; therefore `kernel.jl` initializes the cache at step [6/8],
before the time step at [7/8].
"""
function measure_operator_eigenvalues(mesh, materials, D_max::Float64, C_max::Float64,
                                      T_ref::Float64, active::Dict{Symbol,Bool})
    λ = Dict{Symbol,Float64}()

    ShapeFunctions.shape_funcs === nothing && return λ

    M = zeros(Float64, mesh.num_nodes)
    assemble_lumped_mass_vector!(M, mesh, materials)
    K_elements = assemble_element_stiffness_matrices(mesh)

    for mech in (:diffusion, :advection, :su_stabilization)
        get(active, mech, false) || continue
        coef = build_element_coefficients(mesh, materials, mech, D_max, C_max, T_ref)
        λ[mech] = measure_max_eigenvalue(mesh, M, K_elements, coef)
    end

    if get(active, :conduction, false)
        C_T = zeros(Float64, mesh.num_nodes)
        assemble_thermal_capacity_vector!(C_T, mesh, materials)
        coef = build_element_coefficients(mesh, materials, :conduction, D_max, C_max, T_ref)
        λ[:conduction] = measure_max_eigenvalue(mesh, C_T, K_elements, coef)
    end

    return λ
end


"""
    calculate_critical_time_step(mesh::MeshData, materials, T_ref::Float64) -> Float64

Compute the critical time step Δt_crit [s] from the stability criteria:
1. Diffusive time scale: h_min² / (D_max × τ)
2. Advective time scale: h_min² × (μ_g / (C_g^i × K × T × R))_min
3. Reactive time scale: 1 / ((θ_w/θ_g) × k_T × a × K_H × R × T × (A_s - A_r))

# Arguments
- `T_ref`: Reference temperature [K]

# Formula
```
Δt_crit = min{ h_min² / (D_max × τ),
               h_min² × (μ_g / (C_g^i × K × T × R))_min,
               1 / ((θ_w/θ_g) × k_T × a × K_H × R × T × (A_s - A_r)) }
```
"""
function calculate_critical_time_step(mesh, materials, T_ref::Float64, solver_settings = nothing;
                                      time_functions = nothing, t_end::Float64 = 0.0)
    # Find minimum characteristic length
    h_min = find_minimum_characteristic_length(mesh)

    # Active mechanisms only. Inactive terms must not bound the step (previously a
    # diffusion-only run reported "Advective" as limiting). All on without settings.
    on(key) = solver_settings === nothing || get(solver_settings, key, 1) == 1
    use_diffusion   = on("diffusion")
    use_advection   = on("advection")
    use_reaction    = on("reaction_kinetics")
    use_conduction  = on("heat_conduction")
    # SU stabilization adds a second diffusion-like operator of the same size as the
    # advective one (see build_element_coefficients, :su_stabilization), with its own term
    # in the harmonic sum. Without it dt is ~2x too long for the stabilized operator, and
    # the clamp remains active with advection_stabilization on.
    use_su_stabilization = use_advection && on("advection_stabilization")

    # IMPES: implicit advective term (see IMPES_FORMULATION_NOTES.md), therefore neither
    # advection nor SU bounds the explicit step. Both still MEASURED and reported
    # (justification of the scheme), excluded from the harmonic sum. Replaced by the
    # Courant limit θ_g h / |v|, which requires the velocity field; the IMPES solver
    # recomputes the step each step, and this value is the starting one.
    use_impes = solver_settings !== nothing && get(solver_settings, "impes", false)

    # Get maximum diffusion coefficient
    D_max = get_maximum_diffusion_coefficient(materials)

    # Maximum tortuosity over all soils. D_eff = θ_g τ D, Eq. (effective_diffusion): τ is
    # a multiplier, therefore the largest τ gives the tightest bound.
    τ_max = 0.0
    for soil_name in materials.soil_dictionary
        soil = materials.soils[soil_name]
        τ_max = max(τ_max, soil.granular_tortuosity)
    end

    # Largest total gas concentration (initial field and prescribed boundaries). Scales
    # the advective operator, used by the closed form (get_minimum_permeability_ratio)
    # and the measured one.
    C_max = get_maximum_total_concentration(mesh, length(materials.gas_dictionary);
                                            time_functions = time_functions, t_end = t_end)

    # Gas diffusion.
    #
    #     Δt <= h_min^2 / (2 D_max τ)
    #
    # No theta_g. Discrete equation:
    #     theta_g M_hat dC/dt = -(theta_g D tau) K C
    # theta_g cancels; operator (D tau) M_hat^-1 K, largest eigenvalue 4(D tau)/h^2 on a
    # uniform quad mesh. Forward Euler: dt <= 2/lambda_max.
    #
    # Earlier form h^2 tau/(4 theta_g D) had a spurious 1/theta_g, 4 instead of 2, and tau
    # in the numerator (from D_eff = theta_g D / tau). Same eigenvalue argument with D_eff
    # multiplied by tau gives the expression above.
    dt_diffusion = Inf
    if use_diffusion && D_max > 0.0
        dt_diffusion = h_min^2 / (2 * D_max * τ_max)
    end

    # Gas advection
    dt_advection = Inf
    if use_advection
        permeability_ratio_min = get_minimum_permeability_ratio(mesh, materials, T_ref, C_max)
        if permeability_ratio_min < Inf
            dt_advection = h_min^2 * permeability_ratio_min
        end
    end

    # Closed forms (uniform plane quad mesh) kept for the cross-check with the measured
    # spectrum. Disagreement above a few percent: mesh or formulation outside their
    # validity, measurement preferred.
    dt_diffusion_analytic = dt_diffusion
    dt_advection_analytic = dt_advection

    # SU stabilization, closed form. coef_su == coef_advection (see
    # build_element_coefficients), same closed-form step as advection before the
    # independent measured spectra below.
    dt_su = use_su_stabilization ? dt_advection : Inf
    dt_su_analytic = dt_su

    # Carbonation reaction
    dt_reaction = Inf
    if use_reaction
        reaction_param_max = get_maximum_reaction_parameters(mesh, materials, T_ref)
        if reaction_param_max > 0.0
            dt_reaction = 1.0 / (2 * reaction_param_max)
        end
    end

    # Thermal conduction. Same operator structure as gas diffusion with
    # (theta_g D / tau) -> lambda_e and theta_g -> C_mix, therefore factor 2 instead of 4.
    dt_conduction = Inf
    if use_conduction
        dt_conduction = get_minimum_conduction_time_scale(mesh, materials, h_min)
    end
    dt_conduction_analytic = dt_conduction

    #--------------------------------------------------------------------------
    # Measured spectrum. Each closed form replaced by 2/λ_max of the assembled operator,
    # the only correct version under axisymmetry (see block comment above
    # `apply_element_stiffness!`).
    #
    # Measured λ_max = 0: no operator on this mesh, no limit (Inf, no closed-form
    # fallback). Measurement not possible (uninitialized shape functions): closed forms kept.
    #
    # The measurement can only SHORTEN the step: `dt = min(measured, analytic)`.
    #
    # 1. Power iteration approaches λ_max from below; under-convergence would lengthen
    #    the step (unsafe). The minimum prevents it by construction.
    #
    # 2. Closed forms are EXACT on a uniform plane quad mesh (2 mm mesh: same λ_max as a
    #    dense eigensolve to all printed digits). Where they disagree, the mesh or
    #    formulation is outside their validity and the measurement is larger, i.e. the
    #    minimum selects it.
    #--------------------------------------------------------------------------
    active = Dict(:diffusion  => use_diffusion && D_max > 0.0,
                  :advection  => use_advection,
                  :conduction => use_conduction,
                  :su_stabilization => use_su_stabilization)

    λ_measured = measure_operator_eigenvalues(mesh, materials, D_max, C_max, T_ref, active)
    measured = !isempty(λ_measured)

    dt_from(λ) = λ > 0.0 ? 2.0 / λ : Inf

    if haskey(λ_measured, :diffusion)
        dt_diffusion = min(dt_diffusion, dt_from(λ_measured[:diffusion]))
    end
    if haskey(λ_measured, :advection)
        dt_advection = min(dt_advection, dt_from(λ_measured[:advection]))
    end
    if haskey(λ_measured, :conduction)
        dt_conduction = min(dt_conduction, dt_from(λ_measured[:conduction]))
    end
    if haskey(λ_measured, :su_stabilization)
        dt_su = min(dt_su, dt_from(λ_measured[:su_stabilization]))
    end

    # Diffusion-like operators on the same node: eigenvalues add, therefore the combined
    # limit is a harmonic sum instead of a minimum,
    #     1/dt_combined = 1/dt_a + 1/dt_b + ...
    # (single-term limit when only one is active). advection_1d: molecular and advective
    # Fourier numbers 0.315 and 0.439 (each stable), 0.754 combined, NaN.
    #
    # IMPES: advective and SU operators excluded (not integrated explicitly).
    explicit_dts = use_impes ? (dt_diffusion, dt_conduction) :
                               (dt_diffusion, dt_advection, dt_conduction, dt_su)
    inv_sum = 0.0
    for dtx in explicit_dts
        isfinite(dtx) && dtx > 0.0 && (inv_sum += 1.0 / dtx)
    end
    dt_diffusive_like = inv_sum > 0.0 ? 1.0 / inv_sum : Inf

    dt_min = min(dt_diffusive_like, dt_reaction)

    # Name of the largest single contributor, for the log. su_stabilization never named:
    # dt_su == dt_advection when active (same coefficient and spectrum), and the
    # "Advective" branch wins the tie.
    limiting_scale = "Unknown"
    if dt_min == dt_reaction
        limiting_scale = "Reactive"
    elseif use_impes
        # dt_advection and dt_su implicit, excluded from the name
        smallest = min(dt_diffusion, dt_conduction)
        limiting_scale = smallest == dt_diffusion ? "Diffusive" : "Thermal conduction"
        if inv_sum > 0.0 && dt_diffusive_like < 0.95 * smallest
            limiting_scale *= " (combined)"
        end
    else
        smallest = min(dt_diffusion, dt_advection, dt_conduction, dt_su)
        limiting_scale = smallest == dt_diffusion ? "Diffusive" :
                         smallest == dt_advection ? "Advective" : "Thermal conduction"
        if inv_sum > 0.0 && dt_diffusive_like < 0.95 * smallest
            limiting_scale *= " (combined)"
        end
    end

    # No bound on the starting step: IMPES without diffusion, conduction and reaction
    # (e.g. axisym_darcy_radial) has only the Courant limit, which requires a velocity
    # field. Start from the output interval; adaptive control from the first step.
    if use_impes && !isfinite(dt_min)
        limiting_scale = "Courant (IMPES, set at runtime)"
    end

    # Diagnostics for the log and --report-timestep. λ_analytic = 2/Δt of the closed-form step
    λ_of(dt) = (isfinite(dt) && dt > 0.0) ? 2.0 / dt : NaN

    mechanisms = [
        MechanismStability("Diffusive", use_diffusion,
                           get(λ_measured, :diffusion, NaN),
                           use_diffusion ? λ_of(dt_diffusion_analytic) : NaN,
                           dt_diffusion),
        # IMPES suffix: measured but excluded from the step. Visible on purpose to show the
        # gain over the explicit advective dt.
        MechanismStability(use_impes ? "Advective [implicit]" : "Advective", use_advection,
                           get(λ_measured, :advection, NaN),
                           use_advection ? λ_of(dt_advection_analytic) : NaN,
                           dt_advection),
        MechanismStability("Thermal conduction", use_conduction,
                           get(λ_measured, :conduction, NaN),
                           use_conduction ? λ_of(dt_conduction_analytic) : NaN,
                           dt_conduction),
        MechanismStability(use_impes ? "SU Stabilization [implicit]" : "SU Stabilization",
                           use_su_stabilization,
                           get(λ_measured, :su_stabilization, NaN),
                           use_su_stabilization ? λ_of(dt_su_analytic) : NaN,
                           dt_su),
        # Placeholder for the governing IMPES limit. No closed form or eigenvalue:
        # θ_g h / |v| requires the velocity field of the first step.
        MechanismStability("Courant (IMPES)", use_impes, NaN, NaN, Inf),
        # Reaction: nodal ODE, no spatial operator, no eigenvalue, no mesh dependence
        MechanismStability("Reactive", use_reaction, NaN, NaN, dt_reaction),
    ]

    # is_axisymmetric() throws on an uninitialized cache (`measured` false): report plane
    # instead of an error in a diagnostics field
    axisym = ShapeFunctions.shape_funcs !== nothing && ShapeFunctions.is_axisymmetric()

    report = StabilityReport(mechanisms, measured, axisym,
                             h_min, C_max, dt_min, limiting_scale)

    return dt_min, limiting_scale, report
end


"""
    get_minimum_conduction_time_scale(mesh, materials, h_min) -> Float64

Return the smallest explicit stability limit [s] of the heat conduction operator,
fourth entry of Eq. (time_step).

    Δt <= h_min^2 C_mix / (2 λ_e)

`Inf` if λ_e is zero everywhere.

Bound for conduction only. Combination with diffusion and advection handled by the
harmonic sum in `calculate_critical_time_step`.

# Arguments
- `h_min`: Smallest characteristic element size [m]
"""
function get_minimum_conduction_time_scale(mesh, materials, h_min::Float64)
    dt_min = Inf
    NGases = length(materials.gas_dictionary)

    # Gas contribution to C_mix, from the largest initial concentration present
    C_g_max = get_maximum_initial_concentration(mesh, NGases)
    ρc_g = 0.0
    λ_g = 0.0
    for gas_name in materials.gas_dictionary
        ρc_g += C_g_max * materials.gases[gas_name].molar_heat_capacity
        λ_g += materials.gases[gas_name].thermal_conductivity
    end
    NGases > 0 && (λ_g /= NGases)

    for soil_name in materials.soil_dictionary
        soil = materials.soils[soil_name]
        n = soil.porosity
        θ_w = n * soil.saturation
        θ_g = n * (1.0 - soil.saturation)

        λ_e = (1.0 - n) * soil.thermal_conductivity_solids +
              θ_w * materials.liquid.thermal_conductivity + θ_g * λ_g
        λ_e > 0.0 || continue

        ρ_s = soil.specific_gravity * materials.liquid.density
        C_mix = (1.0 - n) * ρ_s * soil.specific_heat_solids +
                θ_w * materials.liquid.density * materials.liquid.specific_heat +
                θ_g * ρc_g
        C_mix > 0.0 || continue

        dt_min = min(dt_min, h_min^2 * C_mix / (2.0 * λ_e))
    end

    return dt_min
end


"""
    calculate_time_step_info(mesh::MeshData, materials, calc_params::Dict) -> TimeStepData

Compute all time stepping information of the simulation.

# Arguments
- `calc_params`: Calculation parameters, including Courant number and total time

# Example
```julia
time_info = calculate_time_step_info(mesh, materials, calc_params)
println("Critical time step: ", time_info.critical_dt, " s")
println("Actual time step: ", time_info.actual_dt, " s")
println("Number of steps: ", time_info.num_steps)
```
"""
function calculate_time_step_info(mesh, materials, calc_params::Dict)
    time_data = TimeStepData()
    
    # Get reference temperature from initial conditions (use maximum temperature)
    T_ref = 298.15  # Default to 298.15 K (25°C) if not specified

    T_max = 0.0
    for (elem_id, temp) in mesh.initial_temperature
        T_max = max(T_max, temp)
    end
    if T_max > 0.0
        T_ref = T_max
    end

    
    # Critical time step. Time functions and stage length: ramped pressure boundaries
    # bounded over the whole stage instead of at t = 0
    time_data.critical_dt, limiting_scale, stability =
        calculate_critical_time_step(mesh, materials, T_ref,
                                     get(calc_params, "solver_settings", nothing);
                                     time_functions = get(calc_params, "time_functions", nothing),
                                     t_end = Float64(calc_params["time_stepping"]["total_simulation_time"]))
    time_data.stability = stability
    
    # Get Courant number from calculation parameters
    time_data.courant_number = calc_params["time_stepping"]["courant_number"]

    # Get time per load step (when data is saved)
    time_data.time_per_step = calc_params["time_stepping"]["time_per_step"]

    # IMPES with only the Courant limit: no finite starting step (explicit mechanisms off
    # or implicit, no velocity field yet; e.g. axisym_darcy_radial, advection only).
    # Starting step = output interval, shortened by adaptive control on the first step.
    # Not reached by the explicit integrator (would require no active mechanism).
    if !isfinite(time_data.critical_dt)
        time_data.critical_dt = Float64(time_data.time_per_step)
    end

    # Calculate actual time step (apply Courant number)
    time_data.actual_dt = time_data.critical_dt * time_data.courant_number
    
    # Calculate number of time steps per load step
    time_data.num_steps_per_load = ceil(Int, time_data.time_per_step / time_data.actual_dt)
    
    # Get total simulation time
    time_data.total_time = calc_params["time_stepping"]["total_simulation_time"]
    
    # Calculate number of load steps
    time_data.num_load_steps = ceil(Int, time_data.total_time / time_data.time_per_step)
    
    # Calculate total number of time steps for entire simulation
    time_data.num_steps = time_data.num_steps_per_load * time_data.num_load_steps
    
    # Store minimum characteristic length
    time_data.h_min = find_minimum_characteristic_length(mesh)
    
    return time_data, limiting_scale
end


"""
    log_stability_report(report::StabilityReport, log_print)

Print the per-mechanism stability table behind the chosen critical time step.

Measured `λ_max`, closed-form estimate and their ratio, to expose meshes or formulations
outside the closed-form validity. Ratio above 1: closed-form step too long by that factor
(e.g. axisymmetric formulation at the axis).
"""
function log_stability_report(report::StabilityReport, log_print)
    log_print("   Formulation: $(report.axisymmetric ? "axisymmetric" : "plane")")

    if !report.measured
        log_print("   ⚠ Shape functions were not initialized, so λ_max could not be measured.")
        log_print("     Falling back to the closed forms, which assume a uniform plane quad")
        log_print("     mesh and are NOT corrected for the axisymmetric formulation.")
    end

    fmt(v) = isfinite(v) ? @sprintf("%13.4e", v) : @sprintf("%13s", "-")

    log_print(@sprintf("   %-20s %13s %13s %8s %13s",
                       "Mechanism", "λ measured", "λ analytic", "ratio", "dt [s]"))

    for m in report.mechanisms
        if !m.active
            log_print(@sprintf("   %-20s %13s %13s %8s %13s", m.name, "-", "-", "-", "(off)"))
            continue
        end
        ratio = (isfinite(m.λ_measured) && isfinite(m.λ_analytic) && m.λ_analytic > 0.0) ?
                @sprintf("%8.2f", m.λ_measured / m.λ_analytic) : @sprintf("%8s", "-")
        log_print(@sprintf("   %-20s %s %s %s %s",
                           m.name, fmt(m.λ_measured), fmt(m.λ_analytic), ratio, fmt(m.dt)))
    end

    log_print(@sprintf("   Largest total concentration used for the advective bound: %.4g mol/m³",
                       report.C_max))

    # IMPES note: two rows measured but excluded from the step, and the step is only the
    # starting value
    if any(m -> m.name == "Courant (IMPES)" && m.active, report.mechanisms)
        log_print("   IMPES: rows marked [implicit] were measured but do NOT bound the step -")
        log_print("          the advective term is solved implicitly and is unconditionally stable.")
        log_print("          The step shown is the STARTING step. The binding limit at run time is")
        log_print("          the Courant condition θ_g h / |v|, which needs a velocity field; the")
        log_print("          solver recomputes it every step and logs it at each output.")
    end

    # Flag closed-form vs. measurement disagreement above 10%, well outside the power
    # iteration tolerance (1e-3 on the Rayleigh quotient)
    for m in report.mechanisms
        (m.active && isfinite(m.λ_measured) && isfinite(m.λ_analytic) && m.λ_analytic > 0.0) || continue
        r = m.λ_measured / m.λ_analytic
        if r > 1.10
            log_print(@sprintf("   ⚠ %s: measured λ_max is %.2fx the closed-form estimate, so the", m.name, r))
            log_print(@sprintf("     closed form alone would have set a step %.2fx too large here.", r))
        elseif r < 0.90
            log_print(@sprintf("   • %s: measured λ_max is %.2fx the closed-form estimate; the step", m.name, r))
            log_print("     is correspondingly less restrictive than the closed form suggests.")
        end
    end

    return nothing
end


# Export all public functions and types
export TimeStepData, MechanismStability, StabilityReport
export log_stability_report
export calculate_element_characteristic_length, find_minimum_characteristic_length
export get_maximum_diffusion_coefficient, get_minimum_gas_viscosity
export get_maximum_initial_concentration, get_maximum_total_concentration, get_maximum_co2_concentration
export get_minimum_reference_temperature, maximum_over_stage
export measure_max_eigenvalue, measure_operator_eigenvalues, build_element_coefficients
export calculate_critical_time_step, calculate_time_step_info
