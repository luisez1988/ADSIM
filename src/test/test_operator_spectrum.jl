# Regression guard for the Phase 0 stability fix.
#
# The critical time step is Δt <= 2/λ_max of the semi-discrete operator M̂⁻¹K. Before
# Phase 0 that eigenvalue was evaluated from a closed form valid only on a uniform PLANE
# quad mesh, which under-estimated it on the axisymmetric carbonation meshes and produced
# a step long enough for forward Euler to go unstable against the axis of symmetry.
#
# This checks the three properties the fix rests on:
#   1. the power iteration recovers λ_max of the assembled operator
#   2. the axisymmetric operator really does have a larger λ_max than the plane one on a
#      mesh that reaches r = 0
#   3. the step the solver will actually use never exceeds 2/λ_max -- the safety
#      invariant, and the only one whose failure produces a wrong answer rather than a
#      slow one
#
# Run:  julia -t4 src/test/test_operator_spectrum.jl

using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "..", "kernel.jl"))

const DATA = joinpath(@__DIR__, "..", "data")
const PROJECT = "2mm_size"          # axisymmetric, reaches r = 0, 2040 nodes
const T_REF = 293.15

failures = String[]
check(ok, msg) = (ok ? println("   ✓ $msg") : (println("   ✗ $msg"); push!(failures, msg)))

# λ_max of M̂⁻¹ Σ_e coef[e] K_e, computed exactly. K is PSD and M diagonal positive, so
# the generalised problem K u = λ M u maps to the symmetric standard one via M^(-1/2).
function exact_lambda_max(mesh, M, K_elements, coef)
    N = length(M)
    K = zeros(N, N)
    for e in 1:mesh.num_elements
        c = coef[e]
        c == 0.0 && continue
        for i in 1:4, j in 1:4
            K[mesh.elements[e, i], mesh.elements[e, j]] += c * K_elements[e][i, j]
        end
    end
    s = 1.0 ./ sqrt.(M)
    return maximum(eigvals(Symmetric((s .* K) .* s')))
end

function spectrum(mesh, materials, axisym::Bool, mech::Symbol, D_max, C_max)
    initialize_shape_functions!(mesh; axisymmetric = axisym)
    M = zeros(mesh.num_nodes)
    assemble_lumped_mass_vector!(M, mesh, materials)
    K_elements = assemble_element_stiffness_matrices(mesh)
    coef = build_element_coefficients(mesh, materials, mech, D_max, C_max, T_REF)
    return exact_lambda_max(mesh, M, K_elements, coef),
           measure_max_eigenvalue(mesh, M, K_elements, coef)
end

println("="^70)
println("Operator spectrum / stability limit")
println("="^70)

mesh = read_mesh_file(joinpath(DATA, "$(PROJECT).mesh"))
materials = read_materials_file(joinpath(DATA, "$(PROJECT)_mat.toml"))
D_max = get_maximum_diffusion_coefficient(materials)
C_max = get_maximum_total_concentration(mesh, length(materials.gas_dictionary))

println("\nMesh $(PROJECT): $(mesh.num_nodes) nodes, $(mesh.num_elements) elements")

# --- 1. the power iteration finds λ_max, and never overshoots it -------------------
println("\n[1] Power iteration against a dense eigensolve")
for mech in (:diffusion, :advection)
    for axisym in (true, false)
        exact, measured = spectrum(mesh, materials, axisym, mech, D_max, C_max)
        rel = measured / exact
        @printf("   %-10s %-6s exact %.6e   measured %.6e   ratio %.5f\n",
                mech, axisym ? "axi" : "plane", exact, measured, rel)
        check(0.995 <= rel <= 1.0 + 1e-9,
              "$(mech)/$(axisym ? "axi" : "plane"): measured within 0.5% of exact, and not above it")
    end
end

# --- 2. the axisymmetric penalty is real ------------------------------------------
println("\n[2] Axisymmetric operator is stiffer than the plane one on this mesh")
for mech in (:diffusion, :advection)
    e_axi, _ = spectrum(mesh, materials, true, mech, D_max, C_max)
    e_pln, _ = spectrum(mesh, materials, false, mech, D_max, C_max)
    @printf("   %-10s axi %.6e / plane %.6e = %.4f\n", mech, e_axi, e_pln, e_axi / e_pln)
    check(e_axi > e_pln * 1.01,
          "$(mech): axisymmetric λ_max exceeds plane by more than 1% (the axis mode)")
end

# --- 3. the safety invariant ------------------------------------------------------
# The step the solver will use must satisfy dt <= 2/λ_max for every active mechanism.
# This is what the whole phase exists to guarantee; it is checked through the real entry
# point, not through the internals, so a regression anywhere in that path trips it.
println("\n[3] Critical time step respects 2/λ_max through the real code path")
calc_params = get_all_calc_params(joinpath(DATA, "$(PROJECT)_calc.toml"), DATA)
axisym = calc_params["solver_settings"]["axisymmetric"]
zero_variables!(mesh, materials; axisymmetric = axisym)
apply_all_initial_conditions!(mesh, materials, calc_params["time_functions"], 0.0)
initialize_shape_functions!(mesh; axisymmetric = axisym)

time_data, limiting = calculate_time_step_info(mesh, materials, calc_params)
@printf("   critical dt = %.6e s   actual dt = %.6e s   limiting = %s\n",
        time_data.critical_dt, time_data.actual_dt, limiting)

for mech in (:diffusion, :advection)
    exact, _ = spectrum(mesh, materials, axisym, mech, D_max, C_max)
    dt_limit = 2.0 / exact
    @printf("   %-10s 2/λ_exact = %.6e s\n", mech, dt_limit)
    check(time_data.actual_dt <= dt_limit,
          "$(mech): actual dt is within the true stability limit")
end

check(time_data.stability !== nothing, "stability report was populated")
check(time_data.stability.axisymmetric == axisym, "stability report records the formulation")
check(time_data.stability.C_max >= 46.0,
      "C_max includes the pressure boundaries (>= 46 mol/m³, not the initial 41)")

println("\n" * "="^70)
if isempty(failures)
    println("All operator spectrum checks passed")
else
    println("FAILED $(length(failures)) check(s):")
    foreach(f -> println("  - $f"), failures)
end
println("="^70)
exit(isempty(failures) ? 0 : 1)
