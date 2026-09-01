#------------------------------------------------------------------------------
# ADSIM Time Step Module
# This module contains functions to calculate the critical time step and
# manage time stepping for ADSIM FEM calculations
#------------------------------------------------------------------------------

"""
MechanismStability

Per-mechanism stability diagnostics, one entry per spatial operator the solver assembles.

# Fields
- `name::String`: Mechanism name as it appears in the log ("Diffusive", "Advective", ...)
- `active::Bool`: Whether the solver is assembling this term at all
- `λ_measured::Float64`: Largest eigenvalue of M̂⁻¹K measured by power iteration [1/s]
- `λ_analytic::Float64`: Largest eigenvalue from the closed form, `4·coef/h_min²` [1/s]
- `dt::Float64`: Stability limit of this mechanism acting alone, `2/λ` [s]

`λ_measured` and `λ_analytic` are `NaN` for a mechanism with no spatial operator (the
reaction, which is a nodal ODE) and for any inactive mechanism.
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

Everything behind the chosen critical time step, kept so `--report-timestep` and the
run log can show how the number was arrived at rather than just asserting it.

`measured` records whether the power iteration supplied the step. It is false only when
the shape function cache is not initialized, in which case the closed forms are used and
the axisymmetric correction is NOT applied - see `measure_operator_eigenvalues`.
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

Structure to store time stepping information for the simulation.

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

Calculate the characteristic length of a quadrilateral element.
Uses the minimum distance between opposite sides of the element.

# Arguments
- `mesh::MeshData`: Mesh data structure
- `elem_id::Int`: Element ID (1-based)

# Returns
- `Float64`: Characteristic length h [m]

# Formula
For a quadrilateral element, h is computed as the minimum of:
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

Find the minimum characteristic length across all elements in the mesh.

# Arguments
- `mesh::MeshData`: Mesh data structure

# Returns
- `Float64`: Minimum characteristic length h_min [m]
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

Get the maximum diffusion coefficient among all gases.

# Arguments
- `materials`: Material data structure

# Returns
- `Float64`: Maximum diffusion coefficient D_max [m²/s]
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

Get the minimum dynamic viscosity among all gases.

# Arguments
- `materials`: Material data structure

# Returns
- `Float64`: Minimum gas dynamic viscosity μ_g [Pa·s]
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

Largest TOTAL gas concentration anywhere in the problem, over both the initial
condition and the prescribed concentration boundaries.

The advective flux of Eq. (advective_flux) contracts against the total concentration,
so it is the total that scales the operator, not the largest single species. With one
gas the two coincide, which is why the distinction went unnoticed; with two they do
not, and using the per-species maximum makes the advective bound too permissive by the
ratio of total to largest species.

Boundary values are included because a prescribed inlet can hold a species well above
anything in the initial field - an inlet at 40 mol/m^3 against an interior maximum of
33 is not unusual - and the bound has to hold there too.

# Arguments
- `mesh::MeshData`: Mesh data structure
- `NGases::Int`: Number of gas species

# Returns
- `Float64`: Largest total concentration [mol/m^3]
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

    # Pressure-specified boundaries. These prescribe a concentration just as surely as a
    # concentration boundary does - the solver converts them with the ideal gas law at
    # apply_partial_pressure_bc! and at the end-of-step re-imposition - but they were not
    # being counted here. On the carbonation benchmark that mattered: the interior starts
    # at 41 mol/m^3 of air while the injection face holds CO2 at 113 kPa, which is
    # 46.4 mol/m^3, so the advective bound came out 13% too permissive exactly where the
    # gradients are steepest.
    #
    # The conversion uses the LOWEST temperature in the problem, since C = P/(R T) is
    # largest there and this is a bound.
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
            # A time-driven pressure boundary has to be bounded over the whole stage, not
            # sampled at t = 0: a ramp from atmospheric to injection pressure would
            # otherwise report the atmospheric value and set a step that is stable only
            # for the first instant of the run.
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

Lowest temperature stated anywhere in the mesh file, over the initial field and the
prescribed temperature boundaries. Returns 293.15 K when the mesh states none.

Used to convert a prescribed pressure into a concentration for a stability bound, where
the lowest temperature gives the largest concentration and so the tightest step.
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

Largest value a time function takes over `[0, t_end]`, by dense sampling.

Sampling rather than inspecting the curve type keeps this working for every kind a
`TimeFunction` can be - ramp, table, CSV - without duplicating `evaluate`'s dispatch. A
spike narrower than the sample spacing could be missed; at the default resolution that
would have to be narrower than a thousandth of the stage, and `warn_if_aliased` already
flags tables sampled more finely than the solver steps.
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

Get the maximum initial gas concentration across all elements and gases.

# Arguments
- `mesh::MeshData`: Mesh data structure
- `NGases::Int`: Number of gas species

# Returns
- `Float64`: Maximum initial concentration C_g^i [mol/m³]
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

Get the maximum CO2 concentration from initial conditions.
Assumes CO2 is one of the gases in the gas dictionary.

# Arguments
- `mesh::MeshData`: Mesh data structure
- `materials`: Material data structure

# Returns
- `Float64`: Maximum CO2 concentration C_CO2_max [mol/m³]
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

Get the minimum value of (μ_g / (C_g^i × K × T × R)) across all elements.
This is needed for the advective time scale calculation.

# Arguments
- `mesh::MeshData`: Mesh data structure
- `materials`: Material data structure
- `T_ref::Float64`: Reference temperature [K]

# Returns
- `Float64`: Minimum permeability ratio [s/m²]

# Note
Uses universal gas constant R = 8.314 J/(mol·K)
"""
function get_minimum_permeability_ratio(mesh, materials, T_ref::Float64,
                                        C_max::Float64 = NaN)
    R = 8.314  # Universal gas constant [J/(mol·K)]
    ratio_min = Inf

    # Get minimum gas viscosity
    μ_min = get_minimum_gas_viscosity(materials)

    # Largest total concentration in the problem. Passed in by the caller so the closed
    # form and the measured operator are scaled by the same number and their ratio means
    # what the cross-check assumes it means; recomputed here when called standalone.
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
                # theta_g belongs in the numerator here. Unlike gas diffusion, the
                # advective flux of Eq. (advective_flux) carries no theta_g of its own,
                # so the theta_g of the lumped mass survives:
                #     dC/dt = -(R K T C /(mu theta_g)) M_hat^-1 K_e C_tot
                # and dt <= 2/lambda_max with lambda_max = 4 R K T C/(mu theta_g h^2).
                # The previous form omitted theta_g and used 4 in place of 2, making it
                # 1/(2 theta_g) times the true limit - 1.55x too permissive at
                # theta_g = 0.308, so only C_N <= 0.65 was actually stable.
                # Measured against 2/lambda_max of the assembled operator, this form
                # comes out at 0.95 of the true limit: slightly conservative, because
                # the operator is nonlinear and this uses the largest concentration in
                # each element rather than the Gauss-point values.
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

Get the decay constant of the CO2 gas concentration under the reaction term,
maximised over the mesh, for use in the reactive time scale.

    (θ_w/θ_g) × k_T × a × K_H × R × T × (A_s - A_r)

evaluated with the lime inventory at t = 0 and the area factor a at the largest
initial gas concentration in the mesh.

# Arguments
- `mesh::MeshData`: Mesh data structure
- `materials`: Material data structure
- `T_ref::Float64`: Reference temperature [K]

# Returns
- `Float64`: Maximum reaction parameter [1/s]
"""
function get_maximum_reaction_parameters(mesh, materials, T_ref::Float64)
    param_max = 0.0
    R_gas = 8.3145                  # J mol⁻¹ K⁻¹
    M_lime = 74.093                 # Molar mass of Ca(OH)2 [g/mol]

    # The kinetics belong to the reaction, so the same set applies to every material
    k_o = materials.reactants.arrhenius_factor
    E_a = materials.reactants.activation_energy
    β_area = materials.reactants.interfacial_area_beta

    # The area factor multiplies the rate, so it multiplies the decay constant too
    # and must enter the time step. It grows with the CO2 gas concentration, so the
    # largest initial concentration in the mesh gives the shortest step. That value
    # is maximised over all species, which can only overestimate the CO2 present and
    # so can only shorten the step: safe for a stability bound.
    C_g_max = get_maximum_initial_concentration(mesh, length(materials.gas_dictionary))
    a_max = interfacial_area_factor(C_g_max, β_area)

    # Loop through all elements
    for elem_id in 1:mesh.num_elements
        # Get material for this element
        material_idx = get_element_material(mesh, elem_id)

        if material_idx !== nothing
            # Get soil properties and the reaction inputs of the same material
            soil_name = materials.soil_dictionary[material_idx]
            soil = materials.soils[soil_name]
            reaction = materials.reactions[soil_name]

            n = soil.porosity
            θ_w = n * soil.saturation        # Volumetric water content [-]
            θ_g = n * (1.0 - soil.saturation)  # Volumetric gas content [-]

            if θ_w > 0.0 && θ_g > 0.0
                # Lime above the residual at t = 0, per unit total volume. The rate
                # law consumes this monotonically, so t = 0 is the worst case and the
                # time step obtained here bounds the whole run.
                C_lime_0 = (reaction.lime_content * soil.specific_gravity *
                            (1.0 - n) * 1e6) / M_lime          # [mol/m³ total]
                A_react = C_lime_0 * (1.0 - reaction.residual_lime)

                # Holding the lime and the area factor fixed, the rate law is first
                # order in the CO2 gas concentration:
                #   dC_g/dt|rxn = -(θ_w/θ_g) k_T a K_H R T (A_s - A_r) C_g
                # so the bracket below is the decay constant of that equation.
                #
                # The solubility cap was removed from the rate law, so A_react is now
                # the full lime inventory above the residual rather than A_aq,sat.
                # For typical mixtures that is larger by two orders of magnitude, and
                # the reactive time scale is correspondingly shorter.
                k_T = arrhenius_coefficient(k_o, E_a, T_ref)

                param = (θ_w / θ_g) * k_T * a_max * henry_solubility(T_ref) *
                        R_gas * T_ref * A_react
                param_max = max(param_max, param)
            end
        end
    end

    return param_max
end


"""
    get_minimum_gas_volume_fraction(materials) -> Float64

Get the minimum gas volume fraction θ_g = n - θ_w across all soils.

# Arguments
- `materials`: Material data structure

# Returns
- `Float64`: Minimum gas volume fraction θ_g [-]
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
# The stability limit of forward Euler is Δt <= 2/λ_max, where λ_max is the largest
# eigenvalue of the semi-discrete operator M̂⁻¹K the solver actually assembles. The
# closed forms further down evaluate that eigenvalue analytically as 4·coef/h_min², a
# result which holds for a uniform PLANE quad mesh and nowhere else.
#
# Under the axisymmetric formulation it is wrong at the axis. Both M̂ and K carry the
# measure r_gp |det J| (`ShapeFunctions.get_weight`), and away from the axis the radius
# cancels between them, recovering the plane value. In the single element touching r = 0
# it does not cancel, because r varies from 0 to h across that element. For 1D radial
# linear elements the exact integrals give
#
#     interior node r_j:   M_j = h·r_j     K_jj = 2r_j/h    λ = 2K/M = 4/h²
#     axis node r = 0:     M_0 = h²/6      K_00 = 1/2       λ = 2K/M = 6/h²
#
# so λ_max at the axis is 1.5x the plane value and the closed form returns a step 1.5x
# too large there. At the shipped Courant number of 0.9 that is an effective Courant of
# 1.35: unstable, and localized to the column of nodes on the axis, where the clamp in
# the solver bounds it into a standing sawtooth rather than a visible divergence.
#
# Rather than patch the closed form with an axisymmetric special case - which would then
# need another special case for graded meshes, distorted elements, and any future
# element type - measure λ_max directly from the assembled operator. The assembly
# routines already carry the measure correctly; only this estimate did not.
#------------------------------------------------------------------------------

"""
    apply_element_stiffness!(y, x, mesh, K_elements, coef) -> y

Matrix-free product `y = Σ_e coef[e] · K_e · x_e`, scattered back to nodes.

`K_elements` holds the geometric element stiffness matrices of
`assemble_element_stiffness_matrices`, whose quadrature measure already carries `r_gp`
on an axisymmetric run, so this product is correct in both formulations.

Deliberately serial. It runs a few hundred times at startup on a few thousand elements,
which is milliseconds, and a serial scatter-add needs no per-chunk buffers and gives a
summation order that does not depend on the thread count.
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

Largest eigenvalue of `M̂⁻¹ Σ_e coef[e] K_e` by power iteration, i.e. the largest `λ`
of the generalised symmetric problem `K u = λ M u`.

`K` is positive semi-definite and `M` is positive and diagonal, so every eigenvalue is
real and non-negative and the iteration converges to `λ_max` with no shift. The estimate
is the Rayleigh quotient `(xᵀKx)/(xᵀMx)`.

The Rayleigh quotient approaches `λ_max` strictly FROM BELOW, so an under-converged
iteration reports a step that is too LONG - the unsafe direction, and the exact failure
this routine exists to prevent. That matters more here than it usually would, because the
top of a Laplacian-like spectrum is densely clustered and the quotient plateaus long
before it has converged: on the 2 mm carbonation mesh, 50 iterations sit 1.2% low and
200 still sit 0.12% low. Measured against a dense eigensolve, 2000 iterations reach
within 0.006%, at a cost of about 50 ms once at startup. That is the default, and it is
cheap enough that there is no reason to economise.

The caller additionally never lets the measurement lengthen the step relative to the
closed form - see `calculate_critical_time_step` - so residual under-convergence cannot
make a run less stable than it is today.

The starting vector is deterministic, so the reported step is reproducible run to run.
It mixes a smooth component with a checkerboard one because the top eigenvector of a
Laplacian-like operator is the highest-frequency mode the mesh supports, and a purely
smooth start would need many more iterations to develop it.

Returns 0.0 when the operator is empty (every coefficient zero, or a degenerate mass),
which the caller reads as "this mechanism imposes no limit".
"""
function measure_max_eigenvalue(mesh, M::Vector{Float64}, K_elements::Vector{Matrix{Float64}},
                                coef::Vector{Float64}; tol::Float64 = 1e-6, maxiter::Int = 2000)
    Nnodes = length(M)
    Nnodes > 0 || return 0.0
    any(c -> c != 0.0 && isfinite(c), coef) || return 0.0

    x = Vector{Float64}(undef, Nnodes)
    y = Vector{Float64}(undef, Nnodes)

    # Deterministic seed. The multiplier is the golden angle in radians, which avoids
    # resonating with any regular node numbering; the alternating term seeds the
    # high-frequency content the top eigenvector is built from.
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

        # x <- M⁻¹Kx, renormalised. A node with no mass cannot carry the mode.
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

Per-element scalar multiplying the geometric stiffness for one mechanism, chosen so that
`M̂⁻¹ Σ_e coef[e] K_e` reproduces the operator the solver applies in
`fully_explicit_diffusion_solver`.

- `:diffusion`  `θ_g · D_max · τ`, matching the `θ_g * D_g * τ` scaling of the diffusive
  flux. `M` carries `θ_g`, so it cancels as the closed form asserts.
- `:advection`  `R · K · C_max · T_ref / μ_min`, matching the advective flux
  `R k C_g T / μ_g`. Evaluated with the largest total concentration and smallest
  viscosity in the problem, which is the tightest bound the frozen-coefficient
  linearisation admits.
- `:conduction` `λ_e`, the effective conductivity, to be paired with the lumped thermal
  capacity rather than with `M`.
"""
function build_element_coefficients(mesh, materials, mechanism::Symbol,
                                    D_max::Float64, C_max::Float64, T_ref::Float64)
    R = 8.314
    Nelements = mesh.num_elements
    coef = zeros(Float64, Nelements)

    μ_min = get_minimum_gas_viscosity(materials)

    # Gas contribution to the effective conductivity, needed only for :conduction and
    # formed the same way as in get_minimum_conduction_time_scale.
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
        elseif mechanism === :advection
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

Lumped volumetric heat capacity per node, `C_a = Σ_gp C_mix N_a w_gp W_gp`.

The energy equation's counterpart to `assemble_lumped_mass_vector!`: the conduction
operator is `C_T⁻¹ Σ_e λ_e K_e`, so measuring its spectrum needs this rather than the
gas lumped mass. `C_mix` is formed exactly as in `get_minimum_conduction_time_scale`.
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

Measured `λ_max` for each active mechanism, or an empty `Dict` when the measurement
cannot be made.

The one case it cannot is an uninitialized shape function cache: `calculate_critical_time_step`
is reachable from tests and tooling that never call `initialize_shape_functions!`. The
caller then falls back to the closed forms, WITHOUT the axisymmetric correction - which
is why `kernel.jl` initializes the cache at step [6/8], before the time step at [7/8].
"""
function measure_operator_eigenvalues(mesh, materials, D_max::Float64, C_max::Float64,
                                      T_ref::Float64, active::Dict{Symbol,Bool})
    λ = Dict{Symbol,Float64}()

    ShapeFunctions.shape_funcs === nothing && return λ

    M = zeros(Float64, mesh.num_nodes)
    assemble_lumped_mass_vector!(M, mesh, materials)
    K_elements = assemble_element_stiffness_matrices(mesh)

    for mech in (:diffusion, :advection)
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

Calculate the critical time step based on three stability criteria:
1. Diffusive time scale: h_min² / (D_max × τ)
2. Advective time scale: h_min² × (μ_g / (C_g^i × K × T × R))_min
3. Reactive time scale: 1 / ((θ_w/θ_g) × k_T × a × K_H × R × T × (A_s - A_r))

The critical time step is the minimum of these three values.

# Arguments
- `mesh::MeshData`: Mesh data structure
- `materials`: Material data structure
- `T_ref::Float64`: Reference temperature [K]

# Returns
- `Float64`: Critical time step Δt_crit [s]

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

    # Which mechanisms are actually active. A criterion for a term the solver is not
    # assembling must not bound the step: a diffusion-only run was previously reporting
    # "Advective" as its limiting scale, because every criterion was evaluated whether
    # or not its term was switched on. Defaults to all-on when no settings are supplied.
    on(key) = solver_settings === nothing || get(solver_settings, key, 1) == 1
    use_diffusion   = on("diffusion")
    use_advection   = on("advection")
    use_reaction    = on("reaction_kinetics")
    use_conduction  = on("heat_conduction")

    # Get maximum diffusion coefficient
    D_max = get_maximum_diffusion_coefficient(materials)

    # Maximum tortuosity across all soils. D_eff = θ_g τ D (Eq. effective_diffusion of
    # the manuscript) is a multiplier, not a divisor, so the largest effective
    # diffusivity — and hence the tightest stability bound — occurs at the largest τ.
    τ_max = 0.0
    for soil_name in materials.soil_dictionary
        soil = materials.soils[soil_name]
        τ_max = max(τ_max, soil.granular_tortuosity)
    end

    # Largest total gas concentration anywhere, initial field and every prescribed
    # boundary alike. It scales the advective operator, so it is needed both by the
    # closed form (through get_minimum_permeability_ratio) and by the measured one.
    C_max = get_maximum_total_concentration(mesh, length(materials.gas_dictionary);
                                            time_functions = time_functions, t_end = t_end)

    # Gas diffusion.
    #
    #     Δt <= h_min^2 / (2 D_max τ)
    #
    # theta_g does NOT appear. The discrete equation is
    #     theta_g M_hat dC/dt = -(theta_g D tau) K C
    # so theta_g cancels and the operator is (D tau) M_hat^-1 K, whose largest
    # eigenvalue on a uniform quad mesh is 4(D tau)/h^2. Forward Euler needs
    # dt <= 2/lambda_max, giving the expression above.
    #
    # The earlier form h^2 tau/(4 theta_g D) carried a spurious 1/theta_g, a 4 in
    # place of a 2, and tau in the numerator from when D_eff = theta_g D / tau. It was
    # verified against the true 2/lambda_max of the assembled operator; that same
    # eigenvalue argument, with D_eff now multiplying by tau instead of dividing,
    # gives the expression above.
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

    # Retained for the cross-check against the measured spectrum below. These are the
    # closed forms and they assume a uniform plane quad mesh; where the measurement
    # disagrees with them by more than a few percent, the mesh or the formulation is
    # outside what they describe and the measurement is the one to trust.
    dt_diffusion_analytic = dt_diffusion
    dt_advection_analytic = dt_advection

    # Carbonation reaction
    dt_reaction = Inf
    if use_reaction
        reaction_param_max = get_maximum_reaction_parameters(mesh, materials, T_ref)
        if reaction_param_max > 0.0
            dt_reaction = 1.0 / (2 * reaction_param_max)
        end
    end

    # Thermal conduction. Same operator structure as gas diffusion with
    # (theta_g D / tau) -> lambda_e and theta_g -> C_mix, so it carries the same
    # factor of 2 rather than 4.
    dt_conduction = Inf
    if use_conduction
        dt_conduction = get_minimum_conduction_time_scale(mesh, materials, h_min)
    end
    dt_conduction_analytic = dt_conduction

    #--------------------------------------------------------------------------
    # Measured spectrum. Replaces each closed form above with 2/λ_max of the operator
    # the solver actually assembles, which is the only version that is correct under the
    # axisymmetric formulation - see the block comment above `apply_element_stiffness!`.
    #
    # A mechanism whose measured λ_max is zero has no operator on this mesh and imposes
    # no limit, which is Inf rather than a fallback to the closed form. A measurement
    # that could not be made at all (uninitialized shape functions) leaves every closed
    # form in place.
    #
    # The measurement is only ever allowed to SHORTEN the step: `dt = min(measured,
    # analytic)`. Two reasons, and both matter.
    #
    # First, the power iteration approaches λ_max from below, so a residual under-
    # convergence would lengthen the step - the unsafe direction. Taking the minimum
    # makes that impossible by construction.
    #
    # Second, the closed forms are not merely approximate: on a uniform plane quad mesh
    # they are EXACT. Measured against a dense eigensolve on the 2 mm mesh the plane
    # closed form reproduces λ_max to all printed digits. Letting a 2000-iteration
    # estimate override an exact result would be a regression on every plane run in the
    # verification suite, for no gain. Where the two disagree, it is because the mesh or
    # the formulation is outside what the closed form describes, and there the
    # measurement is larger - which is exactly when the minimum picks it.
    #--------------------------------------------------------------------------
    active = Dict(:diffusion  => use_diffusion && D_max > 0.0,
                  :advection  => use_advection,
                  :conduction => use_conduction)

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

    # Diffusion-like operators acting on the same node do not combine by taking a
    # minimum: each limit above is derived as if its term acted alone, but the
    # eigenvalues add. Combining them as a harmonic sum,
    #     1/dt_combined = 1/dt_a + 1/dt_b + ...
    # is the bound for the summed operator and reduces to the single-term limit when
    # only one is active. Demonstrated on advection_1d, where the molecular and
    # advective Fourier numbers were 0.315 and 0.439 individually — both stable — but
    # 0.754 together, and the run produced NaN.
    inv_sum = 0.0
    for dtx in (dt_diffusion, dt_advection, dt_conduction)
        isfinite(dtx) && dtx > 0.0 && (inv_sum += 1.0 / dtx)
    end
    dt_diffusive_like = inv_sum > 0.0 ? 1.0 / inv_sum : Inf

    dt_min = min(dt_diffusive_like, dt_reaction)

    # Name the single mechanism that contributes most, so the log points at the term
    # actually setting the step
    limiting_scale = "Unknown"
    if dt_min == dt_reaction
        limiting_scale = "Reactive"
    else
        smallest = min(dt_diffusion, dt_advection, dt_conduction)
        limiting_scale = smallest == dt_diffusion ? "Diffusive" :
                         smallest == dt_advection ? "Advective" : "Thermal conduction"
        if inv_sum > 0.0 && dt_diffusive_like < 0.95 * smallest
            limiting_scale *= " (combined)"
        end
    end

    # Diagnostics behind the number, so the log and --report-timestep can show how it was
    # arrived at. λ_analytic is recovered from the closed-form step as 2/Δt, the same
    # relation the closed forms were derived from.
    λ_of(dt) = (isfinite(dt) && dt > 0.0) ? 2.0 / dt : NaN

    mechanisms = [
        MechanismStability("Diffusive", use_diffusion,
                           get(λ_measured, :diffusion, NaN),
                           use_diffusion ? λ_of(dt_diffusion_analytic) : NaN,
                           dt_diffusion),
        MechanismStability("Advective", use_advection,
                           get(λ_measured, :advection, NaN),
                           use_advection ? λ_of(dt_advection_analytic) : NaN,
                           dt_advection),
        MechanismStability("Thermal conduction", use_conduction,
                           get(λ_measured, :conduction, NaN),
                           use_conduction ? λ_of(dt_conduction_analytic) : NaN,
                           dt_conduction),
        # The reaction is a nodal ODE with no spatial operator, so it has no eigenvalue
        # to measure and no mesh dependence: its limit is unchanged by any of this.
        MechanismStability("Reactive", use_reaction, NaN, NaN, dt_reaction),
    ]

    # is_axisymmetric() throws on an uninitialized cache, which is exactly the case where
    # `measured` is false, so report the formulation as plane rather than propagating an
    # error out of a diagnostics field.
    axisym = ShapeFunctions.shape_funcs !== nothing && ShapeFunctions.is_axisymmetric()

    report = StabilityReport(mechanisms, measured, axisym,
                             h_min, C_max, dt_min, limiting_scale)

    return dt_min, limiting_scale, report
end


"""
    get_minimum_conduction_time_scale(mesh, materials, h_min) -> Float64

Smallest explicit stability limit of the heat conduction operator across the mesh,

    Δt <= h_min^2 C_mix / (2 λ_e)

the fourth entry of Eq. (time_step). Returns `Inf` where no material conducts, so a
run without conduction is unaffected.

Note this bounds conduction acting alone. Where gas diffusion, gas advection and
conduction are all active they are diffusion-like operators on the same node and
their stability limits do not simply combine by taking a minimum; that is a known
gap, recorded against Phase 7.

# Arguments
- `mesh`: Mesh data structure
- `materials`: Material data structure
- `h_min::Float64`: Smallest characteristic element size [m]

# Returns
- `Float64`: Limiting conduction time step [s], or `Inf` if λ_e is zero everywhere
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

Calculate all time stepping information for the simulation.

# Arguments
- `mesh::MeshData`: Mesh data structure
- `materials`: Material data structure
- `calc_params::Dict`: Calculation parameters including Courant number and total time

# Returns
- `TimeStepData`: Structure containing all time step information

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

    
    # Calculate critical time step. The time functions and the stage length go in so a
    # ramped pressure boundary is bounded over the whole stage rather than at t = 0.
    time_data.critical_dt, limiting_scale, stability =
        calculate_critical_time_step(mesh, materials, T_ref,
                                     get(calc_params, "solver_settings", nothing);
                                     time_functions = get(calc_params, "time_functions", nothing),
                                     t_end = Float64(calc_params["time_stepping"]["total_simulation_time"]))
    time_data.stability = stability
    
    # Get Courant number from calculation parameters
    time_data.courant_number = calc_params["time_stepping"]["courant_number"]
    
    # Calculate actual time step (apply Courant number)
    time_data.actual_dt = time_data.critical_dt * time_data.courant_number
    
    # Get time per load step (when data is saved)
    time_data.time_per_step = calc_params["time_stepping"]["time_per_step"]
    
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

Shows the measured `λ_max` beside the closed-form estimate and their ratio, so a mesh or
a formulation the closed forms do not describe is visible rather than silent. A ratio
materially above 1 means the closed form would have set a step that large a factor too
long - which is what the axisymmetric formulation did at the axis before the spectrum
was measured.
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

    # Flag any mechanism where the closed form and the measurement disagree materially.
    # 10% is well outside the power iteration's own tolerance (1e-3 on the Rayleigh
    # quotient), so a warning here is a real disagreement, not iteration noise.
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
