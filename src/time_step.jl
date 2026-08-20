#------------------------------------------------------------------------------
# ADSIM Time Step Module
# This module contains functions to calculate the critical time step and
# manage time stepping for ADSIM FEM calculations
#------------------------------------------------------------------------------

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
    
    function TimeStepData()
        new(0.0, 0.0, 0.0, 0.0, 0, 0, 0, 0.0, 0.0)
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
function get_maximum_total_concentration(mesh, NGases::Int)
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

    return C_max
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
function get_minimum_permeability_ratio(mesh, materials, T_ref::Float64)
    R = 8.314  # Universal gas constant [J/(mol·K)]
    ratio_min = Inf
    
    # Get minimum gas viscosity
    μ_min = get_minimum_gas_viscosity(materials)
    
    # Get maximum initial concentration
    C_max = get_maximum_total_concentration(mesh, length(materials.gas_dictionary))
    
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


"""
    calculate_critical_time_step(mesh::MeshData, materials, T_ref::Float64) -> Float64

Calculate the critical time step based on three stability criteria:
1. Diffusive time scale: h_min² × τ / (θ_g × D_max)
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
Δt_crit = min{ h_min² × τ / (θ_g × D_max),
               h_min² × (μ_g / (C_g^i × K × T × R))_min,
               1 / ((θ_w/θ_g) × k_T × a × K_H × R × T × (A_s - A_r)) }
```
"""
function calculate_critical_time_step(mesh, materials, T_ref::Float64, solver_settings = nothing)
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

    # Minimum tortuosity across all soils
    τ_min = Inf
    for soil_name in materials.soil_dictionary
        soil = materials.soils[soil_name]
        τ_min = min(τ_min, soil.granular_tortuosity)
    end

    # Gas diffusion.
    #
    #     Δt <= h_min^2 τ / (2 D_max)
    #
    # theta_g does NOT appear. The discrete equation is
    #     theta_g M_hat dC/dt = -(theta_g D / tau) K C
    # so theta_g cancels and the operator is (D/tau) M_hat^-1 K, whose largest
    # eigenvalue on a uniform quad mesh is 4(D/tau)/h^2. Forward Euler needs
    # dt <= 2/lambda_max, giving the expression above.
    #
    # The earlier form h^2 tau/(4 theta_g D) carried a spurious 1/theta_g and a 4 in
    # place of a 2. Verified against the true 2/lambda_max of the assembled operator:
    # this expression matches it to 1.000000 at theta_g = 1.0 and at theta_g = 0.3318,
    # while the old one was 0.5x and 1.51x of it respectively. That is why a Courant
    # number above about 0.4 used to go unstable — at theta_g = 0.2 the old criterion
    # was 2.5x too permissive, so only C_N <= 0.4 was actually stable.
    dt_diffusion = Inf
    if use_diffusion && D_max > 0.0
        dt_diffusion = (h_min^2 * τ_min) / (2 * D_max)
    end

    # Gas advection
    dt_advection = Inf
    if use_advection
        permeability_ratio_min = get_minimum_permeability_ratio(mesh, materials, T_ref)
        if permeability_ratio_min < Inf
            dt_advection = h_min^2 * permeability_ratio_min
        end
    end

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

    return dt_min, limiting_scale
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

    
    # Calculate critical time step
    time_data.critical_dt, limiting_scale = calculate_critical_time_step(mesh, materials, T_ref,
                                                    get(calc_params, "solver_settings", nothing))
    
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


# Export all public functions and types
export TimeStepData
export calculate_element_characteristic_length, find_minimum_characteristic_length
export get_maximum_diffusion_coefficient, get_minimum_gas_viscosity
export get_maximum_initial_concentration, get_maximum_total_concentration, get_maximum_co2_concentration
export calculate_critical_time_step, calculate_time_step_info
