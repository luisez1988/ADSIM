#------------------------------------------------------------------------------
# ADSIM variable initialization Module
# This module contains functions to initialize variables for ADSIM
# FEM calculations
#------------------------------------------------------------------------------

using Base.Threads

#------------------------------------------------------------------------------
# Global simulation variables - shared across all modules
#------------------------------------------------------------------------------

# Dimensional parameters
global NDim::Int = 0
# Geometric formulation of the analysis. Held here beside the other dimension globals
# because it is fixed for the whole run and is needed before the shape functions - which
# also carry it - have been initialized.
global Axisymmetric::Bool = false
global Nnodes::Int = 0
global Nelements::Int = 0
global NSoils::Int = 0
global NGases::Int = 0

# State variables
global C_g::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
global P::Vector{Float64} = Float64[]
global T::Vector{Float64} = Float64[]
global v::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
global P_boundary::Matrix{Int} = Matrix{Int}(undef, 0, 0)
# Thermal Dirichlet mask, the temperature counterpart of P_boundary:
# 1 = free node, 0 = prescribed temperature. EVERY write to T must be gated by this,
# including any correction term. The concentration path was bitten by exactly that:
# the Lagrange correction was added after the P_boundary gate and walked prescribed
# concentrations off their boundary values.
global T_boundary::Vector{Int} = Int[]
# Tributary edge length of each convective-heat boundary node, for the Robin integral
global thermal_node_influences::Dict{Int, Float64} = Dict{Int, Float64}()
global λ_bc::Vector{Float64} = Float64[]

# Boundary node influence lengths
global boundary_node_influences::Dict{Int, Float64} = Dict{Int, Float64}()

# Reactive species
global C_lime::Vector{Float64} = Float64[]
global C_caco3::Vector{Float64} = Float64[]
global C_lime_residual::Vector{Float64} = Float64[]
global Caco3_max::Vector{Float64} = Float64[]

# Time derivatives
global dC_g_dt::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
global dT_dt::Vector{Float64} = Float64[]
global dC_lime_dt::Vector{Float64} = Float64[]

# Energy equation. M_T is the lumped thermal capacity of Eq. (lumped_capacity),
# rebuilt every step because C_mix depends on the current gas concentrations. The
# rest are the nodal thermal fluxes of Eq. (FEM_compact_energy):
#     M^L_T dT/dt = q_ext_T - q_cond_T - q_adv_T + q_react_T
# q_adv_T and q_ext_T are allocated here but stay zero until thermal advection and
# the thermal boundary conditions are implemented.
global M_T::Vector{Float64} = Float64[]
global q_cond_T::Vector{Float64} = Float64[]
global q_adv_T::Vector{Float64} = Float64[]
global q_react_T::Vector{Float64} = Float64[]
global q_ext_T::Vector{Float64} = Float64[]
# Total heat flux q_d + q_a at the nodes, for output only. Projected from the Gauss
# points with the normalised transfer of Eq. (nodal_velocities).
global q_heat::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)

#Analysis variables for soil carbonation
global binder_content::Vector{Float64} = Float64[]
global degree_of_carbonation::Vector{Float64} = Float64[]

#------------------------------------------------------------------------------
# Initialize all variables
#------------------------------------------------------------------------------
"""
    zero_variables!(mesh, materials)

Zero all simulation variables based on mesh and material data.
This function should be called once at the start of the simulation.
The exclamation mark indicates it modifies global variables.

# Arguments
- `mesh`: Mesh data structure containing node and element information
- `materials`: Material data structure containing soil and gas dictionaries
"""
function zero_variables!(mesh, materials; axisymmetric::Bool = false)
    global NDim, Nnodes, Nelements, NSoils, NGases, Axisymmetric
    global C_g, P, T, v, P_boundary, λ_bc, boundary_node_influences
    global T_boundary, thermal_node_influences
    global C_lime, C_caco3, C_lime_residual, binder_content, degree_of_carbonation, Caco3_max
    global dC_g_dt, dT_dt, dC_lime_dt
    global M_T, q_cond_T, q_adv_T, q_react_T, q_ext_T, q_heat
    
    # Set dimensions
    NDim = 2  # Number of spatial dimensions - TODO: generalize for 3D
    Axisymmetric = axisymmetric
    Nnodes = mesh.num_nodes
    Nelements = mesh.num_elements
    NSoils = length(materials.soil_dictionary)
    NGases = length(materials.gas_dictionary)
    
    # Allocate and initialize state variables
    C_g = zeros(Float64, Nnodes, NGases)
    P = zeros(Float64, Nnodes)
    T = zeros(Float64, Nnodes)
    v = zeros(Float64, Nnodes, NDim)
    P_boundary = ones(Int, Nnodes, NGases)  # 1 = free node, 0 = concentration BC node
    T_boundary = ones(Int, Nnodes)          # 1 = free node, 0 = prescribed temperature
    λ_bc = zeros(Float64, Nnodes)  # Lagrange multipliers for pressure BCs
    
    # Calculate and store boundary node influence lengths. Under the axisymmetric
    # formulation these are ring measures l_e * r̄_e rather than arclengths, which is what
    # the Lagrange multiplier of Eq. (ax_Lambda_solution) and the convective heat boundary
    # both require.
    boundary_influences = get_boundary_node_influences(mesh; axisymmetric = axisymmetric)
    boundary_node_influences = boundary_influences.node_influences

    # Same geometry, but over the convective-heat boundary nodes
    thermal_influences = get_boundary_node_influences(mesh, keys(mesh.convective_heat_bc);
                                                     axisymmetric = axisymmetric)
    thermal_node_influences = thermal_influences.node_influences
    
    # Allocate and initialize reactive species
    C_lime = zeros(Float64, Nnodes)
    C_caco3 = zeros(Float64, Nnodes)
    C_lime_residual = zeros(Float64, NSoils)
    
    # Allocate and initialize time derivatives
    dC_g_dt = zeros(Float64, Nnodes, NGases)
    dT_dt = zeros(Float64, Nnodes)
    dC_lime_dt = zeros(Float64, Nnodes)

    # Energy equation arrays
    M_T = zeros(Float64, Nnodes)
    q_cond_T = zeros(Float64, Nnodes)
    q_adv_T = zeros(Float64, Nnodes)
    q_react_T = zeros(Float64, Nnodes)
    q_ext_T = zeros(Float64, Nnodes)
    q_heat = zeros(Float64, Nnodes, NDim)


    # Allocate analysis variables
    binder_content = zeros(Float64, Nnodes)
    degree_of_carbonation = zeros(Float64, Nnodes)
    Caco3_max = zeros(Float64, Nnodes)

end


#------------------------------------------------------------------------------
# Apply initial conditions
#------------------------------------------------------------------------------
"""
    apply_initial_concentrations!(mesh::MeshData)

Apply initial gas concentrations from mesh data to the global C_g array.
This function reads element-based initial concentrations from the mesh and 
assigns them to all nodes within those elements.

# Arguments
- `mesh::MeshData`: Mesh data structure containing initial concentration data

# Note
- Modifies global variable `C_g`
- If a node belongs to multiple elements with different initial conditions,
  the last element's value will be used
"""
function apply_initial_concentrations!(mesh)
    global C_g, NGases
    
    # Apply element-based initial concentrations to nodes.
    # Iterate in ascending element order rather than in the Dict's own order: a node
    # shared by several elements keeps whichever value is written last, so Dict
    # iteration order would decide the initial condition at every shared node. Sorting
    # makes that rule "highest element index wins", which is stable and reproducible.
    for elem_id in sort(collect(keys(mesh.initial_concentrations)))
        concentrations = mesh.initial_concentrations[elem_id]
        # Get nodes of this element
        element_nodes = get_element_nodes(mesh, elem_id)
        
        # Apply concentrations to each node of the element
        for node_id in element_nodes
            @threads for gas_idx in 1:NGases
                C_g[node_id, gas_idx] = concentrations[gas_idx]
            end
        end
    end    
end


"""
    apply_initial_temperature!(mesh::MeshData)

Apply initial temperatures from mesh data to the global T array.
This function reads element-based initial temperatures from the mesh and 
assigns them to all nodes within those elements.

# Arguments
- `mesh::MeshData`: Mesh data structure containing initial temperature data

# Note
- Modifies global variable `T`
- If a node belongs to multiple elements with different initial conditions,
  the last element's value will be used
"""
function apply_initial_temperature!(mesh)
    global T
    
    # Apply element-based initial temperatures to nodes.
    # Ascending element order for the same reason as the concentrations above: the
    # last write to a shared node wins, so the ordering must not be left to the Dict.
    for elem_id in sort(collect(keys(mesh.initial_temperature)))
        temperature = mesh.initial_temperature[elem_id]
        # Get nodes of this element
        element_nodes = get_element_nodes(mesh, elem_id)
        
        # Apply temperature to each node of the element
        for node_id in element_nodes
            T[node_id] = temperature
        end
    end
end


"""
    apply_temperature_bc!(mesh::MeshData, time_functions = Dict{Int, TimeFunction}(), t = 0.0)

Apply prescribed nodal temperatures, Eq. (temperature_BC), and mark those nodes in
`T_boundary` so the solver never writes to them again.

Called after `apply_initial_temperature!`, so a prescribed value overrides the
element-based initial condition at the same node.

A node driven by a time function is seeded at time `t`, which is 0 for a fresh
run and the checkpoint time for a restart, so a stage 2 run picks its curve up
where stage 1 left it rather than at the origin.

# Note
- Modifies global variables `T` and `T_boundary`
"""
function apply_temperature_bc!(mesh, time_functions = Dict{Int, TimeFunction}(), t = 0.0)
    global T, T_boundary

    for (node_id, temperature) in mesh.temperature_bc
        tf_id = get(mesh.temperature_bc_tf, node_id, 0)
        T[node_id] = bc_value(time_functions, tf_id, temperature, t)
        T_boundary[node_id] = 0   # frozen from here on
    end
end


"""
    apply_concentration_bc!(mesh::MeshData, time_functions = Dict{Int, TimeFunction}(), t = 0.0)

Apply concentration boundary conditions from mesh data to the global C_g array.
This function sets fixed gas concentrations at nodes where concentration
boundary conditions are specified.

# Arguments
- `mesh::MeshData`: Mesh data structure containing concentration BC data
- `time_functions`: Curve table from the calculation file, empty for a constant model
- `t`: Time to seed time-driven nodes at, 0 for a fresh run or the checkpoint time

# Note
- Modifies global variable `C_g`
- A node with no time function keeps this value for the whole run, since
  `P_boundary` gates its rate to zero. A time-driven node is rewritten every
  step by the solver, because gating the rate holds the value still but nothing
  would otherwise advance it.
"""
function apply_concentration_bc!(mesh, time_functions = Dict{Int, TimeFunction}(), t = 0.0)
    global C_g, NGases, P_boundary

    # Apply nodal concentration boundary conditions
    for (node_id, concentrations) in mesh.concentration_bc
        tf_ids = get(mesh.concentration_bc_tf, node_id, Int[])
        @threads for gas_idx in 1:NGases
            tf_id = gas_idx <= length(tf_ids) ? tf_ids[gas_idx] : 0
            C_g[node_id, gas_idx] = bc_value(time_functions, tf_id,
                                             concentrations[gas_idx], t)
            P_boundary[node_id, gas_idx] = 0  # Mark node as having concentration BC
        end
    end
end


"""
    apply_pressure_bc!(mesh::MeshData, time_functions = Dict{Int, TimeFunction}(), t = 0.0)

Apply absolute pressure boundary conditions from mesh data to the global P array.
This function sets fixed pressures at nodes where pressure boundary conditions
are specified and restricts all gases at those nodes.

# Arguments
- `mesh::MeshData`: Mesh data structure containing pressure BC data
- `time_functions`: Curve table from the calculation file, empty for a constant model
- `t`: Time to seed time-driven nodes at, 0 for a fresh run or the checkpoint time

# Note
- Modifies global variables `P` and `P_boundary`
- P_boundary is set to 0 for all gases at pressure BC nodes
- These values should be maintained throughout the simulation for BC nodes
"""
function apply_pressure_bc!(mesh, time_functions = Dict{Int, TimeFunction}(), t = 0.0)
    global P, P_boundary

    # Apply nodal pressure boundary conditions
    for (node_id, pressure) in mesh.absolute_pressure_bc
        tf_id = get(mesh.absolute_pressure_tf, node_id, 0)
        P[node_id] = bc_value(time_functions, tf_id, pressure, t)
        # Restrict all gases at this node
        #P_boundary[node_id, :] .= 0  # Mark node as having pressure BC
    end
end


"""
    apply_initial_lime_concentration!(mesh::MeshData, materials)

Apply initial lime concentrations from material properties to the global C_lime array.
This function loops through all elements, gets their assigned material, retrieves
the lime content from that material, and assigns it to all nodes in the element.

# Arguments
- `mesh::MeshData`: Mesh data structure containing element-material assignments
- `materials`: Material data structure containing the per-material reaction inputs
  (lime content, residual lime) and the soil properties they are combined with

# Note
- Modifies global variable `C_lime`
- If a node belongs to multiple elements with different materials,
  the last element's value will be used
"""
function apply_initial_lime_concentration!(mesh, materials)
    global C_lime, C_lime_residual, Caco3_max
    
    # Loop through all elements
    for elem_id in 1:mesh.num_elements
        # Get material index for this element
        material_idx = get_element_material(mesh, elem_id)
        
        if material_idx !== nothing
            # Get the soil name from the soil dictionary
            soil_name = materials.soil_dictionary[material_idx]
            
            # Get the soil properties for this material
            soil_props = get_soil_properties(materials, soil_name)

            # Get the reaction inputs for this material. β_l and the residual
            # lime fraction live alongside the kinetics rather than in the soil
            # block, so they are looked up by the same soil name.
            reaction_props = get_reaction_properties(materials, soil_name)

            if soil_props !== nothing && reaction_props !== nothing
                # Get lime content from the reaction inputs
                β_l = reaction_props.lime_content
                G_s = soil_props.specific_gravity
                n=soil_props.porosity
                M_lime=74.093   # Molar mass of Ca(OH)2 in g/mol
                #Calculate lime concentration in mol/m^3
                lime_concentration= (β_l * G_s * (1 - n) * 1e6 ) / M_lime #Asumes ρ_w= 1000 kg/m^3

                #Calculatte reidual lime
                residual_percent= reaction_props.residual_lime
                C_lime_residual[material_idx] = residual_percent * lime_concentration

                #Calculate Caco3 max for degree of carbonation
                Caco3_max_concentration= lime_concentration  
                
                # Get nodes of this element
                element_nodes = get_element_nodes(mesh, elem_id)
                
                # Assign lime content and Caco3_max to each node of the element
                for node_id in element_nodes
                    C_lime[node_id] = lime_concentration
                    Caco3_max[node_id] = Caco3_max_concentration
                end
            end
        end
    end
end


"""
    apply_partial_pressure_bc!(mesh::MeshData, time_functions = Dict{Int, TimeFunction}(), t = 0.0)

Apply partial pressure boundary conditions from mesh data to the global C_g array.
This function sets gas concentrations at nodes where partial pressure boundary
conditions are specified, using the ideal gas law: C_g[i] = P_partial[i] / (R * T).
It also marks these nodes in P_boundary to prevent the solver from updating them.

# Arguments
- `mesh::MeshData`: Mesh data structure containing partial pressure BC data
- `time_functions`: Curve table from the calculation file, empty for a constant model
- `t`: Time to seed time-driven nodes at, 0 for a fresh run or the checkpoint time

# Note
- Modifies global variables `C_g` and `P_boundary`
- Uses ideal gas law: P_partial = C_g * R * T, where R = 8.314 J/(mol·K)
- P_boundary is set to 0 for all gases at partial pressure BC nodes
- Concentrations will be dynamically updated in solver to maintain partial pressure
"""
function apply_partial_pressure_bc!(mesh, time_functions = Dict{Int, TimeFunction}(), t = 0.0)
    global C_g, NGases, P_boundary, T

    R = 8.314  # Universal gas constant [J/(mol·K)]

    # Apply nodal partial pressure boundary conditions
    for (node_id, partial_pressures) in mesh.partial_pressure_bc
        tf_ids = get(mesh.partial_pressure_bc_tf, node_id, Int[])
        @threads for gas_idx in 1:NGases
            # Calculate concentration from partial pressure using ideal gas law
            # P_partial = C_g * R * T  =>  C_g = P_partial / (R * T)
            tf_id = gas_idx <= length(tf_ids) ? tf_ids[gas_idx] : 0
            p_partial = bc_value(time_functions, tf_id, partial_pressures[gas_idx], t)
            C_g[node_id, gas_idx] = p_partial / (R * T[node_id])
            P_boundary[node_id, gas_idx] = 0  # Mark node as having partial pressure BC
        end
    end
end


"""
    apply_all_initial_conditions!(mesh::MeshData, materials)

Apply all initial conditions and boundary conditions from mesh and material data.
This is a convenience function that calls all individual application functions.

# Arguments
- `mesh::MeshData`: Mesh data structure containing all initial and boundary condition data
- `materials`: Material data structure containing soil and gas properties
- `time_functions`: Curve table from the calculation file, empty for a constant model
- `t`: Time to seed time-driven boundaries at, 0 for a fresh run

# Note
- Modifies global variables: `C_g`, `T`, `P`, `C_lime`
- Call this after `zero_variables!()` to set up the initial state
"""
function apply_all_initial_conditions!(mesh, materials,
                                       time_functions = Dict{Int, TimeFunction}(), t = 0.0)
    apply_initial_concentrations!(mesh)
    apply_initial_temperature!(mesh)
    apply_temperature_bc!(mesh, time_functions, t)
    apply_concentration_bc!(mesh, time_functions, t)
    apply_partial_pressure_bc!(mesh, time_functions, t)
    apply_pressure_bc!(mesh, time_functions, t)
    apply_initial_lime_concentration!(mesh, materials)
    
    println("\nAll initial conditions and BCs applied successfully")
end


