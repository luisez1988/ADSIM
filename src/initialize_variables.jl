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

# Reactive species
global C_lime::Vector{Float64} = Float64[]
global C_caco3::Vector{Float64} = Float64[]
global C_lime_residual::Vector{Float64} = Float64[]

# Time derivatives
global dC_g_dt::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
global dT_dt::Vector{Float64} = Float64[]
global dC_lime_dt::Vector{Float64} = Float64[]

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
function zero_variables!(mesh, materials)
    global NDim, Nnodes, Nelements, NSoils, NGases
    global C_g, P, T, v, P_boundary
    global C_lime, C_caco3, C_lime_residual
    global dC_g_dt, dT_dt, dC_lime_dt
    
    # Set dimensions
    NDim = 2  # Number of spatial dimensions - TODO: generalize for 3D
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
    
    # Allocate and initialize reactive species
    C_lime = zeros(Float64, Nnodes)
    C_caco3 = zeros(Float64, Nnodes)
    C_lime_residual = zeros(Float64, NSoils)
    
    # Allocate and initialize time derivatives
    dC_g_dt = zeros(Float64, Nnodes, NGases)
    dT_dt = zeros(Float64, Nnodes)
    dC_lime_dt = zeros(Float64, Nnodes)    

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
    
    # Apply element-based initial concentrations to nodes
    for (elem_id, concentrations) in mesh.initial_concentrations
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
    
    # Apply element-based initial temperatures to nodes
    for (elem_id, temperature) in mesh.initial_temperature
        # Get nodes of this element
        element_nodes = get_element_nodes(mesh, elem_id)
        
        # Apply temperature to each node of the element
        for node_id in element_nodes
            T[node_id] = temperature
        end
    end
end


"""
    apply_concentration_bc!(mesh::MeshData)

Apply concentration boundary conditions from mesh data to the global C_g array.
This function sets fixed gas concentrations at nodes where concentration 
boundary conditions are specified.

# Arguments
- `mesh::MeshData`: Mesh data structure containing concentration BC data

# Note
- Modifies global variable `C_g`
- These values should be maintained throughout the simulation for BC nodes
"""
function apply_concentration_bc!(mesh)
    global C_g, NGases, P_boundary
    
    # Apply nodal concentration boundary conditions
    for (node_id, concentrations) in mesh.concentration_bc
        @threads for gas_idx in 1:NGases
            C_g[node_id, gas_idx] = concentrations[gas_idx]
            P_boundary[node_id, gas_idx] = 0  # Mark node as having concentration BC
        end
    end
end


"""
    apply_pressure_bc!(mesh::MeshData)

Apply absolute pressure boundary conditions from mesh data to the global P array.
This function sets fixed pressures at nodes where pressure boundary conditions 
are specified and marks the vacating gas as having restricted flow.

# Arguments
- `mesh::MeshData`: Mesh data structure containing pressure BC data

# Note
- Modifies global variables `P` and `P_boundary`
- Sets P_boundary to 0 for the vacating gas at pressure BC nodes
- These values should be maintained throughout the simulation for BC nodes
"""
function apply_pressure_bc!(mesh)
    global P, P_boundary
    
    # Apply nodal pressure boundary conditions
    for (node_id, pressure) in mesh.absolute_pressure_bc
        P[node_id] = pressure
        
        # Get the vacating gas index for this node and mark it as restricted
        if haskey(mesh.vacating_gas_bc, node_id)
            gas_idx = mesh.vacating_gas_bc[node_id]
            P_boundary[node_id, gas_idx] = 0  # Gas can vacate the system
        end
    end
end


"""
    apply_initial_lime_concentration!(mesh::MeshData, materials)

Apply initial lime concentrations from material properties to the global C_lime array.
This function loops through all elements, gets their assigned material, retrieves
the lime content from that material, and assigns it to all nodes in the element.

# Arguments
- `mesh::MeshData`: Mesh data structure containing element-material assignments
- `materials`: Material data structure containing soil properties with lime content

# Note
- Modifies global variable `C_lime`
- If a node belongs to multiple elements with different materials,
  the last element's value will be used
"""
function apply_initial_lime_concentration!(mesh, materials)
    global C_lime, C_lime_residual
    
    # Loop through all elements
    for elem_id in 1:mesh.num_elements
        # Get material index for this element
        material_idx = get_element_material(mesh, elem_id)
        
        if material_idx !== nothing
            # Get the soil name from the soil dictionary
            soil_name = materials.soil_dictionary[material_idx]
            
            # Get the soil properties for this material
            soil_props = get_soil_properties(materials, soil_name)
            
            if soil_props !== nothing
                # Get lime content from material
                β_l = soil_props.lime_content
                G_s = soil_props.specific_gravity
                n=soil_props.porosity                
                M_lime=74.093   # Molar mass of Ca(OH)2 in g/mol
                #Calculate lime concentration in mol/m^3 
                lime_concentration= (β_l * G_s * (1 - n) * 1e6 ) / M_lime #Asumes ρ_w= 1000 kg/m^3  

                #Calculatte reidual lime 
                residual_percent= soil_props.residual_lime
                C_lime_residual[material_idx] = residual_percent * lime_concentration
                
                # Get nodes of this element
                element_nodes = get_element_nodes(mesh, elem_id)
                
                # Assign lime content to each node of the element
                for node_id in element_nodes
                    C_lime[node_id] = lime_concentration
                end
            end
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

# Note
- Modifies global variables: `C_g`, `T`, `P`, `C_lime`
- Call this after `zero_variables!()` to set up the initial state
"""
function apply_all_initial_conditions!(mesh, materials)
    apply_initial_concentrations!(mesh)
    apply_initial_temperature!(mesh)
    apply_concentration_bc!(mesh)
    apply_pressure_bc!(mesh)
    apply_initial_lime_concentration!(mesh, materials)
    
    println("\nAll initial conditions and BCs applied successfully")
end


