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

# Reactive species
global C_lime::Vector{Float64} = Float64[]
global C_caco3::Vector{Float64} = Float64[]

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
    global C_g, P, T, v
    global C_lime, C_caco3
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
    
    # Allocate and initialize reactive species
    C_lime = zeros(Float64, Nnodes)
    C_caco3 = zeros(Float64, Nnodes)
    
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
    global C_g, NGases
    
    # Apply nodal concentration boundary conditions
    for (node_id, concentrations) in mesh.concentration_bc
        @threads for gas_idx in 1:NGases
            C_g[node_id, gas_idx] = concentrations[gas_idx]            
        end
    end
end


"""
    apply_pressure_bc!(mesh::MeshData)

Apply absolute pressure boundary conditions from mesh data to the global P array.
This function sets fixed pressures at nodes where pressure boundary conditions 
are specified.

# Arguments
- `mesh::MeshData`: Mesh data structure containing pressure BC data

# Note
- Modifies global variable `P`
- These values should be maintained throughout the simulation for BC nodes
"""
function apply_pressure_bc!(mesh)
    global P
    
    # Apply nodal pressure boundary conditions
    for (node_id, pressure) in mesh.absolute_pressure_bc
        P[node_id] = pressure
    end
end


"""
    apply_all_initial_conditions!(mesh::MeshData)

Apply all initial conditions and boundary conditions from mesh data.
This is a convenience function that calls all individual application functions.

# Arguments
- `mesh::MeshData`: Mesh data structure containing all initial and boundary condition data

# Note
- Modifies global variables: `C_g`, `T`, `P`
- Call this after `initialize_variables!()` to set up the initial state
"""
function apply_all_initial_conditions!(mesh)
    apply_initial_concentrations!(mesh)
    apply_initial_temperature!(mesh)
    apply_concentration_bc!(mesh)
    apply_pressure_bc!(mesh)
    
    println("\nAll initial conditions and BCs applied successfully")
end


