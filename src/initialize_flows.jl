#------------------------------------------------------------------------------
# ADSIM Flow Initialization Module
# This module contains functions to initialize flow vectors for ADSIM
# FEM calculations following the lumped mass explicit formulation
#------------------------------------------------------------------------------

using Base.Threads
using LinearAlgebra

#------------------------------------------------------------------------------
# Global flow variables - shared across all modules
#------------------------------------------------------------------------------

# Flow vectors for each gas species (Nnodes × NGases)
# Each flow type contributes to the total concentration change
global q_advection::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
global q_gravitational::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
global q_diffusion::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
global q_boundary::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
global q_source_sink::Vector{Float64} = Float64[]  # Only for CO2

# Boundary edge geometry - precomputed during initialization
# Each tuple contains: (element_id, node_i, node_j, edge_length, outward_normal)
global boundary_edges::Vector{Tuple{Int, Int, Int, Float64, Vector{Float64}}} = []

#------------------------------------------------------------------------------
# Initialize flow vectors
#------------------------------------------------------------------------------
"""
    zero_flow_vectors!(Nnodes::Int, NGases::Int)

Initialize all flow vectors to zero based on mesh dimensions and number of gases.
This function should be called after mesh and material data are loaded.

# Arguments
- `Nnodes::Int`: Total number of nodes in the mesh
- `NGases::Int`: Total number of gas species

# Note
- Modifies global flow vectors: `q_advection`, `q_gravitational`, `q_diffusion`, 
  `q_boundary`, `q_source_sink`
- All flow vectors are initialized to zero
- q_source_sink is only a vector (not matrix) as it only applies to CO2
"""
function zero_flow_vectors!(Nnodes::Int, NGases::Int)
    global q_advection, q_gravitational, q_diffusion, q_boundary, q_source_sink
    
    # Initialize flow matrices (Nnodes × NGases)
    q_advection = zeros(Float64, Nnodes, NGases)
    q_gravitational = zeros(Float64, Nnodes, NGases)
    q_diffusion = zeros(Float64, Nnodes, NGases)
    q_boundary = zeros(Float64, Nnodes, NGases)
    
    # Initialize source/sink vector (only for CO2)
    q_source_sink = zeros(Float64, Nnodes)
end


#------------------------------------------------------------------------------
# Apply boundary flow conditions
#------------------------------------------------------------------------------
"""
    apply_boundary_flows!(mesh::MeshData, materials)

Apply uniform flow boundary conditions from mesh data to the global q_boundary array.
The GID interface permits assignment of uniform steady normal flows at nodes.
These values must consider the length of influence associated with each node.

For 1D boundaries:
- Internal boundary nodes: influence length = h (full element length)
- Corner/end nodes: influence length = h/2 (half element length)

# Arguments
- `mesh::MeshData`: Mesh data structure containing uniform flow BC data

# Note
- Modifies global variable `q_boundary`
- Flow values need to be weighted by nodal influence length
- Flow is positive when entering the domain, negative when leaving
"""
function apply_boundary_flows!(mesh)
    global q_boundary, NGases    
    
    # Calculate influence length for each boundary node
    influence_lengths = calculate_boundary_influence_lengths(mesh)
    
    # Apply nodal uniform flow boundary conditions
    for (node_id, flows) in mesh.uniform_flow_bc        
        @threads for gas_idx in 1:NGases
            # Weight flows by nodal influence length
            q_boundary[node_id, gas_idx] = flows[gas_idx] * influence_lengths[node_id]
        end        
    end
end

"""
calculate_boundary_influence_lengths(mesh::MeshData)

Calculate the influence length for each boundary node based on connected boundary elements.

# Returns
- Dictionary mapping node_id to influence length
"""
function calculate_boundary_influence_lengths(mesh)
    influence_lengths = Dict{Int, Float64}()
    #get a list of all nodes with boundary conditions
    boundary_nodes = collect(keys(mesh.uniform_flow_bc))
    for node_id in boundary_nodes
        # Find all elements this node belongs to
        connected_elements = get_node_elements(mesh, node_id)
        total_length = 0.0
        # loop remaining nodes and check if they share an element
        for elem_id in connected_elements
            element_nodes = get_element_nodes(mesh, elem_id)
            for xnode_id in element_nodes
                if xnode_id != node_id && xnode_id in boundary_nodes
                    #calculate distance between nodes
                    x_1 = mesh.coordinates[node_id, :]                    
                    x_2 = mesh.coordinates[xnode_id, :]
                    length = norm(x_2 - x_1)
                    total_length += length
                end
            end
        end
        influence_lengths[node_id] = total_length / 2.0  # divide by 2 to account for half-length at ends
    end
    
    return influence_lengths
end

#------------------------------------------------------------------------------
# Initialize all flows
#------------------------------------------------------------------------------
"""
    initialize_all_flows!(mesh::MeshData, materials, Nnodes::Int, NGases::Int)

Complete flow initialization procedure. This function performs all necessary
steps to initialize flow vectors and lumped mass matrix.

Call sequence:
1. Zero all flow vectors
2. Identify and store boundary edges with geometry
3. Apply boundary flow conditions

# Arguments
- `mesh::MeshData`: Mesh data structure
- `materials`: Material data structure (currently unused, reserved for future use)
- `Nnodes::Int`: Total number of nodes in the mesh
- `NGases::Int`: Total number of gas species

# Note
- Should be called once during initialization after mesh and materials are loaded
- Precomputes boundary edge geometry (length and outward normals) for efficiency
"""
function initialize_all_flows!(mesh, materials, Nnodes::Int, NGases::Int)
    global boundary_edges
    
    zero_flow_vectors!(Nnodes, NGases)
    
    # Precompute boundary edge geometry for pressure BC nodes
    # This is done once during initialization to avoid repeated calculation
    boundary_edges = identify_boundary_edges(mesh)
    
    apply_boundary_flows!(mesh)
end


# Export all public functions
export zero_flow_vectors!, apply_boundary_flows!, initialize_all_flows!
export q_advection, q_gravitational, q_diffusion, q_boundary, q_source_sink
export boundary_edges
