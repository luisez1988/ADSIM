#------------------------------------------------------------------------------
# ADSIM Flow Initialization Module
# Flow vector initialization for ADSIM FEM calculations (lumped mass explicit formulation)
#------------------------------------------------------------------------------

using Base.Threads
using LinearAlgebra

#------------------------------------------------------------------------------
# Global flow variables - shared across all modules
#------------------------------------------------------------------------------

# Flow vectors per gas species (Nnodes × NGases)
global q_advection::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
global q_gravitational::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
global q_diffusion::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
global q_boundary::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
global q_source_sink::Vector{Float64} = Float64[]  # Only for CO2

# Boundary edge geometry, precomputed at initialization
# Tuple: (element_id, node_i, node_j, edge_length, outward_normal)
global boundary_edges::Vector{Tuple{Int, Int, Int, Float64, Vector{Float64}}} = []

# Tributary edge length per uniform flow BC node, computed once at initialization.
# Cached: time-driven flow boundaries refill q_boundary every step, and
# calculate_boundary_influence_lengths loops over the element connectivity of every
# boundary node (too expensive per step).
global flow_node_influences::Dict{Int, Float64} = Dict{Int, Float64}()

#------------------------------------------------------------------------------
# Initialize flow vectors
#------------------------------------------------------------------------------
"""
    zero_flow_vectors!(Nnodes::Int, NGases::Int)

Allocate all flow vectors as zeros. Call after mesh and material data are loaded.

Mutates the globals `q_advection`, `q_gravitational`, `q_diffusion`, `q_boundary`
(Nnodes × NGases) and `q_source_sink` (Nnodes, CO2 only).
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

Apply the mesh uniform flow BCs at time `t`, weighted by nodal influence length.
Mutates the global `q_boundary`.

GID assigns uniform normal flows at nodes. Influence length on 1D boundaries:
- Internal boundary nodes: h (full element length)
- Corner/end nodes: h/2

Sign: positive into the domain, negative out.
"""
function apply_boundary_flows!(mesh, time_functions = Dict{Int, TimeFunction}(), t = 0.0)
    global q_boundary, NGases, flow_node_influences

    # Influence lengths from initialize_all_flows!; computed here for standalone calls
    if isempty(flow_node_influences) && !isempty(mesh.uniform_flow_bc)
        flow_node_influences = calculate_boundary_influence_lengths(mesh)
    end

    # Apply nodal uniform flow boundary conditions
    for (node_id, flows) in mesh.uniform_flow_bc
        le = flow_node_influences[node_id]
        tf_ids = get(mesh.uniform_flow_bc_tf, node_id, Int[])
        @threads for gas_idx in 1:NGases
            tf_id = gas_idx <= length(tf_ids) ? tf_ids[gas_idx] : 0
            # Weight flows by nodal influence length
            q_boundary[node_id, gas_idx] =
                bc_value(time_functions, tf_id, flows[gas_idx], t) * le
        end
    end
end

"""
calculate_boundary_influence_lengths(mesh::MeshData; axisymmetric = false)

Compute the influence length of each uniform flow BC node from connected boundary elements.

# Arguments
- `axisymmetric`: Weight each edge by its midpoint radius, giving the ring measure
  `l_e * r̄_e` instead of the edge length (same as `get_boundary_node_influences` in
  `read_mesh.jl`).

# Returns
- node_id => influence length [m] (plane) or ring measure [m²] (axisymmetric)
"""
function calculate_boundary_influence_lengths(mesh; axisymmetric::Bool = false)
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
                    # Ring measure instead of arclength; r = first coordinate
                    if axisymmetric
                        length *= 0.5 * (x_1[1] + x_2[1])
                    end
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

Initialize all flow vectors and boundary geometry. Call once after mesh and materials
are loaded. Mutates the flow globals, `boundary_edges` and `flow_node_influences`.

Sequence:
1. Zero all flow vectors
2. Identify boundary edges and their geometry (length, outward normal)
3. Apply boundary flow conditions

# Arguments
- `materials`: Currently unused
"""
function initialize_all_flows!(mesh, materials, Nnodes::Int, NGases::Int,
                              time_functions = Dict{Int, TimeFunction}(), t = 0.0;
                              axisymmetric::Bool = false)
    global boundary_edges, flow_node_influences

    zero_flow_vectors!(Nnodes, NGases)

    # Boundary edge geometry of pressure BC nodes, computed once
    boundary_edges = identify_boundary_edges(mesh)

    # Flow boundary influence lengths, computed once (time-driven fluxes rebuild
    # q_boundary every step). Axisymmetric: ring measure l_e * r̄_e (Eq. ax_external_flux).
    flow_node_influences = calculate_boundary_influence_lengths(mesh; axisymmetric = axisymmetric)

    apply_boundary_flows!(mesh, time_functions, t)
end


# Export all public functions
export zero_flow_vectors!, apply_boundary_flows!, initialize_all_flows!
export q_advection, q_gravitational, q_diffusion, q_boundary, q_source_sink
export boundary_edges, flow_node_influences, calculate_boundary_influence_lengths
