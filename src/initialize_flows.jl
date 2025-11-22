#------------------------------------------------------------------------------
# ADSIM Flow Initialization Module
# This module contains functions to initialize flow vectors for ADSIM
# FEM calculations following the lumped mass explicit formulation
#------------------------------------------------------------------------------

using Base.Threads

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

# Lumped mass matrix (same for all gases)
global M_L::Vector{Float64} = Float64[]

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
  `q_boundary`, `q_source_sink`, and `M_L`
- All flow vectors are initialized to zero
- q_source_sink is only a vector (not matrix) as it only applies to CO2
"""
function zero_flow_vectors!(Nnodes::Int, NGases::Int)
    global q_advection, q_gravitational, q_diffusion, q_boundary, q_source_sink, M_L
    
    # Initialize flow matrices (Nnodes × NGases)
    q_advection = zeros(Float64, Nnodes, NGases)
    q_gravitational = zeros(Float64, Nnodes, NGases)
    q_diffusion = zeros(Float64, Nnodes, NGases)
    q_boundary = zeros(Float64, Nnodes, NGases)
    
    # Initialize source/sink vector (only for CO2)
    q_source_sink = zeros(Float64, Nnodes)
    
    # Initialize lumped mass matrix
    M_L = zeros(Float64, Nnodes)

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
- `materials`: Material data structure containing gas dictionary

# Note
- Modifies global variable `q_boundary`
- Flow values need to be weighted by nodal influence length
- Flow is positive when entering the domain, negative when leaving
"""
function apply_boundary_flows!(mesh, materials)
    global q_boundary, NGases
    
    # Initialize boundary flow counter
    num_bc_nodes = 0
    
    # Calculate influence length for each boundary node
    influence_lengths = calculate_boundary_influence_lengths(mesh)
    
    # Apply nodal uniform flow boundary conditions
    for (node_id, flows) in mesh.uniform_flow_bc
        influence_length = get(influence_lengths, node_id, 0.0)
        @threads for gas_idx in 1:NGases
            # Weight flows by nodal influence length
            q_boundary[node_id, gas_idx] = flows[gas_idx] * influence_length
        end
        num_bc_nodes += 1
    end
    
    if num_bc_nodes > 0
        println("\nBoundary flows applied to ", num_bc_nodes, " nodes")
    else
        println("\nNo boundary flow conditions to apply")
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
    
    for (node_id, _) in mesh.uniform_flow_bc
        # Find all boundary elements connected to this node
        connected_lengths = Float64[]
        
        # Search through boundary elements to find connections
        # This assumes boundary elements are edges with 2 nodes
        for boundary_elem in mesh.boundary_elements
            if node_id in boundary_elem
                # Calculate edge length
                node1, node2 = boundary_elem
                x1, y1 = mesh.coordinates[node1, :]
                x2, y2 = mesh.coordinates[node2, :]
                edge_length = sqrt((x2 - x1)^2 + (y2 - y1)^2)
                push!(connected_lengths, edge_length)
            end
        end
        
        # Influence length is sum of half-lengths of connected edges
        influence_lengths[node_id] = sum(connected_lengths) / 2.0
    end
    
    return influence_lengths
end


#------------------------------------------------------------------------------
# Calculate lumped mass matrix
#------------------------------------------------------------------------------
"""
    calculate_lumped_mass!(mesh::MeshData, materials)

Calculate the lumped mass matrix for all nodes based on element areas and 
material porosity. The lumped mass formulation concentrates element mass at nodes.

For 2D quadrilateral elements with 4 nodes:
- Element area is divided equally among the 4 nodes
- Mass at each node = (element_area / 4) × porosity

# Arguments
- `mesh::MeshData`: Mesh data structure containing coordinates and connectivity
- `materials`: Material data structure containing soil properties with porosity

# Note
- Modifies global variable `M_L`
- Uses bilinear quadrilateral element formulation
- Assumes constant porosity within each element
"""
function calculate_lumped_mass!(mesh, materials)
    global M_L
    
    # Loop through all elements
    for elem_id in 1:mesh.num_elements
        # Get element nodes
        element_nodes = get_element_nodes(mesh, elem_id)
        
        # Get material index for this element
        material_idx = get_element_material(mesh, elem_id)
        
        if material_idx !== nothing
            # Get the soil name from the soil dictionary
            soil_name = materials.soil_dictionary[material_idx]
            
            # Get the soil properties for this material
            soil_props = get_soil_properties(materials, soil_name)
            
            if soil_props !== nothing
                # Get porosity from material
                n = soil_props.porosity
                
                # Calculate element area using shoelace formula
                # For quadrilateral: Area = 0.5 * |sum of cross products|
                x = [mesh.coordinates[node, 1] for node in element_nodes]
                y = [mesh.coordinates[node, 2] for node in element_nodes]
                
                # Shoelace formula for quadrilateral (ordered nodes)
                area = 0.5 * abs(
                    (x[1]*y[2] - x[2]*y[1]) + 
                    (x[2]*y[3] - x[3]*y[2]) + 
                    (x[3]*y[4] - x[4]*y[3]) + 
                    (x[4]*y[1] - x[1]*y[4])
                )
                
                # Lumped mass contribution per node (equal distribution)
                nodal_mass = (area / 4.0) * n
                
                # Add contribution to each node of the element
                for node_id in element_nodes
                    M_L[node_id] += nodal_mass
                end
            end
        end
    end
    
    println("\nLumped mass matrix calculated")
    println("  - Total system mass (porous volume): ", sum(M_L))
    println("  - Min nodal mass: ", minimum(M_L))
    println("  - Max nodal mass: ", maximum(M_L))
end


#------------------------------------------------------------------------------
# Reset flow vectors (called at each time step)
#------------------------------------------------------------------------------
"""
    reset_flows!()

Reset all flow vectors to zero. This function should be called at the beginning
of each time step before calculating new flows.

# Note
- Modifies global flow vectors: `q_advection`, `q_gravitational`, `q_diffusion`
- Does NOT reset `q_boundary` (boundary conditions remain constant)
- Does NOT reset `q_source_sink` until reaction module is implemented
- Does NOT reset `M_L` (mass matrix is constant)
"""
function reset_flows!()
    global q_advection, q_gravitational, q_diffusion
    
    # Reset calculated flows (boundary flows remain constant)
    q_advection .= 0.0
    q_gravitational .= 0.0
    q_diffusion .= 0.0
end


#------------------------------------------------------------------------------
# Initialize all flows
#------------------------------------------------------------------------------
"""
    initialize_all_flows!(mesh::MeshData, materials)

Complete flow initialization procedure. This function performs all necessary
steps to initialize flow vectors and lumped mass matrix.

Call sequence:
1. Zero all flow vectors
2. Calculate lumped mass matrix
3. Apply boundary flow conditions

# Arguments
- `mesh::MeshData`: Mesh data structure
- `materials`: Material data structure

# Note
- Should be called once during initialization after mesh and materials are loaded
- Requires global variables `Nnodes` and `NGases` to be set
"""
function initialize_all_flows!(mesh, materials, Nnodes::Int, NGases::Int)
    zero_flow_vectors!(Nnodes, NGases)
    calculate_lumped_mass!(mesh, materials)
    apply_boundary_flows!(mesh, materials)
    
    println("\n" * "="^70)
    println("Flow initialization complete")
    println("="^70)
end


# Export all public functions
export zero_flow_vectors!, apply_boundary_flows!, calculate_lumped_mass!
export reset_flows!, initialize_all_flows!
export q_advection, q_gravitational, q_diffusion, q_boundary, q_source_sink, M_L
