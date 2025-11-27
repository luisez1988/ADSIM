#------------------------------------------------------------------------------
# ADSIM Mesh Reader Module
# This module contains functions to read and parse .mesh files for ADSIM
# FEM calculations
#------------------------------------------------------------------------------

"""
MeshData

Structure to store all mesh data and associated boundary/initial conditions.

# Fields
- `num_nodes::Int`: Total number of nodes in the mesh
- `num_elements::Int`: Total number of elements in the mesh
- `coordinates::Matrix{Float64}`: Nodal coordinates (num_nodes × 2)
- `elements::Matrix{Int}`: Element connectivity (num_elements × 4)
- `concentration_bc::Dict{Int, Vector{Float64}}`: Concentration BC (node_id => [gas1, gas2, ...])
- `uniform_flow_bc::Dict{Int, Vector{Float64}}`: Uniform flow BC (node_id => [gas1, gas2, ...])
- `absolute_pressure_bc::Dict{Int, Float64}`: Absolute pressure BC (node_id => pressure)
# - `vacating_gas_bc::Dict{Int, Int}`: Vacating gas index for pressure BC (node_id => gas_index)
- `initial_concentrations::Dict{Int, Vector{Float64}}`: Initial concentrations (elem_id => [gas1, gas2, ...])
- `initial_temperature::Dict{Int, Float64}`: Initial temperature (elem_id => temperature)
- `materials::Dict{Int, Int}`: Material assignment (elem_id => material_index)
"""

mutable struct MeshData
    num_nodes::Int
    num_elements::Int
    coordinates::Matrix{Float64}
    elements::Matrix{Int}
    concentration_bc::Dict{Int, Vector{Float64}}
    uniform_flow_bc::Dict{Int, Vector{Float64}}
    absolute_pressure_bc::Dict{Int, Float64}
    # vacating_gas_bc::Dict{Int, Int}
    initial_concentrations::Dict{Int, Vector{Float64}}
    initial_temperature::Dict{Int, Float64}
    materials::Dict{Int, Int}
    
    function MeshData()
        new(0, 0, 
            zeros(Float64, 0, 0), 
            zeros(Int, 0, 0),
            Dict{Int, Vector{Float64}}(),
            Dict{Int, Vector{Float64}}(),
            Dict{Int, Float64}(),
            # Dict{Int, Int}(),
            Dict{Int, Vector{Float64}}(),
            Dict{Int, Float64}(),
            Dict{Int, Int}())
    end
end


"""
read_mesh_file(filename::String) -> MeshData

Read an ADSIM .mesh file and return a MeshData structure containing all
geometrical and physical data.

# Arguments
- `filename::String`: Path to the .mesh file

# Returns
- `MeshData`: Structure containing all mesh information

# Example
```julia
mesh = read_mesh_file("problem.mesh")
println("Number of nodes: ", mesh.num_nodes)
println("Number of elements: ", mesh.num_elements)
```
"""
function read_mesh_file(filename::String)
    mesh = MeshData()
    
    open(filename, "r") do file
        lines = readlines(file)
        line_idx = 1
        
        while line_idx <= length(lines)
            line = strip(lines[line_idx])
            
            # Skip header and empty lines
            if startswith(line, "###") || isempty(line)
                line_idx += 1
                continue
            end
            
            # Parse MESH counters
            if startswith(line, "MESH")
                line_idx = parse_mesh_counters!(mesh, lines, line_idx)
                
            # Parse coordinates
            elseif line == "coordinates"
                line_idx = parse_coordinates!(mesh, lines, line_idx + 1)
                
            # Parse elements
            elseif line == "elements"
                line_idx = parse_elements!(mesh, lines, line_idx + 1)
                
            # Parse concentration boundary conditions
            elseif line == "concentration_bc"
                line_idx = parse_concentration_bc!(mesh, lines, line_idx + 1)
                
            # Parse uniform flow boundary conditions
            elseif line == "uniform_flow_bc"
                line_idx = parse_uniform_flow_bc!(mesh, lines, line_idx + 1)
                
            # Parse absolute pressure boundary conditions
            elseif line == "absolute_pressure"
                line_idx = parse_absolute_pressure!(mesh, lines, line_idx + 1)
                
            # Parse initial concentrations
            elseif line == "initial_concentrations"
                line_idx = parse_initial_concentrations!(mesh, lines, line_idx + 1)
                
            # Parse initial temperature
            elseif line == "initial_temperature"
                line_idx = parse_initial_temperature!(mesh, lines, line_idx + 1)
                
            # Parse materials
            elseif line == "materials"
                line_idx = parse_materials!(mesh, lines, line_idx + 1)
                
            else
                line_idx += 1
            end
        end
    end
    
    return mesh
end


"""
parse_mesh_counters!(mesh::MeshData, lines::Vector{String}, line_idx::Int) -> Int

Parse the MESH counter line and store node/element counts.

# Arguments
- `mesh::MeshData`: Mesh data structure to populate
- `lines::Vector{String}`: All lines from the file
- `line_idx::Int`: Current line index

# Returns
- `Int`: Next line index to process
"""
function parse_mesh_counters!(mesh::MeshData, lines::Vector{String}, line_idx::Int)
    parts = split(strip(lines[line_idx]))
    mesh.num_nodes = parse(Int, parts[2])
    mesh.num_elements = parse(Int, parts[3])
    return line_idx + 1
end


"""
parse_coordinates!(mesh::MeshData, lines::Vector{String}, line_idx::Int) -> Int

Parse nodal coordinates from the mesh file.

# Arguments
- `mesh::MeshData`: Mesh data structure to populate
- `lines::Vector{String}`: All lines from the file
- `line_idx::Int`: Current line index

# Returns
- `Int`: Next line index to process
"""
function parse_coordinates!(mesh::MeshData, lines::Vector{String}, line_idx::Int)
    mesh.coordinates = zeros(Float64, mesh.num_nodes, 2)
    
    for i in 1:mesh.num_nodes # loop over node numbers
        line = strip(lines[line_idx])
        coords = split(line)
        mesh.coordinates[i, 1] = parse(Float64, coords[1]) #read x coordinate
        mesh.coordinates[i, 2] = parse(Float64, coords[2]) #read y coordinate
        line_idx += 1 # move to next line
    end
    if strip(lines[line_idx]) == "end coordinates" # check for end marker
        line_idx += 1 # move to next line
    end    
    return line_idx
end


"""
parse_elements!(mesh::MeshData, lines::Vector{String}, line_idx::Int) -> Int

Parse element connectivity from the mesh file.

# Arguments
- `mesh::MeshData`: Mesh data structure to populate
- `lines::Vector{String}`: All lines from the file
- `line_idx::Int`: Current line index

# Returns
- `Int`: Next line index to process
"""
function parse_elements!(mesh::MeshData, lines::Vector{String}, line_idx::Int)
    mesh.elements = zeros(Int, mesh.num_elements, 4)
    
    for i in 1:mesh.num_elements
        line = strip(lines[line_idx])
        nodes = split(line)
        for j in 1:4
            mesh.elements[i, j] = parse(Int, nodes[j])
        end
        line_idx += 1
    end
    
    if strip(lines[line_idx]) == "end elements"
        line_idx += 1
    end
    
    return line_idx
end


"""
parse_concentration_bc!(mesh::MeshData, lines::Vector{String}, line_idx::Int) -> Int

Parse concentration boundary conditions from the mesh file.

# Arguments
- `mesh::MeshData`: Mesh data structure to populate
- `lines::Vector{String}`: All lines from the file
- `line_idx::Int`: Current line index

# Returns
- `Int`: Next line index to process
"""
function parse_concentration_bc!(mesh::MeshData, lines::Vector{String}, line_idx::Int)
    # Read counter
    counter = parse(Int, strip(lines[line_idx]))
    line_idx += 1
    
    # Read boundary conditions
    for i in 1:counter
        line = strip(lines[line_idx])
        parts = split(line)
        node_id = parse(Int, parts[1])
        concentrations = [parse(Float64, parts[i]) for i in 2:length(parts)]
        mesh.concentration_bc[node_id] = concentrations
        line_idx += 1
    end
    
    # Skip end marker
    if strip(lines[line_idx]) == "end concentration_bc"
        line_idx += 1
    end
    
    return line_idx
end


"""
parse_uniform_flow_bc!(mesh::MeshData, lines::Vector{String}, line_idx::Int) -> Int

Parse uniform flow boundary conditions from the mesh file.

# Arguments
- `mesh::MeshData`: Mesh data structure to populate
- `lines::Vector{String}`: All lines from the file
- `line_idx::Int`: Current line index

# Returns
- `Int`: Next line index to process
"""
function parse_uniform_flow_bc!(mesh::MeshData, lines::Vector{String}, line_idx::Int)
    # Read counter
    counter = parse(Int, strip(lines[line_idx]))
    line_idx += 1
    
    # Read boundary conditions
    for _ in 1:counter #counter is 0 if no BCs
        line = strip(lines[line_idx])
        parts = split(line)
        node_id = parse(Int, parts[1])
        flows = [parse(Float64, parts[i]) for i in 2:lastindex(parts)]
        mesh.uniform_flow_bc[node_id] = flows
        line_idx += 1
    end

    if strip(lines[line_idx]) == "end uniform_flow_bc"
        line_idx += 1
    end
    
    return line_idx
end


"""
parse_absolute_pressure!(mesh::MeshData, lines::Vector{String}, line_idx::Int) -> Int

Parse absolute pressure boundary conditions from the mesh file.

# Arguments
- `mesh::MeshData`: Mesh data structure to populate
- `lines::Vector{String}`: All lines from the file
- `line_idx::Int`: Current line index

# Returns
- `Int`: Next line index to process
"""
function parse_absolute_pressure!(mesh::MeshData, lines::Vector{String}, line_idx::Int)
    # Read counter
    counter = parse(Int, strip(lines[line_idx]))
    line_idx += 1
    
    # Read boundary conditions
    for _ in 1:counter
        line = strip(lines[line_idx])
        parts = split(line)
        node_id = parse(Int, parts[1])
        pressure = parse(Float64, parts[2])
        # gas_index = parse(Int, parts[3])
        mesh.absolute_pressure_bc[node_id] = pressure
        # mesh.vacating_gas_bc[node_id] = gas_index
        line_idx += 1
    end
    
    # Skip end marker
    if strip(lines[line_idx]) == "end absolute_pressure"
        line_idx += 1
    end
    
    return line_idx
end


"""
parse_initial_concentrations!(mesh::MeshData, lines::Vector{String}, line_idx::Int) -> Int

Parse initial gas concentrations from the mesh file.

# Arguments
- `mesh::MeshData`: Mesh data structure to populate
- `lines::Vector{String}`: All lines from the file
- `line_idx::Int`: Current line index

# Returns
- `Int`: Next line index to process
"""
function parse_initial_concentrations!(mesh::MeshData, lines::Vector{String}, line_idx::Int)
    # Read counter
    counter = parse(Int, strip(lines[line_idx]))
    line_idx += 1
    
    # Read initial conditions
    for _ in 1:counter
        line = strip(lines[line_idx])
        parts = split(line)
        elem_id = parse(Int, parts[1])
        concentrations = [parse(Float64, parts[i]) for i in 2:length(parts)]
        mesh.initial_concentrations[elem_id] = concentrations
        line_idx += 1
    end
    
    # Skip end marker
    if strip(lines[line_idx]) == "end initial_concentrations"
        line_idx += 1
    end
    
    return line_idx
end


"""
parse_initial_temperature!(mesh::MeshData, lines::Vector{String}, line_idx::Int) -> Int

Parse initial temperature from the mesh file.

# Arguments
- `mesh::MeshData`: Mesh data structure to populate
- `lines::Vector{String}`: All lines from the file
- `line_idx::Int`: Current line index

# Returns
- `Int`: Next line index to process
"""
function parse_initial_temperature!(mesh::MeshData, lines::Vector{String}, line_idx::Int)
    # Read counter
    counter = parse(Int, strip(lines[line_idx]))
    line_idx += 1
    
    # Read initial conditions
    for _ in 1:counter
        line = strip(lines[line_idx])
        parts = split(line)
        elem_id = parse(Int, parts[1])
        temperature = parse(Float64, parts[2])
        mesh.initial_temperature[elem_id] = temperature
        line_idx += 1
    end
    
    # Skip end marker
    if strip(lines[line_idx]) == "end initial_temperature"
        line_idx += 1
    end
    
    return line_idx
end


"""
parse_materials!(mesh::MeshData, lines::Vector{String}, line_idx::Int) -> Int

Parse material assignments from the mesh file.

# Arguments
- `mesh::MeshData`: Mesh data structure to populate
- `lines::Vector{String}`: All lines from the file
- `line_idx::Int`: Current line index

# Returns
- `Int`: Next line index to process
"""
function parse_materials!(mesh::MeshData, lines::Vector{String}, line_idx::Int)
    # Read counter
    counter = parse(Int, strip(lines[line_idx]))
    line_idx += 1
    
    # Read material assignments
    for mat_count in 1:counter
        line = strip(lines[line_idx])

        parts = split(line)
        elem_id = parse(Int, parts[1])
        material_idx = parse(Int, parts[2])
        mesh.materials[elem_id] = material_idx
        
        line_idx += 1
    end
    

    if strip(lines[line_idx]) == "end materials"
        return line_idx + 1
    end
    
    return line_idx
end


"""
get_element_nodes(mesh::MeshData, elem_id::Int) -> Vector{Int}

Get the node IDs that define a specific element.

# Arguments
- `mesh::MeshData`: Mesh data structure
- `elem_id::Int`: Element ID (1-based)

# Returns
- `Vector{Int}`: Vector of 4 node IDs

# Example
```julia
nodes = get_element_nodes(mesh, 5)
println("Element 5 nodes: ", nodes)
```
"""
function get_element_nodes(mesh::MeshData, elem_id::Int)
    return mesh.elements[elem_id, :]
end


"""
get_node_coordinates(mesh::MeshData, node_id::Int) -> Tuple{Float64, Float64}

Get the coordinates of a specific node.

# Arguments
- `mesh::MeshData`: Mesh data structure
- `node_id::Int`: Node ID (1-based)

# Returns
- `Tuple{Float64, Float64}`: (x, y) coordinates

# Example
```julia
x, y = get_node_coordinates(mesh, 10)
println("Node 10 at (", x, ", ", y, ")")
```
"""
function get_node_coordinates(mesh::MeshData, node_id::Int)
    return (mesh.coordinates[node_id, 1], mesh.coordinates[node_id, 2])
end


"""
get_element_material(mesh::MeshData, elem_id::Int) -> Union{Int, Nothing}

Get the material index assigned to an element.

# Arguments
- `mesh::MeshData`: Mesh data structure
- `elem_id::Int`: Element ID

# Returns
- `Union{Int, Nothing}`: Material index or nothing if not assigned
"""
function get_element_material(mesh::MeshData, elem_id::Int)
    return get(mesh.materials, elem_id, nothing)
end


"""
get_node_elements(mesh::MeshData, node_id::Int) -> Vector{Int}

Get all element IDs that contain a specific node.

# Arguments
- `mesh::MeshData`: Mesh data structure
- `node_id::Int`: Node ID (1-based)

# Returns
- `Vector{Int}`: Vector of element IDs that contain this node

# Example
```julia
elements = get_node_elements(mesh, 10)
println("Node 10 belongs to ", length(elements), " elements: ", elements)
```
"""
function get_node_elements(mesh::MeshData, node_id::Int)
    element_list = Int[]
    
    # Loop through all elements
    for elem_id in 1:mesh.num_elements
        # Check each node position directly without creating a slice
        if mesh.elements[elem_id, 1] == node_id || 
           mesh.elements[elem_id, 2] == node_id || 
           mesh.elements[elem_id, 3] == node_id || 
           mesh.elements[elem_id, 4] == node_id
           push!(element_list, elem_id)
        end
    end
    
    return element_list
end

"""
    has_pressure_bc(mesh::MeshData, node_id::Int) -> Bool

Check if a node has a constant pressure boundary condition.

# Arguments
- `mesh::MeshData`: The mesh data structure
- `node_id::Int`: The node ID to check

# Returns
- `Bool`: true if the node has a pressure BC, false otherwise
"""
function has_pressure_bc(mesh::MeshData, node_id::Int)
    return haskey(mesh.absolute_pressure_bc, node_id)
end


"""
    calculate_edge_outward_normal(mesh::MeshData, node_i::Int, node_j::Int) -> Tuple{Float64, Vector{Float64}}

Calculate the length and outward normal vector for a boundary edge between two nodes.

For 2D quadrilateral elements with counter-clockwise node ordering, the outward normal
is obtained by rotating the edge vector 90° clockwise.

# Arguments
- `mesh::MeshData`: The mesh data structure
- `node_i::Int`: First node ID of the edge
- `node_j::Int`: Second node ID of the edge

# Returns
- `Tuple{Float64, Vector{Float64}}`: Edge length l_e [m] and unit outward normal vector n̂ [dimensionless]

# Formula
Given edge vector: `edge_vec = x_j - x_i = [dx, dy]`
- Edge length: `l_e = ||edge_vec||`
- Outward normal (CCW ordering): `n̂ = [dy, -dx] / l_e`

# Example
```julia
mesh = read_mesh_file("mesh.mesh", materials)
l_e, n_hat = calculate_edge_outward_normal(mesh, 1, 2)
# Returns: (0.5, [0.0, -1.0]) for horizontal edge pointing right
```
"""
function calculate_edge_outward_normal(mesh::MeshData, node_i::Int, node_j::Int)
    # Get node coordinates
    x_i = mesh.coordinates[node_i, :]
    x_j = mesh.coordinates[node_j, :]
    
    # Calculate edge vector from node_i to node_j
    edge_vec = x_j - x_i  # [dx, dy]
    
    # Calculate edge length
    l_e = norm(edge_vec)
    
    if l_e < 1e-12
        error("Edge between nodes $node_i and $node_j has zero length (l_e = $l_e)")
    end
    
    # Calculate outward normal by rotating edge vector 90° clockwise
    # For counter-clockwise node ordering: [dx, dy] → [dy, -dx]
    outward_normal = [edge_vec[2], -edge_vec[1]]
    
    # Normalize to unit vector
    n_hat = outward_normal / l_e
    
    return l_e, n_hat
end


"""
    identify_boundary_edges(mesh::MeshData) -> Vector{Tuple{Int, Int, Int, Float64, Vector{Float64}}}

Identify all boundary edges where both nodes have pressure boundary conditions.

This function efficiently identifies boundary edges by only examining elements connected
to nodes with pressure BCs, avoiding unnecessary iteration over the entire mesh.

# Arguments
- `mesh::MeshData`: The mesh data structure

# Returns
- `Vector{Tuple{Int, Int, Int, Float64, Vector{Float64}}}`: Vector of boundary edge information tuples:
  - `element_id::Int`: Element containing the edge
  - `node_i::Int`: First node of the edge
  - `node_j::Int`: Second node of the edge (following element node ordering)
  - `l_e::Float64`: Edge length [m]
  - `n_hat::Vector{Float64}`: Unit outward normal vector [dimensionless]

# Algorithm
1. Extract all nodes with pressure BCs from `mesh.absolute_pressure_bc`
2. For each pressure BC node, find all connected elements
3. Check each element edge: if both nodes have pressure BC, it's a boundary edge
4. Calculate edge geometry (length and outward normal) and store
5. Remove duplicates (same edge may be found from both nodes)

# Example
```julia
mesh = read_mesh_file("mesh.mesh", materials)
boundary_edges = identify_boundary_edges(mesh)
println("Found ", length(boundary_edges), " boundary edges")

# Access edge information
for (elem_id, node_i, node_j, l_e, n_hat) in boundary_edges
    println("Element \$elem_id: edge \$node_i-\$node_j, length=\$l_e, normal=\$n_hat")
end
```
"""
function identify_boundary_edges(mesh::MeshData)
    # Use a Set to avoid duplicate edges
    boundary_edge_set = Set{Tuple{Int, Int, Int}}()  # (element_id, min_node, max_node)
    
    # Get all nodes with pressure boundary conditions
    pressure_bc_nodes = keys(mesh.absolute_pressure_bc)
    
    # For each pressure BC node, examine connected elements
    for node_id in pressure_bc_nodes
        # Get all elements containing this node
        connected_elements = get_node_elements(mesh, node_id)
        
        # Check each connected element for boundary edges
        for elem_id in connected_elements
            elem_nodes = mesh.elements[elem_id, :]
            
            # Check all 4 edges of the quadrilateral element
            # Edge 1: nodes[1] -> nodes[2]
            # Edge 2: nodes[2] -> nodes[3]
            # Edge 3: nodes[3] -> nodes[4]
            # Edge 4: nodes[4] -> nodes[1]
            for i in 1:4
                j = (i % 4) + 1  # Next node in sequence (wraps around)
                
                local_node_i = elem_nodes[i]
                local_node_j = elem_nodes[j]
                
                # Check if both nodes have pressure BC
                if has_pressure_bc(mesh, local_node_i) && has_pressure_bc(mesh, local_node_j)
                    # Store edge with normalized node ordering to avoid duplicates
                    # (element_id, min_node, max_node) ensures uniqueness
                    min_node = min(local_node_i, local_node_j)
                    max_node = max(local_node_i, local_node_j)
                    push!(boundary_edge_set, (elem_id, min_node, max_node))
                end
            end
        end
    end
    
    # Convert set to vector and calculate geometry for each edge
    boundary_edges = Vector{Tuple{Int, Int, Int, Float64, Vector{Float64}}}()
    
    for (elem_id, min_node, max_node) in boundary_edge_set
        # Determine correct node ordering from element connectivity
        # to preserve counter-clockwise orientation for outward normal
        elem_nodes = mesh.elements[elem_id, :]
        
        # Find the edge in the element's node sequence
        node_i, node_j = min_node, max_node
        for i in 1:4
            j = (i % 4) + 1
            if (elem_nodes[i] == min_node && elem_nodes[j] == max_node) ||
               (elem_nodes[i] == max_node && elem_nodes[j] == min_node)
                # Use the actual element ordering
                node_i = elem_nodes[i]
                node_j = elem_nodes[j]
                break
            end
        end
        
        # Calculate edge geometry
        l_e, n_hat = calculate_edge_outward_normal(mesh, node_i, node_j)
        
        # Store complete edge information
        push!(boundary_edges, (elem_id, node_i, node_j, l_e, n_hat))
    end
    
    return boundary_edges
end


# Export all public functions and types
export MeshData, read_mesh_file, get_element_nodes, get_node_coordinates
export get_element_material, get_node_elements, has_pressure_bc
export calculate_edge_outward_normal, identify_boundary_edges
