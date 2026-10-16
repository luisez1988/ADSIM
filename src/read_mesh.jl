#------------------------------------------------------------------------------
# ADSIM Mesh Reader Module
# Reading and parsing of .mesh files for ADSIM FEM calculations
#------------------------------------------------------------------------------

"""
MeshData

Mesh data with boundary and initial conditions.

# Fields
- `num_nodes::Int`: Total number of nodes in the mesh
- `num_elements::Int`: Total number of elements in the mesh
- `coordinates::Matrix{Float64}`: Nodal coordinates (num_nodes × 2)
- `elements::Matrix{Int}`: Element connectivity (num_elements × 4)
- `concentration_bc::Dict{Int, Vector{Float64}}`: Concentration BC (node_id => [gas1, gas2, ...])
- `uniform_flow_bc::Dict{Int, Vector{Float64}}`: Uniform flow BC (node_id => [gas1, gas2, ...])
- `absolute_pressure_bc::Dict{Int, Float64}`: Absolute pressure BC (node_id => pressure)
- `partial_pressure_bc::Dict{Int, Vector{Float64}}`: Partial pressure BC (node_id => [P_gas1, P_gas2, ...])
# - `vacating_gas_bc::Dict{Int, Int}`: Vacating gas index for pressure BC (node_id => gas_index)
- `initial_concentrations::Dict{Int, Vector{Float64}}`: Initial concentrations (elem_id => [gas1, gas2, ...])
- `initial_temperature::Dict{Int, Float64}`: Initial temperature (elem_id => temperature)
- `temperature_bc::Dict{Int, Float64}`: Prescribed temperature (node_id => T), Eq. (temperature_BC)
- `convective_heat_bc::Dict{Int, Vector{Float64}}`: Convective heat BC (node_id => [h_c, T_inf]),
  Eq. (conduction_flux_BC). h_c = 0 or unset boundary: adiabatic (natural condition of
  the weak form). The two thermal boundary sets must be disjoint.
- `materials::Dict{Int, Int}`: Material assignment (elem_id => material_index)

Time-dependent boundaries: parallel `*_tf` fields with the time function id of each
value, 0 = constant. Stored separately instead of as extra columns, since
`parse_concentration_bc!` reads every remaining field as a gas value without knowing
the number of gases. Without `*_tf` blocks every id defaults to 0 (constant boundaries).

- `concentration_bc_tf::Dict{Int, Vector{Int}}`: Curve ids (node_id => [gas1, gas2, ...])
- `uniform_flow_bc_tf::Dict{Int, Vector{Int}}`: Curve ids (node_id => [gas1, gas2, ...])
- `absolute_pressure_tf::Dict{Int, Int}`: Curve id (node_id => tf_id)
- `temperature_bc_tf::Dict{Int, Int}`: Curve id (node_id => tf_id)
- `partial_pressure_bc_tf::Dict{Int, Vector{Int}}`: Curve ids (node_id => [gas1, gas2, ...]),
  for `partial_pressure_bc`. The solver looks up the curve inline each step instead of a
  pre-filtered `transient_*` list, since partial_pressure_bc is re-applied every step
  (follows T via the ideal gas law).
"""

mutable struct MeshData
    num_nodes::Int
    num_elements::Int
    coordinates::Matrix{Float64}
    elements::Matrix{Int}
    concentration_bc::Dict{Int, Vector{Float64}}
    uniform_flow_bc::Dict{Int, Vector{Float64}}
    absolute_pressure_bc::Dict{Int, Float64}
    partial_pressure_bc::Dict{Int, Vector{Float64}}
    initial_concentrations::Dict{Int, Vector{Float64}}
    initial_temperature::Dict{Int, Float64}
    temperature_bc::Dict{Int, Float64}
    convective_heat_bc::Dict{Int, Vector{Float64}}
    materials::Dict{Int, Int}

    # Time function ids, parallel to the boundary fields above (0 = constant)
    concentration_bc_tf::Dict{Int, Vector{Int}}
    uniform_flow_bc_tf::Dict{Int, Vector{Int}}
    absolute_pressure_tf::Dict{Int, Int}
    temperature_bc_tf::Dict{Int, Int}
    partial_pressure_bc_tf::Dict{Int, Vector{Int}}

    function MeshData()
        new(0, 0,
            zeros(Float64, 0, 0),
            zeros(Int, 0, 0),
            Dict{Int, Vector{Float64}}(),
            Dict{Int, Vector{Float64}}(),
            Dict{Int, Float64}(),
            Dict{Int, Vector{Float64}}(),
            Dict{Int, Vector{Float64}}(),
            Dict{Int, Float64}(),
            Dict{Int, Float64}(),
            Dict{Int, Vector{Float64}}(),
            Dict{Int, Int}(),
            Dict{Int, Vector{Int}}(),
            Dict{Int, Vector{Int}}(),
            Dict{Int, Int}(),
            Dict{Int, Int}(),
            Dict{Int, Vector{Int}}())
    end
end


"""
BoundaryNodeInfluence

Boundary node influence lengths.

# Fields
- `node_influences::Dict{Int, Float64}`: node_id => total influence length [m]

Influence length: portion of boundary edges assigned to each node. Edges with the
condition on both nodes give half the edge length to each node.
"""
struct BoundaryNodeInfluence
    node_influences::Dict{Int, Float64}
    
    function BoundaryNodeInfluence()
        new(Dict{Int, Float64}())
    end
end


"""
read_mesh_file(filename::String) -> MeshData

Read an ADSIM .mesh file into a `MeshData` structure (geometry and physical data).

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
                
            # Parse partial pressure boundary conditions
            elseif line == "partial_pressure_bc"
                line_idx = parse_partial_pressure_bc!(mesh, lines, line_idx + 1)
                
            # Parse initial concentrations
            elseif line == "initial_concentrations"
                line_idx = parse_initial_concentrations!(mesh, lines, line_idx + 1)
                
            # Parse initial temperature
            elseif line == "initial_temperature"
                line_idx = parse_initial_temperature!(mesh, lines, line_idx + 1)
                
            # Parse prescribed temperature (thermal Dirichlet)
            elseif line == "temperature_bc"
                line_idx = parse_temperature_bc!(mesh, lines, line_idx + 1)

            # Parse convective heat boundary (thermal Robin)
            elseif line == "convective_heat_bc"
                line_idx = parse_convective_heat_bc!(mesh, lines, line_idx + 1)

            # Time function ids of the boundary conditions above.
            # Optional blocks; without them boundaries are constant.
            elseif line == "concentration_bc_tf"
                line_idx = parse_bc_tf_vector!(mesh.concentration_bc_tf, lines,
                                               line_idx + 1, "concentration_bc_tf")

            elseif line == "uniform_flow_bc_tf"
                line_idx = parse_bc_tf_vector!(mesh.uniform_flow_bc_tf, lines,
                                               line_idx + 1, "uniform_flow_bc_tf")

            elseif line == "absolute_pressure_tf"
                line_idx = parse_bc_tf_scalar!(mesh.absolute_pressure_tf, lines,
                                               line_idx + 1, "absolute_pressure_tf")

            elseif line == "temperature_bc_tf"
                line_idx = parse_bc_tf_scalar!(mesh.temperature_bc_tf, lines,
                                               line_idx + 1, "temperature_bc_tf")

            elseif line == "partial_pressure_bc_tf"
                line_idx = parse_bc_tf_vector!(mesh.partial_pressure_bc_tf, lines,
                                               line_idx + 1, "partial_pressure_bc_tf")

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

Mutates `mesh`. Returns the next line index.
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

Mutates `mesh`. Returns the next line index.
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

Mutates `mesh`. Returns the next line index.
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

Mutates `mesh`. Returns the next line index.
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

Mutates `mesh`. Returns the next line index.
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

Mutates `mesh`. Returns the next line index.
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
    parse_temperature_bc!(mesh::MeshData, lines::Vector{String}, line_idx::Int)

Parse the `temperature_bc` block: a prescribed nodal temperature, Eq. (temperature_BC).

Same format as the other nodal blocks: counter line, `node_id value` lines,
`end temperature_bc`.

Mutates `mesh`. Returns the index of the first line after the block.
"""
function parse_temperature_bc!(mesh::MeshData, lines::Vector{String}, line_idx::Int)
    counter = parse(Int, strip(lines[line_idx]))
    line_idx += 1

    for _ in 1:counter
        parts = split(strip(lines[line_idx]))
        node_id = parse(Int, parts[1])
        mesh.temperature_bc[node_id] = parse(Float64, parts[2])
        line_idx += 1
    end

    if strip(lines[line_idx]) == "end temperature_bc"
        line_idx += 1
    end

    return line_idx
end


"""
    parse_convective_heat_bc!(mesh::MeshData, lines::Vector{String}, line_idx::Int)

Parse the `convective_heat_bc` block, Eq. (conduction_flux_BC):

    q_d . n = h_c (T - T_inf)

Data lines: `node_id h_c T_inf`. `h_c = 0` or omitted node: adiabatic (natural
condition of the weak form).

Mutates `mesh`. Returns the index of the first line after the block.
"""
function parse_convective_heat_bc!(mesh::MeshData, lines::Vector{String}, line_idx::Int)
    counter = parse(Int, strip(lines[line_idx]))
    line_idx += 1

    for _ in 1:counter
        parts = split(strip(lines[line_idx]))
        node_id = parse(Int, parts[1])
        h_c = parse(Float64, parts[2])
        T_inf = length(parts) >= 3 ? parse(Float64, parts[3]) : 0.0
        mesh.convective_heat_bc[node_id] = [h_c, T_inf]
        line_idx += 1
    end

    if strip(lines[line_idx]) == "end convective_heat_bc"
        line_idx += 1
    end

    return line_idx
end


"""
parse_bc_tf_vector!(target::Dict{Int, Vector{Int}}, lines, line_idx::Int, block::String) -> Int

Parse a block of per-gas time function ids, one line per node:

    concentration_bc_tf
    2
    40 1 0
    42 1 0
    end concentration_bc_tf

Shared by the concentration and flow blocks (same layout). Ids: 1-based references to
the `[[time_function]]` section of the calculation file; 0 = constant value from the
matching value block. Mutates `target`, returns the next line index.
"""
function parse_bc_tf_vector!(target::Dict{Int, Vector{Int}}, lines::Vector{String},
                             line_idx::Int, block::String)
    counter = parse(Int, strip(lines[line_idx]))
    line_idx += 1

    for _ in 1:counter
        parts = split(strip(lines[line_idx]))
        node_id = parse(Int, parts[1])
        target[node_id] = [parse(Int, parts[i]) for i in 2:length(parts)]
        line_idx += 1
    end

    if strip(lines[line_idx]) == "end $(block)"
        line_idx += 1
    end

    return line_idx
end

"""
parse_bc_tf_scalar!(target::Dict{Int, Int}, lines, line_idx::Int, block::String) -> Int

Parse a block of single time function ids, one line per node:

    temperature_bc_tf
    2
    1 2
    4 2
    end temperature_bc_tf

Scalar version of `parse_bc_tf_vector!`, for absolute pressure and prescribed
temperature boundaries. Mutates `target`, returns the next line index.
"""
function parse_bc_tf_scalar!(target::Dict{Int, Int}, lines::Vector{String},
                             line_idx::Int, block::String)
    counter = parse(Int, strip(lines[line_idx]))
    line_idx += 1

    for _ in 1:counter
        parts = split(strip(lines[line_idx]))
        target[parse(Int, parts[1])] = parse(Int, parts[2])
        line_idx += 1
    end

    if strip(lines[line_idx]) == "end $(block)"
        line_idx += 1
    end

    return line_idx
end


"""
parse_partial_pressure_bc!(mesh::MeshData, lines::Vector{String}, line_idx::Int) -> Int

Parse partial pressure boundary conditions from the mesh file.

Mutates `mesh`. Returns the next line index.
"""
function parse_partial_pressure_bc!(mesh::MeshData, lines::Vector{String}, line_idx::Int)
    # Read counter
    counter = parse(Int, strip(lines[line_idx]))
    line_idx += 1
    
    # Read boundary conditions
    for i in 1:counter
        line = strip(lines[line_idx])
        parts = split(line)
        node_id = parse(Int, parts[1])
        partial_pressures = [parse(Float64, parts[i]) for i in 2:length(parts)]
        mesh.partial_pressure_bc[node_id] = partial_pressures
        line_idx += 1
    end
    
    # Skip end marker
    if strip(lines[line_idx]) == "end partial_pressure_bc"
        line_idx += 1
    end
    
    return line_idx
end


"""
parse_initial_concentrations!(mesh::MeshData, lines::Vector{String}, line_idx::Int) -> Int

Parse initial gas concentrations from the mesh file.

Mutates `mesh`. Returns the next line index.
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

Mutates `mesh`. Returns the next line index.
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

Mutates `mesh`. Returns the next line index.
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

Return the 4 node IDs of an element (1-based `elem_id`).

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

Return the (x, y) coordinates of a node (1-based `node_id`).

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

Return the material index of an element; `nothing` if unassigned.
"""
function get_element_material(mesh::MeshData, elem_id::Int)
    return get(mesh.materials, elem_id, nothing)
end


"""
get_node_elements(mesh::MeshData, node_id::Int) -> Vector{Int}

Return the IDs of all elements containing a node (1-based `node_id`).

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
        # Direct comparison, no slice allocation
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

Return `true` if the node has an absolute pressure boundary condition.
"""
function has_pressure_bc(mesh::MeshData, node_id::Int)
    return haskey(mesh.absolute_pressure_bc, node_id)
end


"""
    calculate_edge_outward_normal(mesh::MeshData, node_i::Int, node_j::Int) -> Tuple{Float64, Vector{Float64}}

Compute the length l_e [m] and unit outward normal n̂ of the boundary edge `node_i -> node_j`.

Counter-clockwise node ordering: outward normal = edge vector rotated 90° clockwise.

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

Identify all boundary edges with the condition on both nodes (default: absolute pressure BC).

Only elements connected to BC nodes are examined, no full-mesh iteration.

# Arguments
- `bc_nodes`: Node set of the condition; defaults to the absolute-pressure boundary

# Returns
- Vector of edge tuples:
  - `element_id::Int`: Element containing the edge
  - `node_i::Int`: First node of the edge
  - `node_j::Int`: Second node of the edge (following element node ordering)
  - `l_e::Float64`: Edge length [m]
  - `n_hat::Vector{Float64}`: Unit outward normal vector [dimensionless]

# Algorithm
1. For each BC node, find all connected elements
2. Element edge with both nodes in the BC set: boundary edge
3. Compute edge geometry (length and outward normal)
4. Remove duplicates (same edge found from both nodes)

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
function identify_boundary_edges(mesh::MeshData, bc_nodes = keys(mesh.absolute_pressure_bc))
    # Use a Set to avoid duplicate edges
    boundary_edge_set = Set{Tuple{Int, Int, Int}}()  # (element_id, min_node, max_node)

    # Node set as a parameter (default: absolute-pressure boundary), reused by the thermal
    # boundary. An edge belongs to a condition when BOTH nodes have it.
    bc_node_set = Set(bc_nodes)

    # For each BC node, examine connected elements
    for node_id in bc_node_set
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
                
                # Both nodes in the BC set
                if (local_node_i in bc_node_set) && (local_node_j in bc_node_set)
                    # Normalized ordering (element_id, min_node, max_node), no duplicates
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


"""
    get_boundary_node_influences(mesh::MeshData) -> BoundaryNodeInfluence

Compute the influence length of each node on the boundary edges of a condition.

Boundary edge (both nodes in the BC set): l_e/2 to each node, accumulated over edges.

# Arguments
- `bc_nodes`: Node set of the condition; defaults to the absolute-pressure boundary
- `axisymmetric`: Weight each edge by its midpoint radius, giving the ring measure
  `l_e * r̄_e` of Eq. (ax_nodal_Lagrange_System1) instead of the edge length. `2π`
  cancelled, as in the rest of the axisymmetric formulation.

# Returns
- `BoundaryNodeInfluence`: node_id => influence length [m] (plane strain) or ring
  measure [m²] (axisymmetric)

# Algorithm
1. For each BC node, find all connected elements
2. Element edge with both nodes in the BC set: boundary edge
3. Half the edge length to each node
4. Accumulate over nodes on multiple edges

# Axis of symmetry
Edges on the axis have `r̄_e = 0` and contribute nothing, without a special case:
Eq. (ax_external_flux) vanishes there (zero ring area), i.e. no flux through the axis
without an imposed Neumann condition.

# Example
```julia
mesh = read_mesh_file("mesh.mesh")
node_influences = get_boundary_node_influences(mesh)

# Access influence length for a specific node
for (node_id, influence_length) in node_influences.node_influences
    println("Node \$node_id has influence length: \$influence_length m")
end
```
"""
function get_boundary_node_influences(mesh::MeshData, bc_nodes = keys(mesh.absolute_pressure_bc);
                                      axisymmetric::Bool = false)
    # Initialize the result structure
    influences = BoundaryNodeInfluence()

    # Use a Set to avoid processing duplicate edges
    processed_edges = Set{Tuple{Int, Int}}()  # (min_node, max_node)

    # Node set as a parameter (as in identify_boundary_edges), reused by the convective
    # heat boundary for its tributary edge lengths
    bc_node_set = Set(bc_nodes)

    # For each BC node, examine connected elements
    for node_id in bc_node_set
        # Get all elements containing this node
        connected_elements = get_node_elements(mesh, node_id)
        
        # Check each connected element for boundary edges
        for elem_id in connected_elements
            elem_nodes = mesh.elements[elem_id, :]
            
            # Check all 4 edges of the quadrilateral element
            for i in 1:4
                j = (i % 4) + 1  # Next node in sequence (wraps around)
                
                local_node_i = elem_nodes[i]
                local_node_j = elem_nodes[j]
                
                # Both nodes in the BC set
                if (local_node_i in bc_node_set) && (local_node_j in bc_node_set)
                    # Normalized edge identifier, no duplicates
                    min_node = min(local_node_i, local_node_j)
                    max_node = max(local_node_i, local_node_j)
                    edge_key = (min_node, max_node)
                    
                    # Only process each edge once
                    if !(edge_key in processed_edges)
                        push!(processed_edges, edge_key)
                        
                        # Calculate edge length
                        l_e, _ = calculate_edge_outward_normal(mesh, local_node_i, local_node_j)

                        # Axisymmetric boundary measure: ring area l_e * r̄_e, r̄_e = midpoint
                        # radius (Eq. ax_nodal_Lagrange_System1). r = first coordinate.
                        if axisymmetric
                            r_bar = 0.5 * (mesh.coordinates[local_node_i, 1] +
                                           mesh.coordinates[local_node_j, 1])
                            l_e *= r_bar
                        end

                        # Distribute half the edge measure to each node
                        half_length = l_e / 2.0
                        
                        # Accumulate influence length for each node
                        influences.node_influences[local_node_i] = 
                            get(influences.node_influences, local_node_i, 0.0) + half_length
                        influences.node_influences[local_node_j] = 
                            get(influences.node_influences, local_node_j, 0.0) + half_length
                    end
                end
            end
        end
    end
    
    return influences
end


# Export all public functions and types
export MeshData, read_mesh_file, get_element_nodes, get_node_coordinates
export get_element_material, get_node_elements, has_pressure_bc
export calculate_edge_outward_normal, BoundaryNodeInfluence, get_boundary_node_influences
export identify_boundary_edges
