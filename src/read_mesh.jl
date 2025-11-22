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
- `vacating_gas_bc::Dict{Int, Int}`: Vacating gas index for pressure BC (node_id => gas_index)
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
    vacating_gas_bc::Dict{Int, Int}
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
            Dict{Int, Int}(),
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
        mesh.absolute_pressure_bc[node_id] = pressure
        # Optional third field: vacating gas index
        if length(parts) >= 3
            gas_index = parse(Int, parts[3])
            mesh.vacating_gas_bc[node_id] = gas_index
        end
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
has_concentration_bc(mesh::MeshData, node_id::Int) -> Bool

Check if a node has concentration boundary conditions applied.

# Arguments
- `mesh::MeshData`: Mesh data structure
- `node_id::Int`: Node ID to check

# Returns
- `Bool`: true if node has concentration BC, false otherwise
"""
function has_concentration_bc(mesh::MeshData, node_id::Int)
    return haskey(mesh.concentration_bc, node_id)
end


"""
has_flow_bc(mesh::MeshData, node_id::Int) -> Bool

Check if a node has flow boundary conditions applied.

# Arguments
- `mesh::MeshData`: Mesh data structure
- `node_id::Int`: Node ID to check

# Returns
- `Bool`: true if node has flow BC, false otherwise
"""
function has_flow_bc(mesh::MeshData, node_id::Int)
    return haskey(mesh.uniform_flow_bc, node_id)
end


"""
has_pressure_bc(mesh::MeshData, node_id::Int) -> Bool

Check if a node has pressure boundary conditions applied.

# Arguments
- `mesh::MeshData`: Mesh data structure
- `node_id::Int`: Node ID to check

# Returns
- `Bool`: true if node has pressure BC, false otherwise
"""
function has_pressure_bc(mesh::MeshData, node_id::Int)
    return haskey(mesh.absolute_pressure_bc, node_id)
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


# Export all public functions and types
export MeshData, read_mesh_file, get_element_nodes, get_node_coordinates
export has_concentration_bc, has_flow_bc, has_pressure_bc, get_element_material
