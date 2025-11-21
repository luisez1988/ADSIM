"""
    write_vtk.jl

Module for writing VTK output files for ADSIM visualization in ParaView.

This module provides functionality to export simulation results in VTK (Visualization Toolkit)
format, including mesh data with nodal scalar and vector fields.

Note: This module requires read_mesh.jl to be included first to access the MeshData type.
"""

module WriteVTK

export write_vtk_file, write_pvd_file

"""
    write_vtk_file(
        filename::String,
        time_step::Int,
        time::Float64,
        mesh,
        concentrations::Matrix{Float64},
        gas_names::Vector{String},
        total_concentration::Vector{Float64},
        absolute_pressure::Vector{Float64},
        concentration_rates::Matrix{Float64},
        reaction_rates::Vector{Float64},
        lime_concentration::Vector{Float64},
        co2_concentration::Vector{Float64},
        caco3_concentration::Vector{Float64},
        degree_of_carbonation::Vector{Float64},
        volumetric_binder_content::Vector{Float64},
        gas_seepage_velocity::Matrix{Float64},
        temperature::Vector{Float64},
        temperature_rate::Vector{Float64}
    )

Write VTK file for ADSIM simulation results at a specific time step.

# Arguments
- `filename::String`: Base filename (without extension) for the VTK file
- `time_step::Int`: Current time step number
- `time::Float64`: Physical time value
- `mesh`: Mesh data structure from read_mesh.jl (MeshData type)
- `concentrations::Matrix{Float64}`: Gas concentrations (Nnodes × Ngases)
- `gas_names::Vector{String}`: Names of gas species
- `total_concentration::Vector{Float64}`: Total concentration at nodes
- `absolute_pressure::Vector{Float64}`: Absolute pressure at nodes
- `concentration_rates::Matrix{Float64}`: Rate of concentration change (Nnodes × Ngases)
- `reaction_rates::Vector{Float64}`: Rate of reactions at nodes
- `lime_concentration::Vector{Float64}`: Lime concentration at nodes
- `co2_concentration::Vector{Float64}`: CO2 concentration at nodes
- `caco3_concentration::Vector{Float64}`: CaCO3 concentration at nodes
- `degree_of_carbonation::Vector{Float64}`: Degree of carbonation (DoC) at nodes
- `volumetric_binder_content::Vector{Float64}`: Volumetric binder content at nodes
- `gas_seepage_velocity::Matrix{Float64}`: Gas velocity vectors (Nnodes × ndim)
- `temperature::Vector{Float64}`: Temperature at nodes
- `temperature_rate::Vector{Float64}`: Rate of temperature change at nodes

# Returns
- `nothing`: Writes VTK file to disk

# Example
```julia
mesh = read_mesh_file("problem.mesh")
write_vtk_file(
    "output/result",
    10,
    0.5,
    mesh,
    concentrations,
    ["N2", "O2", "CO2"],
    total_conc,
    pressure,
    conc_rates,
    reaction_rates,
    lime_conc,
    co2_conc,
    caco3_conc,
    doc,
    binder,
    velocity,
    temp,
    temp_rate
)
```
"""
function write_vtk_file(
    filename::String,
    time_step::Int,
    time::Float64,
    mesh,
    concentrations::Matrix{Float64},
    gas_names::Vector{String},
    total_concentration::Vector{Float64},
    absolute_pressure::Vector{Float64},
    concentration_rates::Matrix{Float64},
    reaction_rates::Vector{Float64},
    lime_concentration::Vector{Float64},
    co2_concentration::Vector{Float64},
    caco3_concentration::Vector{Float64},
    degree_of_carbonation::Vector{Float64},
    volumetric_binder_content::Vector{Float64},
    gas_seepage_velocity::Matrix{Float64},
    temperature::Vector{Float64},
    temperature_rate::Vector{Float64}
)
    # Construct filename with time step
    vtk_filename = "$(filename)_$(lpad(time_step, 6, '0')).vtk"
    
    # Get mesh dimensions from MeshData structure
    Nnodes = mesh.num_nodes
    Nelements = mesh.num_elements
    ndim = size(mesh.coordinates, 2)
    nodes_per_element = size(mesh.elements, 2)
    Ngases = size(concentrations, 2)
    
    # Open file for writing
    open(vtk_filename, "w") do io
        # Write VTK header
        write_vtk_header(io, time, time_step)
        
        # Write mesh data using MeshData structure
        write_vtk_mesh(io, mesh.coordinates, mesh.elements, Nnodes, Nelements, ndim, nodes_per_element)
        
        # Write point data header
        println(io, "POINT_DATA $Nnodes")
        
        # Write scalar fields
        write_scalar_field(io, "Total_Concentration", total_concentration)
        write_scalar_field(io, "Absolute_Pressure", absolute_pressure)
        write_scalar_field(io, "Reaction_Rate", reaction_rates)
        write_scalar_field(io, "Lime_Concentration", lime_concentration)
        write_scalar_field(io, "CO2_Concentration", co2_concentration)
        write_scalar_field(io, "CaCO3_Concentration", caco3_concentration)
        write_scalar_field(io, "Degree_of_Carbonation", degree_of_carbonation)
        write_scalar_field(io, "Volumetric_Binder_Content", volumetric_binder_content)
        write_scalar_field(io, "Temperature", temperature)
        write_scalar_field(io, "Temperature_Rate", temperature_rate)
        
        # Write individual gas concentrations
        for i in 1:Ngases
            field_name = "$(gas_names[i])_Concentration"
            write_scalar_field(io, field_name, concentrations[:, i])
        end
        
        # Write individual gas concentration rates
        for i in 1:Ngases
            field_name = "$(gas_names[i])_Concentration_Rate"
            write_scalar_field(io, field_name, concentration_rates[:, i])
        end
        
        # Write vector field
        write_vector_field(io, "Gas_Seepage_Velocity", gas_seepage_velocity, ndim)
    end
    
    println("VTK file written: $vtk_filename")
    return nothing
end

"""
    write_vtk_header(io::IO, time::Float64, time_step::Int)

Write VTK file header with metadata.
"""
function write_vtk_header(io::IO, time::Float64, time_step::Int)
    println(io, "# vtk DataFile Version 3.0")
    println(io, "ADSIM Output - Time: $time, Step: $time_step")
    println(io, "ASCII")
    println(io, "DATASET UNSTRUCTURED_GRID")
end

"""
    write_vtk_mesh(io::IO, nodes::Matrix{Float64}, elements::Matrix{Int}, 
                   Nnodes::Int, Nelements::Int, ndim::Int, nodes_per_element::Int)

Write mesh geometry and topology to VTK file.
"""
function write_vtk_mesh(io::IO, nodes::Matrix{Float64}, elements::Matrix{Int}, 
                        Nnodes::Int, Nelements::Int, ndim::Int, nodes_per_element::Int)
    # Write points (nodes)
    println(io, "POINTS $Nnodes float")
    for i in 1:Nnodes
        if ndim == 2
            # For 2D problems, add z=0
            println(io, "$(nodes[i, 1]) $(nodes[i, 2]) 0.0")
        elseif ndim == 3
            println(io, "$(nodes[i, 1]) $(nodes[i, 2]) $(nodes[i, 3])")
        else
            # 1D case
            println(io, "$(nodes[i, 1]) 0.0 0.0")
        end
    end
    
    # Write cells (elements)
    # Calculate total size needed for connectivity data
    size_connectivity = Nelements * (nodes_per_element + 1)
    println(io, "CELLS $Nelements $size_connectivity")
    
    for i in 1:Nelements
        print(io, "$nodes_per_element")
        for j in 1:nodes_per_element
            # VTK uses 0-based indexing
            print(io, " $(elements[i, j] - 1)")
        end
        println(io)
    end
    
    # Write cell types
    println(io, "CELL_TYPES $Nelements")
    vtk_cell_type = get_vtk_cell_type(ndim, nodes_per_element)
    for i in 1:Nelements
        println(io, vtk_cell_type)
    end
end

"""
    get_vtk_cell_type(ndim::Int, nodes_per_element::Int)

Get VTK cell type identifier based on dimensionality and nodes per element.

VTK Cell Types:
- 1: Vertex
- 3: Line
- 5: Triangle
- 8: Pixel
- 9: Quad
- 10: Tetrahedron
- 12: Hexahedron
- 13: Wedge
- 14: Pyramid
"""
function get_vtk_cell_type(ndim::Int, nodes_per_element::Int)
    if ndim == 1
        if nodes_per_element == 2
            return 3  # Line
        else
            return 3  # Default to line for 1D
        end
    elseif ndim == 2
        if nodes_per_element == 3
            return 5  # Triangle
        elseif nodes_per_element == 4
            return 9  # Quad
        else
            error("Unsupported 2D element with $nodes_per_element nodes")
        end
    elseif ndim == 3
        if nodes_per_element == 4
            return 10  # Tetrahedron
        elseif nodes_per_element == 8
            return 12  # Hexahedron
        elseif nodes_per_element == 6
            return 13  # Wedge
        elseif nodes_per_element == 5
            return 14  # Pyramid
        else
            error("Unsupported 3D element with $nodes_per_element nodes")
        end
    else
        error("Unsupported dimensionality: $ndim")
    end
end

"""
    write_scalar_field(io::IO, field_name::String, data::Vector{Float64})

Write a scalar field to VTK file.
"""
function write_scalar_field(io::IO, field_name::String, data::Vector{Float64})
    Nnodes = length(data)
    println(io, "SCALARS $field_name float 1")
    println(io, "LOOKUP_TABLE default")
    for i in 1:Nnodes
        println(io, data[i])
    end
end

"""
    write_vector_field(io::IO, field_name::String, data::Matrix{Float64}, ndim::Int)

Write a vector field to VTK file.
"""
function write_vector_field(io::IO, field_name::String, data::Matrix{Float64}, ndim::Int)
    Nnodes = size(data, 1)
    println(io, "VECTORS $field_name float")
    
    for i in 1:Nnodes
        if ndim == 2
            # For 2D, add z=0 component
            println(io, "$(data[i, 1]) $(data[i, 2]) 0.0")
        elseif ndim == 3
            println(io, "$(data[i, 1]) $(data[i, 2]) $(data[i, 3])")
        else
            # For 1D, add y=0 and z=0 components
            println(io, "$(data[i, 1]) 0.0 0.0")
        end
    end
end

"""
    write_pvd_file(filename::String, time_steps::Vector{Int}, times::Vector{Float64})

Write a ParaView Data (PVD) collection file to enable time series visualization.

This function creates a master file that references all VTK files in the time series,
allowing ParaView to load and animate the entire simulation sequence.

# Arguments
- `filename::String`: Base filename (without extension) for the PVD file
- `time_steps::Vector{Int}`: Vector of time step numbers
- `times::Vector{Float64}`: Vector of corresponding physical time values

# Example
```julia
write_pvd_file("output/result", [0, 1, 2, 3], [0.0, 0.1, 0.2, 0.3])
```
"""
function write_pvd_file(filename::String, time_steps::Vector{Int}, times::Vector{Float64})
    pvd_filename = "$(filename).pvd"
    
    open(pvd_filename, "w") do io
        println(io, "<?xml version=\"1.0\"?>")
        println(io, "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">")
        println(io, "  <Collection>")
        
        for (step, time) in zip(time_steps, times)
            vtk_file = basename("$(filename)_$(lpad(step, 6, '0')).vtk")
            println(io, "    <DataSet timestep=\"$time\" file=\"$vtk_file\"/>")
        end
        
        println(io, "  </Collection>")
        println(io, "</VTKFile>")
    end
    
    println("PVD collection file written: $pvd_filename")
    return nothing
end

end # module WriteVTK
