"""
    write_vtk.jl

Module for writing VTK output files for ADSIM visualization in ParaView.

This module provides functionality to export simulation results in VTK (Visualization Toolkit)
format, including mesh data with nodal scalar and vector fields.

For 2D simulations, velocity vectors have 2 components (x, y).
For 3D simulations, velocity vectors have 3 components (x, y, z).

Note: This module requires read_mesh.jl to be included first to access the MeshData type.
"""

module WriteVTK

export write_vtk_file

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
- `nothing`: Writes VTU file to disk

# Output
- Writes a VTU (XML VTK Unstructured Grid) file compatible with ParaView
- For 2D simulations, velocity has 2 components (x, y)
- For 3D simulations, velocity has 3 components (x, y, z)

# Example
```julia
mesh = read_mesh_file("problem.mesh")
write_vtk_file(
    "output/simulation",
    0,
    0.0,
    mesh,
    concentrations,
    gas_names,
    total_concentration,
    absolute_pressure,
    concentration_rates,
    reaction_rates,
    lime_concentration,
    co2_concentration,
    caco3_concentration,
    degree_of_carbonation,
    volumetric_binder_content,
    gas_seepage_velocity,
    temperature,
    temperature_rate
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
    # Create output filename
    output_file = "$(filename)_$(lpad(time_step, 6, '0')).vtk"
    
    # Determine mesh dimensionality
    ndim = size(mesh.coordinates, 2)
    Nnodes = mesh.num_nodes
    Nelements = mesh.num_elements
    nodes_per_element = size(mesh.elements, 2)
    
    open(output_file, "w") do io
        # Write VTK header
        println(io, "# vtk DataFile Version 3.0")
        println(io, "ADSIM Simulation Results - Time Step $time_step, Time = $time")
        println(io, "ASCII")
        println(io, "DATASET UNSTRUCTURED_GRID")
        
        # Write mesh geometry and topology
        write_vtk_mesh(io, mesh.coordinates, mesh.elements, Nnodes, Nelements, ndim, nodes_per_element)
        
        # Write point data section
        println(io, "POINT_DATA $Nnodes")
        
        # Write scalar fields
        write_vtk_scalar_field(io, "Total_Concentration", total_concentration)
        write_vtk_scalar_field(io, "Absolute_Pressure", absolute_pressure)
        write_vtk_scalar_field(io, "Reaction_Rate", reaction_rates)
        write_vtk_scalar_field(io, "Lime_Concentration", lime_concentration)
        write_vtk_scalar_field(io, "CO2_Concentration", co2_concentration)
        write_vtk_scalar_field(io, "CaCO3_Concentration", caco3_concentration)
        write_vtk_scalar_field(io, "Degree_of_Carbonation", degree_of_carbonation)
        write_vtk_scalar_field(io, "Volumetric_Binder_Content", volumetric_binder_content)
        write_vtk_scalar_field(io, "Temperature", temperature)
        write_vtk_scalar_field(io, "Temperature_Rate", temperature_rate)
        
        # Write individual gas concentrations
        for i in 1:length(gas_names)
            field_name = "$(gas_names[i])_Concentration"
            write_vtk_scalar_field(io, field_name, concentrations[:, i])
        end
        
        # Write concentration rates
        for i in 1:length(gas_names)
            field_name = "$(gas_names[i])_Concentration_Rate"
            write_vtk_scalar_field(io, field_name, concentration_rates[:, i])
        end
        
        # Write velocity field (2D or 3D based on mesh)
        write_vtk_vector_field(io, "Gas_Seepage_Velocity", gas_seepage_velocity, ndim)
    end
    
    println("Wrote VTK file: $output_file")
    return nothing
end

"""
    write_vtk_mesh(io::IO, nodes::Matrix{Float64}, elements::Matrix{Int}, 
                   Nnodes::Int, Nelements::Int, ndim::Int, nodes_per_element::Int)

Write mesh geometry and topology to legacy VTK file.
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
        end
    end
    
    # Write cells (elements)
    vtk_cell_type = get_vtk_cell_type(ndim, nodes_per_element)
    total_connectivity_size = Nelements * (nodes_per_element + 1)  # +1 for count
    
    println(io, "CELLS $Nelements $total_connectivity_size")
    for i in 1:Nelements
        print(io, "$nodes_per_element")
        for j in 1:nodes_per_element
            # VTK uses 0-based indexing
            print(io, " $(elements[i, j] - 1)")
        end
        println(io)
    end
    
    # Cell types
    println(io, "CELL_TYPES $Nelements")
    for i in 1:Nelements
        println(io, "$vtk_cell_type")
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
    if  ndim == 2
        if nodes_per_element == 3
            return 5  # Triangle
        elseif nodes_per_element == 4 #only available option for now
            return 9  # Quad
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
        end
    end
end

"""
    write_vtk_scalar_field(io::IO, field_name::String, data::Vector{Float64})

Write a scalar field to legacy VTK file.
"""
function write_vtk_scalar_field(io::IO, field_name::String, data::Vector{Float64})
    println(io, "SCALARS $field_name float 1")
    println(io, "LOOKUP_TABLE default")
    for i in 1:length(data)
        println(io, data[i])
    end
end

"""
    write_vtk_vector_field(io::IO, field_name::String, data::Matrix{Float64}, ndim::Int)

Write a vector field to legacy VTK file.
Legacy VTK format requires 3 components for vectors, so z=0 is added for 2D.
"""
function write_vtk_vector_field(io::IO, field_name::String, data::Matrix{Float64}, ndim::Int)
    nnodes = size(data, 1)
    
    println(io, "VECTORS $field_name float")
    for i in 1:nnodes
        if ndim == 2
            # Legacy VTK requires 3 components even for 2D
            println(io, "$(data[i, 1]) $(data[i, 2]) 0.0")
        else
            # 3D: write x, y, and z components
            println(io, "$(data[i, 1]) $(data[i, 2]) $(data[i, 3])")
        end
    end
end



end # module WriteVTK
