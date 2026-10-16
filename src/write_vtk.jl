"""
    write_vtk.jl

Legacy VTK output of ADSIM results (mesh, nodal scalar and vector fields) for ParaView.

Velocity vectors: 2 components (x, y) in 2D, 3 components (x, y, z) in 3D.

Requires read_mesh.jl to be included first (MeshData type).
"""

module WriteVTK

using ..ADSIMVersion: get_version

export write_vtk_file, write_vtk_series

"""
    write_vtk_file(filename, time_step, time, mesh, gas_names, fields::NamedTuple)

Write one VTK snapshot of the ADSIM state.

# Arguments
- `filename`: Base filename (without extension or step number)
- `time_step`: Output counter, for the `_NNNNNN` suffix
- `time`: Physical time of the snapshot
- `gas_names`: Gas species names, solver order
- `fields`: Nodal fields by name:
  `concentrations`, `total_concentration`, `absolute_pressure`,
  `concentration_rates`, `reaction_rates`, `lime_concentration`,
  `caco3_concentration`, `degree_of_carbonation`, `volumetric_binder_content`,
  `gas_seepage_velocity`, `temperature`, `temperature_rate`, `heat_flux`

Named instead of positional fields: most share the type `Vector{Float64}` or
`Matrix{Float64}`, and a swapped pair would be written under the wrong label without error.

# Note
Vector fields (`gas_seepage_velocity`, `heat_flux`) padded to three components (VTK
requirement) for any mesh dimensionality.
"""
function write_vtk_file(filename::String, time_step::Int, time::Float64, mesh,
                        gas_names::Vector{String}, fields::NamedTuple)
    # NamedTuple instead of sixteen positional arguments of the same types: a mismatched
    # name is an error at the call site instead of a mislabelled output field.
    concentrations            = fields.concentrations
    total_concentration       = fields.total_concentration
    absolute_pressure         = fields.absolute_pressure
    concentration_rates       = fields.concentration_rates
    reaction_rates            = fields.reaction_rates
    lime_concentration        = fields.lime_concentration
    caco3_concentration       = fields.caco3_concentration
    degree_of_carbonation     = fields.degree_of_carbonation
    volumetric_binder_content = fields.volumetric_binder_content
    gas_seepage_velocity      = fields.gas_seepage_velocity
    temperature               = fields.temperature
    temperature_rate          = fields.temperature_rate
    heat_flux                 = fields.heat_flux
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
        println(io, "ADSIM v$(get_version()) - Time Step $time_step, Time = $time")
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

        # Total heat flux q_d + q_a, Eqs. (fourier) and (heat_advection), projected
        # from the Gauss points. Zero with heat transport disabled.
        write_vtk_vector_field(io, "Heat_Flux", heat_flux, ndim)
    end
    
    println("      Wrote VTK file: $output_file")
    return nothing
end

"""
    write_vtk_mesh(io::IO, nodes::Matrix{Float64}, elements::Matrix{Int}, 
                   Nnodes::Int, Nelements::Int, ndim::Int, nodes_per_element::Int)

Write mesh geometry and topology to a legacy VTK file.
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

Return the VTK cell type identifier for the dimensionality and nodes per element.

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
    sanitize_value(x::Float64, threshold::Float64=1e-30)

Return 0.0 for |x| < threshold, otherwise x. Avoids ParaView parsing issues with
extreme exponents.
"""
function sanitize_value(x::Float64, threshold::Float64=1e-30)
    return abs(x) < threshold ? 0.0 : x
end

"""
    write_vtk_scalar_field(io::IO, field_name::String, data::Vector{Float64})

Write a scalar field to a legacy VTK file. |x| < 1e-30 written as zero.
"""
function write_vtk_scalar_field(io::IO, field_name::String, data::Vector{Float64})
    println(io, "SCALARS $field_name float 1")
    println(io, "LOOKUP_TABLE default")
    for i in 1:length(data)
        println(io, sanitize_value(data[i]))
    end
end

"""
    write_vtk_vector_field(io::IO, field_name::String, data::Matrix{Float64}, ndim::Int)

Write a vector field to a legacy VTK file. Three components required, z = 0 for 2D.
|x| < 1e-30 written as zero.
"""
function write_vtk_vector_field(io::IO, field_name::String, data::Matrix{Float64}, ndim::Int)
    nnodes = size(data, 1)
    
    println(io, "VECTORS $field_name float")
    for i in 1:nnodes
        if ndim == 2
            # Legacy VTK requires 3 components even for 2D
            println(io, "$(sanitize_value(data[i, 1])) $(sanitize_value(data[i, 2])) 0.0")
        else
            # 3D: x, y, z components
            println(io, "$(sanitize_value(data[i, 1])) $(sanitize_value(data[i, 2])) $(sanitize_value(data[i, 3]))")
        end
    end
end



end # module WriteVTK
