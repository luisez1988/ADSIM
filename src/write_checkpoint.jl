#______________________________________________________
# ADSIM: Advection-Diffusion for Soil Improvement and 
# Modification
# Author: Luis Zambrano-Cruzatty
#______________________________________________________

#______________________________________________________
# Checkpoint Writer Module
# Saves simulation state to binary file for multi-stage calculations
#______________________________________________________

using JLD2

"""
    write_checkpoint(project_name::String, current_stage::Int, current_time::Float64, 
                    output_counter::Int, next_output_time::Float64)

Write a binary checkpoint file containing all essential state variables needed to 
continue the simulation in a subsequent stage.

# Essential State Variables Saved:
- Primary state: C_g, P, T, v (gas concentration, pressure, temperature, velocity)
- Reactive species: C_lime, C_caco3, binder_content, degree_of_carbonation
- Material properties: C_lime_residual (per soil type)
- Dimensions: NDim, Nnodes, Nelements, NSoils, NGases
- Time tracking: current_time, output_counter, next_output_time

# Variables NOT saved (recalculated or reapplied from input files):
- Time derivatives: dC_g_dt, dT_dt, dC_lime_dt (recomputed each time step)
- Flow terms: q_boundary, q_advection, q_diffusion, etc. (recomputed/reapplied)
- Boundary data: P_boundary, λ_bc, boundary_node_influences (from mesh file)
- Derived values: Caco3_max (recalculated from initial material properties)

# Arguments
- `project_name::String`: Base name of the project (e.g., "Test")
- `current_stage::Int`: Current stage number (checkpoint will be for this stage)
- `current_time::Float64`: Current simulation time
- `output_counter::Int`: Current output file counter (for VTK numbering)
- `next_output_time::Float64`: Next scheduled output time

# Output
Creates file: `output/{project_name}_stage{current_stage}.jld2`

# Returns
- `String`: Path to the created checkpoint file
"""
function write_checkpoint(project_name::String, current_stage::Int, current_time::Float64, output_counter::Int, next_output_time::Float64)
    
    # Construct checkpoint file path
    output_dir = "output"
    checkpoint_file = joinpath(output_dir, "$(project_name)_stage$(current_stage).jld2")
    
    # Ensure output directory exists
    if !isdir(output_dir)
        mkdir(output_dir)
    end
    
    # Save all essential state variables using JLD2
    jldsave(checkpoint_file;
        # Primary state variables
        C_g = C_g,
        P = P,
        T = T,
        v = v,
        
        # Reactive species state
        C_lime = C_lime,
        C_caco3 = C_caco3,
        binder_content = binder_content,
        degree_of_carbonation = degree_of_carbonation,
        
        # Material properties (per soil type)
        C_lime_residual = C_lime_residual,
        
        # Dimension parameters
        NDim = NDim,
        Nnodes = Nnodes,
        Nelements = Nelements,
        NSoils = NSoils,
        NGases = NGases,
        
        # Time tracking
        current_time = current_time,
        output_counter = output_counter,
        next_output_time = next_output_time
    )
    
    return checkpoint_file
end


"""
    get_checkpoint_file_size(checkpoint_file::String)

Get human-readable file size of checkpoint file.

# Arguments
- `checkpoint_file::String`: Path to checkpoint file

# Returns
- `String`: Formatted file size (e.g., "2.5 MB", "1.3 KB")
"""
function get_checkpoint_file_size(checkpoint_file::String)
    if !isfile(checkpoint_file)
        return "0 B"
    end
    
    size_bytes = filesize(checkpoint_file)
    
    if size_bytes < 1024
        return "$(size_bytes) B"
    elseif size_bytes < 1024^2
        return @sprintf("%.1f KB", size_bytes / 1024)
    elseif size_bytes < 1024^3
        return @sprintf("%.1f MB", size_bytes / 1024^2)
    else
        return @sprintf("%.1f GB", size_bytes / 1024^3)
    end
end
