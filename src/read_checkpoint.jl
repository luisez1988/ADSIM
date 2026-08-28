#______________________________________________________
# ADSIM: Advection-Diffusion for Soil Improvement and 
# Modification
# Author: Luis Zambrano-Cruzatty
#______________________________________________________

#______________________________________________________
# Checkpoint Reader Module
# Loads simulation state from binary file for multi-stage calculations
#______________________________________________________

using JLD2

"""
    find_latest_checkpoint(project_name::String, output_dir::String="output")

Find the most recent checkpoint file for a given project.

# Arguments
- `project_name::String`: Base name of the project (e.g., "Test")
- `output_dir::String`: Directory containing checkpoint files (default: "output")

# Returns
- `(checkpoint_file::String, stage_number::Int)`: Path to checkpoint and its stage number
- `(nothing, 0)`: If no checkpoint files are found
"""
function find_latest_checkpoint(project_name::String, output_dir::String="output")
    
    # Check if output directory exists
    if !isdir(output_dir)
        return (nothing, 0)
    end
    
    # Look for checkpoint files matching pattern: {project_name}_stage*.jld2
    all_files = readdir(output_dir)
    checkpoint_files = filter(f -> occursin(r"_stage\d+\.jld2$", f) && startswith(f, project_name), all_files)
    
    if isempty(checkpoint_files)
        return (nothing, 0)
    end
    
    # Extract stage numbers from checkpoint filenames
    stage_data = Tuple{String, Int}[]
    for filename in checkpoint_files
        # Extract number between "stage" and ".jld2"
        m = match(r"_stage(\d+)\.jld2$", filename)
        if m !== nothing
            stage_num = parse(Int, m.captures[1])
            filepath = joinpath(output_dir, filename)
            push!(stage_data, (filepath, stage_num))
        end
    end
    
    if isempty(stage_data)
        return (nothing, 0)
    end
    
    # Return the checkpoint with the highest stage number
    max_stage = 0
    latest_file = ""
    for (filepath, stage_num) in stage_data
        if stage_num > max_stage
            max_stage = stage_num
            latest_file = filepath
        end
    end
    return (latest_file, max_stage)
end


"""
    validate_checkpoint_dimensions(checkpoint_data::Dict, mesh, materials)

Validate that checkpoint dimensions match current mesh and materials.

# Arguments
- `checkpoint_data::Dict`: Loaded checkpoint data
- `mesh`: Current mesh data structure
- `materials`: Current material data structure

# Returns
- `(valid::Bool, message::String)`: Validation result and error/warning message
"""
function validate_checkpoint_dimensions(checkpoint_data::Dict, mesh, materials)
    
    errors = String[]
    
    # Check node count
    if checkpoint_data["Nnodes"] != mesh.num_nodes
        push!(errors, "Node count mismatch: checkpoint has $(checkpoint_data["Nnodes"]), mesh has $(mesh.num_nodes)")
    end
    
    # Check element count
    if checkpoint_data["Nelements"] != mesh.num_elements
        push!(errors, "Element count mismatch: checkpoint has $(checkpoint_data["Nelements"]), mesh has $(mesh.num_elements)")
    end
    
    # Check number of gases
    NGases_current = length(materials.gas_dictionary)
    if checkpoint_data["NGases"] != NGases_current
        push!(errors, "Gas count mismatch: checkpoint has $(checkpoint_data["NGases"]), materials have $(NGases_current)")
    end
    
    # Check number of soils
    NSoils_current = length(materials.soil_dictionary)
    if checkpoint_data["NSoils"] != NSoils_current
        push!(errors, "Soil count mismatch: checkpoint has $(checkpoint_data["NSoils"]), materials have $(NSoils_current)")
    end
    
    # Check spatial dimensions
    if checkpoint_data["NDim"] != 2
        push!(errors, "Dimension mismatch: checkpoint has $(checkpoint_data["NDim"])D, expected 2D")
    end

    # Check the geometric formulation. Plane strain and axisymmetric are different physics
    # on the same mesh, so resuming across them would continue from a state the current
    # equations never produced. Checkpoints written before this field existed are plane by
    # construction, since that was the only formulation available.
    ckpt_axi = get(checkpoint_data, "axisymmetric", false)
    current_axi = Axisymmetric
    if ckpt_axi != current_axi
        name(a) = a ? "2D-Axisymmetric" : "2D-Plane"
        push!(errors, "Geometry mismatch: checkpoint was written as $(name(ckpt_axi)), " *
                      "this run is $(name(current_axi)). Set solver_type back to " *
                      "\"$(name(ckpt_axi))\" to resume, or start a new analysis.")
    end


    if isempty(errors)
        return (true, "All dimensions validated successfully")
    else
        error_msg = "Checkpoint validation failed:\n  " * join(errors, "\n  ")
        return (false, error_msg)
    end
end


"""
    load_checkpoint(checkpoint_file::String, mesh, materials)

Load simulation state from checkpoint file and restore global variables.

# Arguments
- `checkpoint_file::String`: Path to the checkpoint file
- `mesh`: Current mesh data structure (for validation)
- `materials`: Current material data structure (for validation)

# Returns
- Named tuple with:
  - `success::Bool`: Whether loading was successful
  - `current_time::Float64`: Simulation time from checkpoint
  - `output_counter::Int`: Output counter from checkpoint
  - `next_output_time::Float64`: Next output time from checkpoint
  - `message::String`: Status or error message
"""
function load_checkpoint(checkpoint_file::String, mesh, materials)
    
    if !isfile(checkpoint_file)
        return (success=false, current_time=0.0, output_counter=0, 
                next_output_time=0.0, message="Checkpoint file not found: $checkpoint_file")
    end
    
    try
        # Load checkpoint data
        checkpoint_data = load(checkpoint_file)
        
        # Validate dimensions
        is_valid, validation_msg = validate_checkpoint_dimensions(checkpoint_data, mesh, materials)
        if !is_valid
            return (success=false, current_time=0.0, output_counter=0, 
                    next_output_time=0.0, message=validation_msg)
        end
        
        # Restore global state variables
        # Primary state variables
        global C_g = checkpoint_data["C_g"]
        global P = checkpoint_data["P"]
        global T = checkpoint_data["T"]
        global v = checkpoint_data["v"]
        
        # Reactive species state
        global C_lime = checkpoint_data["C_lime"]
        global C_caco3 = checkpoint_data["C_caco3"]
        global binder_content = checkpoint_data["binder_content"]
        global degree_of_carbonation = checkpoint_data["degree_of_carbonation"]
        
        # Material properties (per soil type)
        global C_lime_residual = checkpoint_data["C_lime_residual"]
        
        # Dimension parameters (validate but don't overwrite as they're set in initialize_variables)
        global NDim = checkpoint_data["NDim"]
        global Nnodes = checkpoint_data["Nnodes"]
        global Nelements = checkpoint_data["Nelements"]
        global NSoils = checkpoint_data["NSoils"]
        global NGases = checkpoint_data["NGases"]
        
        # Extract time tracking info
        current_time = checkpoint_data["current_time"]
        output_counter = checkpoint_data["output_counter"]
        next_output_time = checkpoint_data["next_output_time"]
        
        return (success=true, current_time=current_time, output_counter=output_counter,
                next_output_time=next_output_time, message="Checkpoint loaded successfully")
        
    catch e
        error_msg = "Error loading checkpoint: $(sprint(showerror, e))"
        return (success=false, current_time=0.0, output_counter=0, 
                next_output_time=0.0, message=error_msg)
    end
end


"""
    get_checkpoint_summary(checkpoint_file::String)

Get summary information from a checkpoint file without fully loading it.

# Arguments
- `checkpoint_file::String`: Path to the checkpoint file

# Returns
- `String`: Formatted summary of checkpoint contents
"""
function get_checkpoint_summary(checkpoint_file::String)
    if !isfile(checkpoint_file)
        return "Checkpoint file not found"
    end
    
    try
        checkpoint_data = load(checkpoint_file)
        
        summary = """
        Checkpoint Summary:
          Nodes: $(checkpoint_data["Nnodes"])
          Elements: $(checkpoint_data["Nelements"])
          Gases: $(checkpoint_data["NGases"])
          Soils: $(checkpoint_data["NSoils"])
          Time: $(checkpoint_data["current_time"])
          Output counter: $(checkpoint_data["output_counter"])
        """
        
        return summary
    catch e
        return "Error reading checkpoint: $(sprint(showerror, e))"
    end
end
