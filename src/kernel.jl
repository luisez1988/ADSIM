#______________________________________________________
# ADSIM: Advection-Diffusion for Soil Improvement and 
# Modification
# v0.1.0
# Author: Luis Zambrano-Cruzatty
#______________________________________________________

#______________________________________________________
# Kernel functions for ADSIM
# Calls procedures in order
#______________________________________________________

# Load required packages
using Dates
using Printf

# Load version information. get_version_string is used by the --version flag in
# main(), so it is imported here rather than relied on from an enclosing module,
# which would leave it undefined when kernel.jl runs directly as a script.
include("version.jl")
using .ADSIMVersion: get_version, get_version_string

# Include data reading modules
include("read_mesh.jl")
include("read_materials.jl")
include("time_functions.jl")
include("read_calc_params.jl")
include("initialize_variables.jl")
include("initialize_flows.jl")
include("time_step.jl")
include("shape_functions.jl")
include("write_vtk.jl")
include("write_probes.jl")
include("fully_explicit_solver.jl")
# IMPES solver. Included after the explicit one because it reuses its free functions
# (mass assembly, element properties, kinetics closures, stab_tau, VTK output); the
# explicit solver itself is unchanged by its presence.
include("impes_solver.jl")
include("write_checkpoint.jl")
include("read_checkpoint.jl")

using .ShapeFunctions
using .WriteVTK
using .WriteProbes


#______________________________________________________
# Main execution script with error handling
#______________________________________________________

"""
    main()

Main execution function with comprehensive error handling.
All errors are logged to the log file before exiting.
"""
function main()
    # Check for --version flag
    if length(ARGS) > 0 && (ARGS[1] == "--version" || ARGS[1] == "-v")
        println(get_version_string())
        exit(0)
    end

    # --report-timestep runs steps [1/8] through [7/8] and stops before the solver,
    # printing the stability spectrum behind the critical time step. It answers "what
    # will my step be, and which term sets it" in seconds on any mesh, without
    # committing to a run that may take hours. Accepted in any argument position.
    report_only = any(a -> a == "--report-timestep", ARGS)
    positional = filter(a -> !startswith(a, "--"), ARGS)

    # Check if project name is provided as command-line argument
    #For debugging only
    #ARGS = ["Advection_test"]
    if isempty(positional)
        println("Error: No project name provided")
        println("Usage: julia kernel.jl <project_name>")
        println("       julia kernel.jl <project_name> --report-timestep")
        println("       julia kernel.jl --version")
        println("Example: julia kernel.jl Test")
        exit(1)
    end

    project_name = positional[1]

    # Construct file paths from project name
    # Assuming data files are in the data/ directory relative to src/
    #data_dir = "src\\data"
    data_dir = "data"
    mesh_file = joinpath(data_dir, "$(project_name).mesh")
    calc_file = joinpath(data_dir, "$(project_name)_calc.toml")
    mat_file = joinpath(data_dir, "$(project_name)_mat.toml")

    # Verify mesh file exists
    if !isfile(mesh_file)
        println("Error: Mesh file not found: ", mesh_file)
        exit(1)
    end

    # Setup output directory
    output_dir = "output"
    if !isdir(output_dir)
        mkdir(output_dir)
    end

    # Determine current stage by checking for existing checkpoint files
    # Look for checkpoint files matching pattern: {project_name}_stage*.jld2
    if isdir(output_dir)
        all_files = readdir(output_dir)
        checkpoint_files = filter(f -> occursin(r"_stage\d+\.jld2$", f) && startswith(f, project_name), all_files)
        
        if isempty(checkpoint_files)
            # No checkpoints found - this is stage 1
            current_stage = 1
        else
            # Extract stage numbers from checkpoint filenames
            stage_numbers = Int[]
            for filename in checkpoint_files
                # Extract number between "stage" and ".jld2"
                m = match(r"_stage(\d+)\.jld2$", filename)
                if m !== nothing
                    push!(stage_numbers, parse(Int, m.captures[1]))
                end
            end
            # Current stage is one more than the highest existing stage
            current_stage = isempty(stage_numbers) ? 1 : maximum(stage_numbers) + 1
        end
    else
        # Output directory doesn't exist yet - this is stage 1
        current_stage = 1
    end

    # Setup log file with stage number. A report-only run writes to its own file: it
    # neither advances nor completes a stage, and clobbering the log of a real run with a
    # diagnostic would lose the record of that run.
    log_file_path = report_only ?
        joinpath(output_dir, "$(project_name)_timestep_report.log") :
        joinpath(output_dir, "$(project_name)_stage$(current_stage).log")

    # Delete existing log file if it exists
    if isfile(log_file_path)
        rm(log_file_path)
    end

    # Open log file for writing
    log_file = open(log_file_path, "w")

    # Custom print function that writes to both console and log file
    function log_print(msg::String)
        println(msg)
        println(log_file, msg)
        flush(log_file)
    end

    try
        # Start timer
        start_time = now()

        # Note: calc and mat files will be checked when needed in future steps
        log_print("Project: $(project_name)")
        log_print("Stage: $(current_stage)")

        log_print("="^64)
        log_print("ADSIM: Advection-Diffusion for Soil Improvement and Modification")
        log_print("Version: $(get_version())")
        log_print("="^64)

        # Step 1: Read mesh data
        log_print("\n[1/8] Reading mesh file: $(mesh_file)")
        mesh = read_mesh_file(mesh_file)
        log_print("   ✓ Loaded $(mesh.num_nodes) nodes and $(mesh.num_elements) elements")
        log_print("   ✓ Loaded initial and boundary conditions")

        # Step 2: Read Material properties
        log_print("\n[2/8] Reading material properties file: $(mat_file)")
        materials = read_materials_file(mat_file)
        log_print("   ✓ Loaded $(length(materials.gas_dictionary)) gases")
        log_print("   ✓ Loaded $(length(materials.soil_dictionary)) soils")

        # Step 3: Read calculation parameters
        log_print("\n[3/8] Reading calculation parameters file: $(calc_file)")
        calc_params = get_all_calc_params(calc_file, data_dir)
        log_print(log_analysis_type(calc_params["solver_settings"]))

        # Time functions drive the boundary conditions in time. A model without
        # them keeps every boundary constant, exactly as before.
        time_functions = calc_params["time_functions"]
        validate_time_function_ids(time_functions, mesh)
        if !isempty(time_functions)
            log_print("   ✓ Loaded $(length(time_functions)) time function(s) for " *
                      "time-dependent boundary conditions")
        end

        log_print("   ✓ Total simulation time: $(calc_params["time_stepping"]["total_simulation_time"]) $(calc_params["units"]["time_unit"])")

        # Axisymmetry is a property of the whole analysis, so it is resolved once here and
        # threaded into the routines whose geometry depends on it.
        axisymmetric = calc_params["solver_settings"]["axisymmetric"]

        # Step 3.4: Gravity must act along the axis of revolution (Eq. ax_gravity). A radial
        # component would vary with the azimuthal angle and so contradict the independence
        # from that angle assumed for every field - a geometric requirement, not a
        # simplification, which is why this is an error rather than a warning.
        if axisymmetric && calc_params["solver_settings"]["gravity"] == 1
            g_r = calc_params["gravity"]["x_component"]
            if abs(g_r) > 0.0
                error("Axisymmetric analysis requires gravity along the axis of revolution, " *
                      "so gravity_x_component must be 0.0, but it is $g_r. A radial gravity " *
                      "component varies with the azimuthal angle and is incompatible with " *
                      "axisymmetry (Appendix: Axisymmetric formulation, Eq. ax_gravity).")
            end
        end

        # Step 3.5: Validate reaction kinetics requirements
        if calc_params["solver_settings"]["reaction_kinetics"] == 1
            log_print("\nValidating reaction kinetics requirements")
            validate_reaction_kinetics_requirements(calc_params["solver_settings"], materials)
            log_print("   ✓ CO2 gas is defined in materials")
        end

        # Step 3.6: Check for existing checkpoint from previous stage
        checkpoint_file, prev_stage = find_latest_checkpoint(project_name, output_dir)
        checkpoint_loaded = false
        initial_state = nothing
        
        if checkpoint_file !== nothing && current_stage > 1
            log_print("\n    Loading checkpoint from previous stage")
            log_print("   Found checkpoint: $(basename(checkpoint_file)) (Stage $(prev_stage))")
            
            # Initialize arrays first (dimensions only)
            zero_variables!(mesh, materials; axisymmetric = axisymmetric)
            
            # Load checkpoint data
            checkpoint_result = load_checkpoint(checkpoint_file, mesh, materials)
            
            if checkpoint_result.success
                checkpoint_loaded = true
                initial_state = (current_time=checkpoint_result.current_time, 
                               output_counter=checkpoint_result.output_counter,
                               next_output_time=checkpoint_result.next_output_time)
                
                checkpoint_size = get_checkpoint_file_size(checkpoint_file)
                log_print("   ✓ Checkpoint loaded successfully ($(checkpoint_size))")
                log_print("   ✓ Restored state at time: $(checkpoint_result.current_time) $(calc_params["units"]["time_unit"])")
                log_print("   ✓ Continuing from output counter: $(checkpoint_result.output_counter)")
            else
                log_print("   ⚠ Warning: $(checkpoint_result.message)")
                log_print("   ⚠ Proceeding with normal initialization instead")
            end
        end

        # Step 4: Initialize simulation variables
        if !checkpoint_loaded
            log_print("\n[4/8] Initializing simulation variables")
            zero_variables!(mesh, materials; axisymmetric = axisymmetric)
            log_print("   ✓ Allocated arrays for $(Nnodes) nodes")
            log_print("   ✓ Tracking $(NGases) gas species in $(NSoils) soil types")
        else
            log_print("\n[4/8] Simulation variables initialized from checkpoint")
            log_print("   ✓ Using $(Nnodes) nodes from checkpoint")
            log_print("   ✓ Tracking $(NGases) gas species in $(NSoils) soil types")
        end

        # Step 5: Apply initial conditions and initialize flows
        #
        # Time-driven boundaries are seeded at the time the run actually starts:
        # 0 for a fresh run, the checkpoint time for a restart, so a later stage
        # picks its curves up where the previous stage left them.
        bc_start_time = checkpoint_loaded ? initial_state.current_time : 0.0

        if !checkpoint_loaded
            log_print("\n[5/8] Applying initial conditions and initializing flows")
            apply_all_initial_conditions!(mesh, materials, time_functions, bc_start_time)
            initialize_all_flows!(mesh, materials, Nnodes, NGases, time_functions, bc_start_time;
                                  axisymmetric = axisymmetric)
            log_print("   ✓ Initial and boundary conditions applied")
        else
            log_print("\n[5/8] Applying boundary conditions from mesh file")
            # Apply boundary conditions from mesh file (may have changed between stages)
            # This sets P_boundary and T_boundary, and applies the prescribed values.
            # apply_temperature_bc! belongs here too: without it T_boundary stays
            # all ones after a restart and prescribed temperature nodes are never
            # frozen.
            apply_concentration_bc!(mesh, time_functions, bc_start_time)
            apply_pressure_bc!(mesh, time_functions, bc_start_time)
            apply_temperature_bc!(mesh, time_functions, bc_start_time)

            # Save evolved state variables that will be overwritten
            # when we calculate Caco3_max
            C_lime_checkpoint = copy(C_lime)
            C_caco3_checkpoint = copy(C_caco3)
            binder_content_checkpoint = copy(binder_content)
            degree_of_carbonation_checkpoint = copy(degree_of_carbonation)
            
            # Recalculate Caco3_max and C_lime_residual from material properties
            # These are derived from initial material properties, not evolved state
            # Must be recalculated because they're not saved in checkpoint
            # Note: This will overwrite C_lime with initial values temporarily
            apply_initial_lime_concentration!(mesh, materials)
            
            # Restore evolved state variables from checkpoint
            global C_lime = C_lime_checkpoint
            global C_caco3 = C_caco3_checkpoint
            global binder_content = binder_content_checkpoint
            global degree_of_carbonation = degree_of_carbonation_checkpoint
            
            # Note: Do NOT reapply initial concentrations or temperature
            # as those come from the checkpoint state
            
            # Initialize flow arrays and boundary influences
            initialize_all_flows!(mesh, materials, Nnodes, NGases, time_functions, bc_start_time;
                                  axisymmetric = axisymmetric)
            log_print("   ✓ Boundary conditions reapplied from mesh file")
            log_print("   ✓ Caco3_max and C_lime_residual recalculated from materials")
            log_print("   ✓ Flow arrays and boundary influences initialized")
        end

        # Step 6: Initialize shape functions and calculate time step information
        log_print("\n[6/8] Initializing shape functions")
        initialize_shape_functions!(mesh; axisymmetric = axisymmetric)
        log_print("   ✓ Shape functions and Jacobians precomputed")
        
        log_print("\n[7/8] Calculating time step information")
        time_data, limiting_scale = calculate_time_step_info(mesh, materials, calc_params)
        log_print(@sprintf("   ✓ Minimum characteristic length: %.3g %s", time_data.h_min, calc_params["units"]["geometry_unit"]))
        log_print(@sprintf("   ✓ Critical time step: %.4g %s", time_data.critical_dt, calc_params["units"]["time_unit"]))
        log_print("   ✓ Limiting time scale: $(limiting_scale)")
        log_print("   ✓ Courant number: $(time_data.courant_number)")
        log_print(@sprintf("   ✓ Actual time step: %.4g %s", time_data.actual_dt, calc_params["units"]["time_unit"]))
        log_print("   ✓ Number of time steps: $(time_data.num_steps)")

        # How that step was arrived at: the measured spectrum of every active operator
        # beside the closed form it replaced.
        if time_data.stability !== nothing
            log_print("")
            log_stability_report(time_data.stability, log_print)
        end

        if report_only
            log_print("\n" * "="^64)
            log_print("Time step report complete - solver not run (--report-timestep)")
            log_print("="^64)
            close(log_file)
            return nothing
        end

        # A table sampled more finely than the time step is aliased: the solver
        # steps straight over samples the file recorded. Say so once, here, now
        # that the step is known.
        warn_if_aliased(time_functions, time_data.actual_dt, log_print)

        # Step 8: Run the solver selected by [solver] time_integration.
        #
        # Everything above this line is solver-agnostic, and so is the checkpoint written
        # below: the IMPES solver's state is the same C_g / T / lime field the explicit one
        # carries (C_t is derived from C_g, not stored), so a multi-stage run may switch
        # integrator between stages.
        final_state = if calc_params["solver_settings"]["impes"]
            impes_solver(mesh, materials, calc_params, time_data, project_name,
                         log_print, initial_state, current_stage)
        else
            fully_explicit_diffusion_solver(mesh, materials, calc_params, time_data,
                                            project_name, log_print, initial_state, current_stage)
        end

        # Write checkpoint file for multi-stage calculations
        log_print("\nWriting checkpoint file for stage $(current_stage)...")
        checkpoint_file = write_checkpoint(project_name, current_stage, 
                                          final_state.current_time, 
                                          final_state.output_counter, 
                                          final_state.next_output_time)
        checkpoint_size = get_checkpoint_file_size(checkpoint_file)
        log_print("   ✓ Checkpoint saved: $(basename(checkpoint_file)) ($(checkpoint_size))")

        # Print total run time
        end_time = now()
        total_time = (end_time - start_time).value / 1000.0  # Convert milliseconds to seconds
        log_print("\nTotal run time: $(total_time) seconds")

        log_print("\n" * "="^64)
        log_print("Calculation completed successfully")
        log_print("="^64)

    catch e
        # Log the error with full details
        log_print("\n" * "="^64)
        log_print("FATAL ERROR OCCURRED")
        log_print("="^64)
        
        # Log error type
        log_print("\nError Type: $(typeof(e))")
        
        # Log error message
        log_print("\nError Message:")
        log_print(sprint(showerror, e))
        
        # Log stack trace
        log_print("\nStack Trace:")
        for (exc, bt) in Base.catch_stack()
            showerror(log_file, exc, bt)
            println(log_file)
            flush(log_file)
        end
        
        log_print("\n" * "="^64)
        log_print("Program terminated due to error")
        log_print("="^64)
        
        # Also print to console
        println("\n" * "="^64)
        println("FATAL ERROR - Check log file for details: $(log_file_path)")
        println("="^64)
        
        # Close log file and exit with error code
        close(log_file)
        exit(1)
        
    finally
        # Ensure log file is closed even if error occurs
        if isopen(log_file)
            close(log_file)
        end
    end
end

function print_banner()
    println("\n" * "="^70)
    println(" "^20 * "ADSIM - Adsorption Simulator")
    println(" "^25 * "Version: $(get_version())")
    println(" "^15 * "Advanced Modeling of Adsorption Processes")
    println("="^70 * "\n")
end

# Execute main function when script is run directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

