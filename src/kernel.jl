#______________________________________________________
# ADSIM: Advection-Diffusion for Soil Improvement and 
# Modification
# v0.x.x
# Author: Luis Zambrano-Cruzatty
#______________________________________________________

#______________________________________________________
# Kernel functions for ADSIM
# Calls procedures in order
#______________________________________________________

# Load required packages
using Dates

# Include data reading modules
include("read_mesh.jl")
include("read_materials.jl")
include("read_calc_params.jl")
include("initialize_variables.jl")
include("initialize_flows.jl")


#______________________________________________________
# Main execution script with error handling
#______________________________________________________

"""
    main()

Main execution function with comprehensive error handling.
All errors are logged to the log file before exiting.
"""
function main()
    # Check if project name is provided as command-line argument
    #For debugging only
    #ARGS = ["Test"]
    if length(ARGS) < 1
        println("Error: No project name provided")
        println("Usage: julia kernel.jl <project_name>")
        println("Example: julia kernel.jl Test")
        exit(1)
    end

    project_name = ARGS[1]

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

    # Setup output directory and log file
    output_dir = "output"
    if !isdir(output_dir)
        mkdir(output_dir)
    end

    log_file_path = joinpath(output_dir, "$(project_name).log")

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

        log_print("="^64)
        log_print("ADSIM: Advection-Diffusion for Soil Improvement and Modification")
        log_print("="^64)

        # Step 1: Read mesh data
        log_print("\n[1/N] Reading mesh file: $(mesh_file)")
        mesh = read_mesh_file(mesh_file)
        log_print("   ✓ Loaded $(mesh.num_nodes) nodes and $(mesh.num_elements) elements")
        log_print("   ✓ Loaded initial and boundary conditions")

        # Step 2: Read Material properties
        log_print("\n[2/N] Reading material properties file: $(mat_file)")
        materials = read_materials_file(mat_file)
        log_print("   ✓ Loaded $(length(materials.gas_dictionary)) gases")
        log_print("   ✓ Loaded $(length(materials.soil_dictionary)) soils")

        # Step 3: Read calculation parameters
        log_print("\n[3/N] Reading calculation parameters file: $(calc_file)")
        calc_params = get_all_calc_params(calc_file)
        log_print("   ✓ Solver type: $(calc_params["solver_type"])")
        log_print("   ✓ Total simulation time: $(calc_params["time_stepping"]["total_simulation_time"]) $(calc_params["units"]["time_unit"])")

        # Step 4: Initialize simulation variables
        log_print("\n[4/N] Initializing simulation variables")
        zero_variables!(mesh, materials)
        log_print("   ✓ Allocated arrays for $(Nnodes) nodes")
        log_print("   ✓ Tracking $(NGases) gas species in $(NSoils) soil types")

        # Step 5: Apply initial conditions and initialize flows
        log_print("\n[5/N] Applying initial conditions and initializing flows")
        apply_all_initial_conditions!(mesh, materials)
        initialize_all_flows!(mesh, materials, Nnodes, NGases)
        log_print("   ✓ System ready for simulation")

        # Step 6: Additional processing steps (to be implemented)
        # - Material properties
        # - Assembly of system matrices
        # - Time stepping
        # - Solution
        # - Post-processing

        # Print total run time
        end_time = now()
        total_time = end_time - start_time
        log_print("\nTotal run time: $(total_time)")

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

# Execute main function
main()

