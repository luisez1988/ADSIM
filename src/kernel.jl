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

# Include data reading modules
include("read_mesh.jl")
include("read_materials.jl")

#______________________________________________________
# Main execution script
#______________________________________________________

#starts timer
using Dates
start_time = now()

# Check if project name is provided as command-line argument
#For debugging only
ARGS = ["Test"]
if length(ARGS) < 1
    println("Error: No project name provided")
    println("Usage: julia kernel.jl <project_name>")
    println("Example: julia kernel.jl Test")
    exit(1)
end

project_name = ARGS[1]

# Construct file paths from project name
# Assuming data files are in the data/ directory relative to src/
#data_dir = "src//data"
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

# Step 2: Additional processing steps (to be implemented)
# - Material properties
# - Assembly of system matrices
# - Time stepping
# - Solution
# - Post-processing


#Print total run Time
end_time = now()
total_time = end_time - start_time
log_print("Total run time: $(total_time)")

log_print("\n" * "="^64)
log_print("Mesh reading completed successfully")
log_print("="^64)

# Close log file
close(log_file)

