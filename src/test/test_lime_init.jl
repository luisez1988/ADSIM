#------------------------------------------------------------------------------
# Test script to verify lime concentration initialization
#------------------------------------------------------------------------------

# Include required modules
include("read_mesh.jl")
include("read_materials.jl")
include("initialize_variables.jl")

# Load test data
project_name = "Test"
data_dir = "data"
mesh_file = joinpath(data_dir, "$(project_name).mesh")
mat_file = joinpath(data_dir, "$(project_name)_mat.toml")

println("Reading mesh and materials...")
mesh = read_mesh_file(mesh_file)
materials = read_materials_file(mat_file)

println("\nInitializing variables...")
zero_variables!(mesh, materials)

println("\nApplying initial conditions...")
apply_all_initial_conditions!(mesh, materials)

println("\n" * "="^64)
println("LIME CONCENTRATION VERIFICATION")
println("="^64)

# Get lime content from material 1
soil_name = materials.soil_dictionary[1]
soil_props = get_soil_properties(materials, soil_name)
expected_lime = soil_props.lime_content

println("\nExpected lime content from material '$(soil_name)': $(expected_lime)")
println("\nChecking nodal lime concentrations (C_lime):")
println("  First 10 nodes:")
for i in 1:min(10, length(C_lime))
    println("    Node $(i): C_lime = $(C_lime[i])")
end

println("\n  Last 10 nodes:")
for i in max(1, length(C_lime)-9):length(C_lime)
    println("    Node $(i): C_lime = $(C_lime[i])")
end

# Check if all nodes have the expected value
all_correct = all(C_lime .== expected_lime)
min_val = minimum(C_lime)
max_val = maximum(C_lime)

println("\nStatistics:")
println("  Total nodes: $(length(C_lime))")
println("  Minimum C_lime: $(min_val)")
println("  Maximum C_lime: $(max_val)")
println("  All values correct: $(all_correct)")

if all_correct
    println("\n✓ SUCCESS: All nodes correctly initialized with lime content = $(expected_lime)")
else
    println("\n✗ ERROR: Some nodes have incorrect lime concentration!")
end

println("="^64)
