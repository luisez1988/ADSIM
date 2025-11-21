#------------------------------------------------------------------------------
# Test script to verify lime concentration initialization with multiple materials
#------------------------------------------------------------------------------

# Include required modules
include("read_mesh.jl")
include("read_materials.jl")
include("initialize_variables.jl")

# Load test data
mesh_file = "data/Test.mesh"
mat_file = "data/Test_multi_mat.toml"

println("Reading mesh and materials...")
mesh = read_mesh_file(mesh_file)
materials = read_materials_file(mat_file)

println("Materials available:")
for (idx, soil_name) in enumerate(materials.soil_dictionary)
    soil_props = get_soil_properties(materials, soil_name)
    println("  Material $(idx): $(soil_name) - Lime content: $(soil_props.lime_content)")
end

# Modify material assignments for testing
# Assign different materials to different elements
println("\nModifying material assignments for test...")
for elem_id in 1:mesh.num_elements
    if elem_id <= 38
        mesh.materials[elem_id] = 1  # Soil 1 (lime = 0.1)
    elseif elem_id <= 76
        mesh.materials[elem_id] = 2  # Soil 2 (lime = 0.25)
    else
        mesh.materials[elem_id] = 3  # Soil 3 (lime = 0.05)
    end
end

println("\nInitializing variables...")
zero_variables!(mesh, materials)

println("\nApplying initial conditions...")
apply_all_initial_conditions!(mesh, materials)

println("\n" * "="^64)
println("LIME CONCENTRATION VERIFICATION WITH MULTIPLE MATERIALS")
println("="^64)

# Get expected lime contents
lime_soil1 = get_soil_properties(materials, "Soil 1").lime_content
lime_soil2 = get_soil_properties(materials, "Soil 2").lime_content
lime_soil3 = get_soil_properties(materials, "Soil 3").lime_content

println("\nExpected lime contents:")
println("  Soil 1: $(lime_soil1)")
println("  Soil 2: $(lime_soil2)")
println("  Soil 3: $(lime_soil3)")

# Count nodes with each lime concentration
count_0_1 = sum(C_lime .≈ 0.1)
count_0_25 = sum(C_lime .≈ 0.25)
count_0_05 = sum(C_lime .≈ 0.05)

println("\nNode distribution by lime content:")
println("  Nodes with C_lime ≈ 0.1:  $(count_0_1)")
println("  Nodes with C_lime ≈ 0.25: $(count_0_25)")
println("  Nodes with C_lime ≈ 0.05: $(count_0_05)")
println("  Total nodes: $(length(C_lime))")

# Check for unique values
unique_values = unique(C_lime)
println("\nUnique C_lime values found: $(length(unique_values))")
for val in sort(unique_values)
    println("  $(val): $(sum(C_lime .== val)) nodes")
end

# Verify some nodes belong to elements with different materials
println("\nSample verification:")
println("  Element 10 (Material 1) nodes: ", get_element_nodes(mesh, 10))
println("    Expected lime: $(lime_soil1)")
for node in get_element_nodes(mesh, 10)
    println("    Node $(node): C_lime = $(C_lime[node])")
end

println("\n  Element 50 (Material 2) nodes: ", get_element_nodes(mesh, 50))
println("    Expected lime: $(lime_soil2)")
for node in get_element_nodes(mesh, 50)
    println("    Node $(node): C_lime = $(C_lime[node])")
end

println("\n  Element 90 (Material 3) nodes: ", get_element_nodes(mesh, 90))
println("    Expected lime: $(lime_soil3)")
for node in get_element_nodes(mesh, 90)
    println("    Node $(node): C_lime = $(C_lime[node])")
end

# Check statistics
min_val = minimum(C_lime)
max_val = maximum(C_lime)

println("\nStatistics:")
println("  Minimum C_lime: $(min_val)")
println("  Maximum C_lime: $(max_val)")
println("  Expected range: [$(lime_soil3), $(lime_soil2)]")

success = (min_val ≈ lime_soil3) && (max_val ≈ lime_soil2)
if success
    println("\n✓ SUCCESS: Lime concentrations correctly initialized for multiple materials!")
else
    println("\n✗ ERROR: Unexpected lime concentration values!")
end

println("="^64)
