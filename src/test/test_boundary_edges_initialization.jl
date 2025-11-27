# Test for boundary edge initialization in global variables
# Verifies that boundary_edges are correctly precomputed during initialization

using Printf

# Include required modules in correct order
include("../read_mesh.jl")
include("../read_materials.jl")
include("../initialize_variables.jl")  # Defines NGases global
include("../initialize_flows.jl")

"""
Test that boundary edges are properly initialized as a global variable.
"""
function test_boundary_edges_initialization()
    println("="^70)
    println("Testing Boundary Edge Initialization in Global Variables")
    println("="^70)
    
    # Read mesh and materials
    mesh_file = "../data/Advection_test.mesh"
    mat_file = "../data/Advection_test_mat.toml"
    
    println("\nReading mesh file: $mesh_file")
    mesh = read_mesh_file(mesh_file)
    
    println("Reading materials file: $mat_file")
    materials = read_materials_file(mat_file)
    
    println("\nMesh info:")
    println("  Nodes: ", mesh.num_nodes)
    println("  Elements: ", mesh.num_elements)
    println("  Pressure BC nodes: ", length(mesh.absolute_pressure_bc))
    println("  Gases: ", length(materials.gas_dictionary))
    
    # Initialize flows (which should also initialize boundary_edges)
    println("\nInitializing flows...")
    initialize_all_flows!(mesh, materials, mesh.num_nodes, length(materials.gas_dictionary))
    
    # Access the global boundary_edges variable
    println("\n" * "="^70)
    println("TEST: Verify boundary_edges global variable")
    println("="^70)
    
    println("\nGlobal boundary_edges contains $(length(boundary_edges)) edge(s)")
    
    if length(boundary_edges) > 0
        println("\n✓ PASS: boundary_edges is populated")
        
        println("\nBoundary edge details:")
        for (idx, (elem_id, node_i, node_j, l_e, n_hat)) in enumerate(boundary_edges)
            xi, yi = get_node_coordinates(mesh, node_i)
            xj, yj = get_node_coordinates(mesh, node_j)
            @printf("  Edge %d:\n", idx)
            @printf("    Element: %d\n", elem_id)
            @printf("    Nodes: %d → %d\n", node_i, node_j)
            @printf("    Length: %.6f m\n", l_e)
            @printf("    Normal: [%.6f, %.6f]\n", n_hat[1], n_hat[2])
        end
    else
        println("\n✗ FAIL: boundary_edges is empty")
    end
    
    # Verify other flow arrays are also initialized
    println("\n" * "="^70)
    println("TEST: Verify flow arrays are initialized")
    println("="^70)
    
    tests_passed = 0
    tests_total = 5
    
    if size(q_advection) == (mesh.num_nodes, length(materials.gas_dictionary))
        println("✓ q_advection: ", size(q_advection))
        tests_passed += 1
    else
        println("✗ q_advection has wrong size: ", size(q_advection))
    end
    
    if size(q_diffusion) == (mesh.num_nodes, length(materials.gas_dictionary))
        println("✓ q_diffusion: ", size(q_diffusion))
        tests_passed += 1
    else
        println("✗ q_diffusion has wrong size: ", size(q_diffusion))
    end
    
    if size(q_gravitational) == (mesh.num_nodes, length(materials.gas_dictionary))
        println("✓ q_gravitational: ", size(q_gravitational))
        tests_passed += 1
    else
        println("✗ q_gravitational has wrong size: ", size(q_gravitational))
    end
    
    if size(q_boundary) == (mesh.num_nodes, length(materials.gas_dictionary))
        println("✓ q_boundary: ", size(q_boundary))
        tests_passed += 1
    else
        println("✗ q_boundary has wrong size: ", size(q_boundary))
    end
    
    if length(q_source_sink) == mesh.num_nodes
        println("✓ q_source_sink: length = ", length(q_source_sink))
        tests_passed += 1
    else
        println("✗ q_source_sink has wrong length: ", length(q_source_sink))
    end
    
    println("\n" * "="^70)
    println("SUMMARY: $tests_passed / $tests_total flow array tests passed")
    println("="^70)
    
    if tests_passed == tests_total && length(boundary_edges) > 0
        println("\n✓✓✓ ALL TESTS PASSED ✓✓✓")
    else
        println("\n✗✗✗ SOME TESTS FAILED ✗✗✗")
    end
end

# Run the test
test_boundary_edges_initialization()
