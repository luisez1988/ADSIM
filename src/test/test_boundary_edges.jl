# Test for boundary edge identification and normal calculation
# This test verifies the identify_boundary_edges() and calculate_edge_outward_normal() functions

using Printf
using LinearAlgebra

# Include required modules
include("../read_mesh.jl")

"""
Test boundary edge identification and normal calculation.

This test uses the Advection_test.mesh file which has:
- 22 nodes arranged in a 2-column grid (0.0 to 1.0 height, 0.1 width)
- 10 elements (vertical stack)
- Pressure BCs at nodes 1 and 2 (bottom edge)
- Concentration BCs at nodes 21 and 22 (top edge)
"""
function test_boundary_edges()
    println("="^70)
    println("Testing Boundary Edge Identification and Normal Calculation")
    println("="^70)
    
    # Read the mesh file
    mesh_file = "../data/Advection_test.mesh"
    println("\nReading mesh file: $mesh_file")
    mesh = read_mesh_file(mesh_file)
    
    println("Number of nodes: ", mesh.num_nodes)
    println("Number of elements: ", mesh.num_elements)
    
    # Display boundary conditions
    println("\nPressure BCs applied at nodes:")
    for (node_id, pressure) in mesh.absolute_pressure_bc
        x, y = get_node_coordinates(mesh, node_id)
        @printf("  Node %2d: P = %.1f Pa at (%.3f, %.3f)\n", node_id, pressure, x, y)
    end
    
    println("\nConcentration BCs applied at nodes:")
    for (node_id, conc) in mesh.concentration_bc
        x, y = get_node_coordinates(mesh, node_id)
        @printf("  Node %2d: C = %.1f mol/m³ at (%.3f, %.3f)\n", node_id, conc[1], x, y)
    end
    
    # Test 1: Calculate edge normal for the bottom edge
    println("\n" * "="^70)
    println("TEST 1: Manual edge normal calculation (bottom edge: nodes 1→2)")
    println("="^70)
    
    l_e, n_hat = calculate_edge_outward_normal(mesh, 1, 2)
    x1, y1 = get_node_coordinates(mesh, 1)
    x2, y2 = get_node_coordinates(mesh, 2)
    
    @printf("  Node 1: (%.3f, %.3f)\n", x1, y1)
    @printf("  Node 2: (%.3f, %.3f)\n", x2, y2)
    @printf("  Edge vector: [%.3f, %.3f]\n", x2-x1, y2-y1)
    @printf("  Edge length: %.6f m\n", l_e)
    @printf("  Outward normal: [%.6f, %.6f]\n", n_hat[1], n_hat[2])
    @printf("  Normal magnitude: %.6f (should be 1.0)\n", norm(n_hat))
    
    # Expected: horizontal edge pointing right [0.1, 0.0]
    # Outward normal (down) should be [0.0, -1.0]
    expected_normal = [0.0, -1.0]
    normal_error = norm(n_hat - expected_normal)
    @printf("  Expected normal: [%.1f, %.1f]\n", expected_normal[1], expected_normal[2])
    @printf("  Error: %.6e\n", normal_error)
    
    if normal_error < 1e-10
        println("  ✓ PASS: Normal vector is correct")
    else
        println("  ✗ FAIL: Normal vector is incorrect")
    end
    
    # Test 2: Identify all boundary edges
    println("\n" * "="^70)
    println("TEST 2: Automatic boundary edge identification")
    println("="^70)
    
    boundary_edges = identify_boundary_edges(mesh)
    
    println("\nFound $(length(boundary_edges)) boundary edge(s):")
    for (idx, (elem_id, node_i, node_j, l_e, n_hat)) in enumerate(boundary_edges)
        xi, yi = get_node_coordinates(mesh, node_i)
        xj, yj = get_node_coordinates(mesh, node_j)
        @printf("\n  Edge %d:\n", idx)
        @printf("    Element: %d\n", elem_id)
        @printf("    Nodes: %d → %d\n", node_i, node_j)
        @printf("    Node %d: (%.3f, %.3f)\n", node_i, xi, yi)
        @printf("    Node %d: (%.3f, %.3f)\n", node_j, xj, yj)
        @printf("    Length: %.6f m\n", l_e)
        @printf("    Normal: [%.6f, %.6f]\n", n_hat[1], n_hat[2])
        @printf("    |Normal|: %.6f\n", norm(n_hat))
    end
    
    # Expected: 1 boundary edge (nodes 1-2) from element 1
    expected_num_edges = 1
    if length(boundary_edges) == expected_num_edges
        println("\n  ✓ PASS: Found expected number of boundary edges ($expected_num_edges)")
    else
        println("\n  ✗ FAIL: Expected $expected_num_edges edges, found $(length(boundary_edges))")
    end
    
    # Test 3: Verify edge belongs to correct element
    println("\n" * "="^70)
    println("TEST 3: Verify edge-element association")
    println("="^70)
    
    if length(boundary_edges) > 0
        elem_id, node_i, node_j, l_e, n_hat = boundary_edges[1]
        elem_nodes = get_element_nodes(mesh, elem_id)
        
        println("\n  Element $elem_id nodes: ", elem_nodes)
        
        # Check if both boundary nodes are in the element
        nodes_in_element = (node_i in elem_nodes) && (node_j in elem_nodes)
        
        if nodes_in_element
            println("  ✓ PASS: Both boundary nodes ($node_i, $node_j) belong to element $elem_id")
        else
            println("  ✗ FAIL: Boundary nodes do not belong to reported element")
        end
        
        # Check if nodes are consecutive in element connectivity
        found_consecutive = false
        for i in 1:4
            j = (i % 4) + 1
            if (elem_nodes[i] == node_i && elem_nodes[j] == node_j) ||
               (elem_nodes[i] == node_j && elem_nodes[j] == node_i)
                found_consecutive = true
                println("  ✓ PASS: Nodes are consecutive in element connectivity")
                break
            end
        end
        
        if !found_consecutive
            println("  ✗ FAIL: Nodes are not consecutive in element")
        end
    end
    
    # Test 4: Verify normal unit magnitude for all edges
    println("\n" * "="^70)
    println("TEST 4: Verify all normals are unit vectors")
    println("="^70)
    
    all_unit = true
    for (idx, (elem_id, node_i, node_j, l_e, n_hat)) in enumerate(boundary_edges)
        mag = norm(n_hat)
        error = abs(mag - 1.0)
        @printf("  Edge %d: |n_hat| = %.10f, error = %.2e\n", idx, mag, error)
        if error > 1e-10
            all_unit = false
        end
    end
    
    if all_unit
        println("  ✓ PASS: All normals are unit vectors")
    else
        println("  ✗ FAIL: Some normals are not unit vectors")
    end
    
    # Test 5: Create a simple 2-element test case in memory
    println("\n" * "="^70)
    println("TEST 5: Simple square domain test (2 elements)")
    println("="^70)
    
    # Create a simple mesh: 1x1 square divided into 2 triangular quads
    # Nodes at corners: 1:(0,0), 2:(1,0), 3:(0,1), 4:(1,1)
    # Pressure BC on bottom edge (nodes 1,2) and top edge (nodes 3,4)
    
    test_mesh = MeshData()
    test_mesh.num_nodes = 6
    test_mesh.num_elements = 2
    test_mesh.coordinates = [0.0 0.0;   # Node 1
                             1.0 0.0;   # Node 2
                             0.0 1.0;   # Node 3
                             1.0 1.0;   # Node 4
                             0.0 0.5;   # Node 5 (middle left)
                             1.0 0.5]   # Node 6 (middle right)
    
    # Element 1: bottom quad (1,2,6,5)
    # Element 2: top quad (5,6,4,3)
    test_mesh.elements = [1 2 6 5;
                          5 6 4 3]
    
    # Pressure BC on bottom (1,2) and top (3,4) edges
    test_mesh.absolute_pressure_bc[1] = 101325.0
    test_mesh.absolute_pressure_bc[2] = 101325.0
    test_mesh.absolute_pressure_bc[3] = 100000.0
    test_mesh.absolute_pressure_bc[4] = 100000.0
    
    println("\n  Test mesh created:")
    println("    2 elements, 6 nodes")
    println("    Pressure BCs at nodes: 1, 2 (bottom), 3, 4 (top)")
    
    test_edges = identify_boundary_edges(test_mesh)
    
    println("\n  Found $(length(test_edges)) boundary edges:")
    for (idx, (elem_id, node_i, node_j, l_e, n_hat)) in enumerate(test_edges)
        xi, yi = test_mesh.coordinates[node_i, 1], test_mesh.coordinates[node_i, 2]
        xj, yj = test_mesh.coordinates[node_j, 1], test_mesh.coordinates[node_j, 2]
        @printf("    Edge %d: Elem %d, Nodes %d→%d, Length=%.3f, Normal=[%.2f, %.2f]\n",
                idx, elem_id, node_i, node_j, l_e, n_hat[1], n_hat[2])
    end
    
    # Expected: 2 edges (bottom: 1-2, top: 3-4 or 4-3)
    if length(test_edges) == 2
        println("  ✓ PASS: Found 2 boundary edges as expected")
    else
        println("  ✗ FAIL: Expected 2 edges, found $(length(test_edges))")
    end
    
    # Summary
    println("\n" * "="^70)
    println("TEST SUMMARY")
    println("="^70)
    println("All tests completed. Review results above.")
    println("="^70)
end

# Run the test
test_boundary_edges()
