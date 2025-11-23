# Test for shape functions module
# This test verifies that the sum of all Jacobian determinants equals the domain area

using Printf

# Include required modules
include("../read_mesh.jl")
include("../shape_functions.jl")

using .ShapeFunctions

"""
Test that the sum of all Jacobian determinants equals the domain area.

For 2x2 Gaussian quadrature, the area integral is:
Area = Σ_{elements} Σ_{gauss points} det(J) * weight

Since all Gauss weights = 1.0 for 2x2 quadrature, this simplifies to:
Area = Σ_{elements} Σ_{gauss points} det(J)
"""
function test_determinant_sum()
    println("="^70)
    println("Testing Shape Functions - Determinant Sum = Domain Area")
    println("="^70)
    
    # Read the mesh file
    mesh_file = "../data/TEST.mesh"
    println("\nReading mesh file: $mesh_file")
    mesh = read_mesh_file(mesh_file)
    
    println("Number of nodes: ", mesh.num_nodes)
    println("Number of elements: ", mesh.num_elements)
    
    # Find domain boundaries from coordinates
    x_min = minimum(mesh.coordinates[:, 1])
    x_max = maximum(mesh.coordinates[:, 1])
    y_min = minimum(mesh.coordinates[:, 2])
    y_max = maximum(mesh.coordinates[:, 2])
    
    expected_area = (x_max - x_min) * (y_max - y_min)
    
    println("\nDomain boundaries:")
    println("  x: [$x_min, $x_max]")
    println("  y: [$y_min, $y_max]")
    @printf("  Expected area (rectangular): %.10f\n", expected_area)
    
    # Initialize shape functions
    println("\nInitializing shape functions...")
    initialize_shape_functions!(mesh)
    
    println("Shape functions initialized successfully!")
    
    # Sum all determinants
    println("\nCalculating area by summing determinants...")
    total_area = 0.0
    
    for e in 1:mesh.num_elements
        for p in 1:4  # 4 Gauss points
            detJ = get_detJ(e, p)
            total_area += detJ * 1.0  # weight = 1.0 for 2x2 quadrature
        end
    end
    
    @printf("\nCalculated area (from determinants): %.10f\n", total_area)
    @printf("Expected area (rectangular domain):  %.10f\n", expected_area)
    
    # Calculate relative error
    relative_error = abs(total_area - expected_area) / expected_area * 100.0
    @printf("Relative error: %.6f%%\n", relative_error)
    
    # Verify the test passes
    tolerance = 1e-6  # 0.0001% tolerance
    if relative_error < tolerance
        println("\n✓ TEST PASSED: Sum of determinants equals domain area!")
    else
        println("\n✗ TEST FAILED: Significant difference detected!")
    end
    
    # Additional verification: check individual determinants
    println("\n" * "="^70)
    println("Additional Diagnostics")
    println("="^70)
    
    # Check if any determinants are negative (would indicate bad element orientation)
    min_detJ = Inf
    max_detJ = -Inf
    negative_count = 0
    
    for e in 1:mesh.num_elements
        for p in 1:4
            detJ = get_detJ(e, p)
            if detJ < 0
                negative_count += 1
                println("WARNING: Negative determinant in element $e, Gauss point $p: $detJ")
            end
            min_detJ = min(min_detJ, detJ)
            max_detJ = max(max_detJ, detJ)
        end
    end
    
    @printf("\nDeterminant statistics:\n")
    @printf("  Min det(J): %.10f\n", min_detJ)
    @printf("  Max det(J): %.10f\n", max_detJ)
    @printf("  Negative determinants: %d\n", negative_count)
    
    if negative_count > 0
        println("\n⚠ WARNING: Negative determinants detected - check element connectivity!")
    else
        println("\n✓ All determinants positive - element orientations correct!")
    end
    
    # Show sample shape function values
    println("\n" * "="^70)
    println("Sample Shape Function Values at Gauss Point 1")
    println("="^70)
    N = get_N(1)
    B = get_B(1)
    println("N (shape functions): ", N)
    println("Sum of N: ", sum(N), " (should be 1.0)")
    println("\nB (derivatives in isoparametric space):")
    display(B)
    
    return total_area, expected_area, relative_error
end

# Run the test
test_determinant_sum()
