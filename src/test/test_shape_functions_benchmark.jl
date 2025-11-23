# Benchmark test for shape functions module with threading
# This test compares performance and verifies threading is working

using Printf

# Include required modules
include("../read_mesh.jl")
include("../shape_functions.jl")

using .ShapeFunctions

"""
Benchmark the shape function initialization with timing
"""
function benchmark_shape_functions()
    println("="^70)
    println("Shape Functions Performance Benchmark")
    println("="^70)
    
    # Check threading
    num_threads = Base.Threads.nthreads()
    println("\nJulia Threads Available: $num_threads")
    if num_threads == 1
        println("⚠️  Running single-threaded. Use JULIA_NUM_THREADS or --threads to enable parallelism")
    else
        println("✓ Multi-threading enabled!")
    end
    
    # Read the mesh file
    mesh_file = "../data/TEST.mesh"
    println("\nReading mesh file: $mesh_file")
    mesh = read_mesh_file(mesh_file)
    
    println("Number of nodes: ", mesh.num_nodes)
    println("Number of elements: ", mesh.num_elements)
    
    # Create mesh data structure
    mesh_data = (
        Nelements = mesh.num_elements,
        connectivity = mesh.elements,
        coordinates = mesh.coordinates
    )
    
    println("\n" * "="^70)
    println("Timing Test: Initialize Shape Functions")
    println("="^70)
    
    # Warm-up run (JIT compilation)
    println("\nWarm-up run (compiling)...")
    sf_data = initialize_shape_functions!(mesh_data)
    
    # Timed runs
    num_runs = 10
    times = Float64[]
    
    println("\nPerforming $num_runs timed runs...")
    for i in 1:num_runs
        t_start = time()
        sf_data = initialize_shape_functions!(mesh_data)
        t_end = time()
        elapsed = (t_end - t_start) * 1000  # Convert to milliseconds
        push!(times, elapsed)
        @printf("  Run %2d: %.4f ms\n", i, elapsed)
    end
    
    # Statistics
    mean_time = sum(times) / length(times)
    min_time = minimum(times)
    max_time = maximum(times)
    
    println("\n" * "-"^70)
    @printf("Average time: %.4f ms\n", mean_time)
    @printf("Min time:     %.4f ms\n", min_time)
    @printf("Max time:     %.4f ms\n", max_time)
    println("-"^70)
    
    # Verify correctness
    println("\n" * "="^70)
    println("Verification: Sum of Determinants = Area")
    println("="^70)
    
    total_area = 0.0
    for e in 1:mesh.num_elements
        for p in 1:4
            detJ = get_detJ(sf_data, e, p)
            weight = sf_data.gauss_weights[p]
            total_area += detJ * weight
        end
    end
    
    x_min = minimum(mesh.coordinates[:, 1])
    x_max = maximum(mesh.coordinates[:, 1])
    y_min = minimum(mesh.coordinates[:, 2])
    y_max = maximum(mesh.coordinates[:, 2])
    expected_area = (x_max - x_min) * (y_max - y_min)
    
    @printf("\nCalculated area: %.10f\n", total_area)
    @printf("Expected area:   %.10f\n", expected_area)
    
    relative_error = abs(total_area - expected_area) / expected_area * 100.0
    @printf("Relative error:  %.6f%%\n", relative_error)
    
    if relative_error < 1e-6
        println("\n✓ Correctness verification PASSED!")
    else
        println("\n✗ Correctness verification FAILED!")
    end
    
    # Performance estimate for larger meshes
    println("\n" * "="^70)
    println("Performance Scaling Estimate")
    println("="^70)
    
    time_per_element = mean_time / mesh.num_elements
    @printf("\nTime per element: %.6f ms\n", time_per_element)
    
    # Estimate for different mesh sizes
    mesh_sizes = [100, 500, 1000, 5000, 10000]
    println("\nEstimated initialization times for larger meshes:")
    for size in mesh_sizes
        est_time = time_per_element * size
        if est_time < 1000
            @printf("  %5d elements: %.2f ms\n", size, est_time)
        else
            @printf("  %5d elements: %.2f s\n", size, est_time / 1000)
        end
    end
    
    return mean_time, num_threads
end

# Run the benchmark
mean_time, num_threads = benchmark_shape_functions()

println("\n" * "="^70)
println("Summary")
println("="^70)
@printf("Threads used: %d\n", num_threads)
@printf("Mean initialization time: %.4f ms\n", mean_time)
println("\nTo run with different thread counts:")
println("  Windows: \$env:JULIA_NUM_THREADS=8; julia test_shape_functions_benchmark.jl")
println("  Linux:   JULIA_NUM_THREADS=8 julia test_shape_functions_benchmark.jl")
println("="^70)
