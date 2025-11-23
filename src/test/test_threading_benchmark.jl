# Benchmark test with synthetic larger mesh for better timing resolution

using Printf

# Include required modules
include("../shape_functions.jl")

using .ShapeFunctions

"""
Generate a synthetic mesh for testing
"""
function generate_test_mesh(nx::Int, ny::Int)
    # Generate a structured quad mesh
    Lx = 0.1  # domain size
    Ly = 0.1
    
    dx = Lx / nx
    dy = Ly / ny
    
    # Generate nodes
    nnodes = (nx + 1) * (ny + 1)
    coordinates = zeros(nnodes, 2)
    
    node_id = 1
    for j in 0:ny
        for i in 0:nx
            coordinates[node_id, 1] = i * dx
            coordinates[node_id, 2] = j * dy
            node_id += 1
        end
    end
    
    # Generate elements (quad connectivity)
    nelements = nx * ny
    connectivity = zeros(Int, nelements, 4)
    
    elem_id = 1
    for j in 0:(ny-1)
        for i in 0:(nx-1)
            # Node numbering: bottom-left, bottom-right, top-right, top-left
            n1 = j * (nx + 1) + i + 1
            n2 = n1 + 1
            n3 = n2 + (nx + 1)
            n4 = n1 + (nx + 1)
            
            connectivity[elem_id, 1] = n1
            connectivity[elem_id, 2] = n2
            connectivity[elem_id, 3] = n3
            connectivity[elem_id, 4] = n4
            elem_id += 1
        end
    end
    
    return (Nelements = nelements, connectivity = connectivity, coordinates = coordinates)
end

"""
Benchmark with different mesh sizes and thread counts
"""
function benchmark_threading()
    println("="^70)
    println("Threading Performance Benchmark - Shape Functions")
    println("="^70)
    
    num_threads = Base.Threads.nthreads()
    println("\nThreads available: $num_threads")
    
    # Test different mesh sizes
    mesh_configs = [
        (10, 10, "Small"),    # 100 elements
        (20, 20, "Medium"),   # 400 elements
        (40, 40, "Large"),    # 1600 elements
        (60, 60, "X-Large"),  # 3600 elements
    ]
    
    results = []
    
    for (nx, ny, label) in mesh_configs
        println("\n" * "="^70)
        println("Mesh: $label ($nx × $ny = $(nx*ny) elements)")
        println("="^70)
        
        # Generate mesh - convert to proper structure
        mesh_gen = generate_test_mesh(nx, ny)
        # Create mesh object compatible with read_mesh structure
        mesh = (num_elements = mesh_gen.Nelements, 
                elements = mesh_gen.connectivity, 
                coordinates = mesh_gen.coordinates)
        println("Nodes: ", size(mesh.coordinates, 1))
        println("Elements: ", mesh.num_elements)
        
        # Warm-up
        initialize_shape_functions!(mesh)
        
        # Time multiple runs
        num_runs = 20
        times = Float64[]
        
        for i in 1:num_runs
            t_start = time()
            initialize_shape_functions!(mesh)
            t_end = time()
            elapsed = (t_end - t_start) * 1000  # ms
            push!(times, elapsed)
        end
        
        mean_time = sum(times) / length(times)
        min_time = minimum(times)
        std_time = sqrt(sum((times .- mean_time).^2) / length(times))
        
        @printf("\nTiming results (%d runs):\n", num_runs)
        @printf("  Mean:   %.4f ms\n", mean_time)
        @printf("  Min:    %.4f ms\n", min_time)
        @printf("  Std:    %.4f ms\n", std_time)
        
        # Verify correctness
        total_area = 0.0
        for e in 1:mesh.num_elements
            for p in 1:4
                detJ = get_detJ(e, p)
                total_area += detJ
            end
        end
        expected_area = 0.1 * 0.1  # 0.01
        error = abs(total_area - expected_area) / expected_area * 100
        
        @printf("\nArea verification:\n")
        @printf("  Calculated: %.10f\n", total_area)
        @printf("  Expected:   %.10f\n", expected_area)
        @printf("  Error:      %.6f%%\n", error)
        
        if error < 1e-6
            println("  ✓ PASSED")
        else
            println("  ✗ FAILED")
        end
        
        push!(results, (nx*ny, mean_time, label))
    end
    
    # Summary table
    println("\n" * "="^70)
    println("Performance Summary")
    println("="^70)
    @printf("Threads: %d\n\n", num_threads)
    @printf("%-15s %10s %15s %15s\n", "Mesh", "Elements", "Time (ms)", "Time/Element (μs)")
    println("-"^70)
    
    for (nelems, time_ms, label) in results
        time_per_elem = (time_ms * 1000) / nelems  # Convert to microseconds
        @printf("%-15s %10d %15.4f %15.4f\n", label, nelems, time_ms, time_per_elem)
    end
    
    return results
end

# Run benchmark
results = benchmark_threading()

println("\n" * "="^70)
println("Note: For comparison, run with different thread counts:")
println("  1 thread:  \$env:JULIA_NUM_THREADS=1; julia test_threading_benchmark.jl")
println("  4 threads: \$env:JULIA_NUM_THREADS=4; julia test_threading_benchmark.jl")
println("  8 threads: \$env:JULIA_NUM_THREADS=8; julia test_threading_benchmark.jl")
println("="^70)
