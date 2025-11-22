#------------------------------------------------------------------------------
# Example: Using write_vtk.jl module with MeshData from read_mesh.jl
#------------------------------------------------------------------------------

# Include required modules
include("read_mesh.jl")
include("write_vtk.jl")

using .WriteVTK

"""
Example usage of the VTK writer module with ADSIM mesh data.

This example demonstrates how to:
1. Load mesh data using read_mesh_file
2. Prepare simulation results arrays
3. Write VTK output files at different time steps
"""

# Read mesh file
mesh = read_mesh_file("data/Test.mesh")

println("Mesh loaded successfully:")
println("  Number of nodes: ", mesh.num_nodes)
println("  Number of elements: ", mesh.num_elements)

# Initialize simulation result arrays at nodes
Nnodes = mesh.num_nodes
Ngases = 3  # Example: N2, O2, CO2
gas_names = ["N2", "O2", "CO2"]

# Allocate arrays for simulation results
concentrations = zeros(Float64, Nnodes, Ngases)
total_concentration = zeros(Float64, Nnodes)
absolute_pressure = zeros(Float64, Nnodes)
concentration_rates = zeros(Float64, Nnodes, Ngases)
reaction_rates = zeros(Float64, Nnodes)
lime_concentration = zeros(Float64, Nnodes)
co2_concentration = zeros(Float64, Nnodes)
caco3_concentration = zeros(Float64, Nnodes)
degree_of_carbonation = zeros(Float64, Nnodes)
volumetric_binder_content = zeros(Float64, Nnodes)
gas_seepage_velocity = zeros(Float64, Nnodes, 2)  # 2D velocity (x, y)
temperature = zeros(Float64, Nnodes)
temperature_rate = zeros(Float64, Nnodes)

# Example: Initialize with some test values
for i in 1:Nnodes
    x, y = get_node_coordinates(mesh, i)
    
    # Example concentrations
    concentrations[i, 1] = 0.78  # N2
    concentrations[i, 2] = 0.21  # O2
    concentrations[i, 3] = 0.01  # CO2
    
    total_concentration[i] = sum(concentrations[i, :])
    absolute_pressure[i] = 101325.0  # Pa (atmospheric)
    temperature[i] = 293.15 + 10.0 * sin(x) * cos(y)  # K (example spatial variation)
    
    # Example velocity field
    gas_seepage_velocity[i, 1] = 0.001 * cos(x)  # m/s
    gas_seepage_velocity[i, 2] = 0.001 * sin(y)  # m/s
    
    # Carbonation-related fields (example values)
    lime_concentration[i] = 500.0  # kg/m³
    co2_concentration[i] = concentrations[i, 3]
    degree_of_carbonation[i] = 0.1 * (1.0 - exp(-0.1 * i / Nnodes))
    volumetric_binder_content[i] = 0.3
end

# Simulate multiple time steps
time_steps = [0, 1, 2, 3, 4, 5]
times = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]  # Physical time in seconds/hours/days

println("\nWriting VTK files for time series...")

for (step, time) in zip(time_steps, times)
    # Update time-dependent fields (example: concentration evolves with time)
    for i in 1:Nnodes
        # Example: CO2 concentration increases with time
        concentrations[i, 3] = 0.01 * (1.0 + 0.1 * time)
        co2_concentration[i] = concentrations[i, 3]
        
        # Update total concentration
        total_concentration[i] = sum(concentrations[i, :])
        
        # Example: Degree of carbonation increases
        degree_of_carbonation[i] = min(1.0, 0.1 * (1.0 + time) * (1.0 - exp(-0.1 * i / Nnodes)))
        
        # Example: Temperature changes slightly
        x, y = get_node_coordinates(mesh, i)
        temperature[i] = 293.15 + 10.0 * sin(x + 0.1 * time) * cos(y)
        
        if step > 0
            temperature_rate[i] = (temperature[i] - 293.15) / time
            concentration_rates[i, 3] = (concentrations[i, 3] - 0.01) / time
        end
    end
    
    # Update CaCO3 based on carbonation
    global caco3_concentration = lime_concentration .* degree_of_carbonation
    
    # Write VTK file for this time step
    write_vtk_file(
        "output/simulation",
        step,
        time,
        mesh,
        concentrations,
        gas_names,
        total_concentration,
        absolute_pressure,
        concentration_rates,
        reaction_rates,
        lime_concentration,
        co2_concentration,
        caco3_concentration,
        degree_of_carbonation,
        volumetric_binder_content,
        gas_seepage_velocity,
        temperature,
        temperature_rate
    )
end

println("\nSimulation complete!")
println("VTK files written for time steps 0-5")
println("Individual VTK files are available as 'output/simulation_XXXXXX.vtk'")
println("Open any .vtk file in ParaView to visualize the results")
