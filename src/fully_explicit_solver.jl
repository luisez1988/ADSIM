#______________________________________________________
# ADSIM: Advection-Diffusion for Soil Improvement and 
# Modification
# v0.x.x
# Author: Luis Zambrano-Cruzatty
#______________________________________________________

#______________________________________________________
# Fully explicit solver functions for ADSIM
# Solves all governing equations using an explicit time-stepping scheme
#______________________________________________________

#=
IMPLEMENTATION NOTES:
=====================

This module implements a fully explicit finite element solver for transient 
gas diffusion in porous media.

GOVERNING EQUATION:
-------------------
∂(θ_g C_g)/∂t = ∇ · (D_eff ∇C_g)

where:
- θ_g = gas volume fraction = n(1 - S_r) [-]
- C_g = gas concentration [mol/m³]
- D_eff = D_g × τ = effective diffusion coefficient [m²/s]
- D_g = gas diffusion coefficient [m²/s]
- τ = granular tortuosity [-]
- n = porosity [-]
- S_r = degree of saturation [-]

FINITE ELEMENT DISCRETIZATION:
-------------------------------
Semi-discrete form (after spatial discretization):
M dC/dt = F

where:
- M = lumped mass matrix (diagonal, M_i = ∫ θ_g N_i dΩ)
- C = nodal concentration vector
- F = diffusion flow vector = -K × C
- K = stiffness matrix (K_ij = ∫ D_eff ∇N_i · ∇N_j dΩ)

TIME INTEGRATION (Forward Euler):
----------------------------------
C^(n+1) = C^n + Δt × M^(-1) × F^n
C^(n+1) = C^n + Δt × M^(-1) × (-K × C^n)

ALGORITHM:
----------
1. Assemble lumped mass vector M (once, same for all gases)
2. Assemble stiffness matrix K for each gas (once, time-independent)
3. Time stepping loop:
   For each time step:
     For each gas (parallelized):
       a. Compute flow: F = -K × C
       b. Compute rate: dC/dt = F / M (element-wise)
       c. Apply boundary conditions (zero rate at fixed nodes)
       d. Update: C = C + dt × dC/dt
       e. Enforce non-negativity: C = max(C, 0)
4. Write output at specified intervals

STABILITY:
----------
The explicit scheme requires Δt ≤ Δt_crit for stability.
Δt_crit is calculated based on:
- Diffusive time scale: h²τ/(θ_g D_max)
- Advective time scale (for future implementation)
- Reactive time scale (for future implementation)

The actual time step is: Δt = C_N × Δt_crit, where C_N is the Courant number (≤ 1).

PARALLELIZATION:
----------------
- Gas species loop is parallelized using @threads
- Each gas can be solved independently in the same time step
- Element assembly loops are also parallelized where possible

FUTURE EXTENSIONS:
------------------
- Add advection terms (Darcy flow)
- Add chemical reactions (lime carbonation)
- Add heat transfer
- Implement adaptive time stepping
=#

using Base.Threads
using Printf
using LinearAlgebra

"""
    assemble_lumped_mass_vector!(M::Vector{Float64}, mesh, materials)

Assemble the lumped mass vector for all nodes.
Mass lumping sums element contributions to nodes.

# Arguments
- `M::Vector{Float64}`: Lumped mass vector to be filled [Nnodes]
- `mesh`: Mesh data structure
- `materials`: Material data structure

# Formula
For each element, compute M_e = ∫ θ_g N dΩ ≈ θ_g × A_e / 4
where θ_g is the gas volume fraction and A_e is the element area.
"""
function assemble_lumped_mass_vector!(M::Vector{Float64}, mesh, materials)
    fill!(M, 0.0)
    
    # Loop over all elements
    @threads for e in 1:mesh.num_elements
        # Get element nodes
        nodes = mesh.elements[e, :]
        
        # Get material properties for this element
        material_idx = get_element_material(mesh, e)
        if material_idx === nothing
            continue
        end
        
        soil_name = materials.soil_dictionary[material_idx]
        soil = materials.soils[soil_name]
        
        # Calculate gas volume fraction θ_g = n - θ_w = n(1 - S_r)
        θ_g = soil.porosity * (1.0 - soil.saturation)
        
        # Calculate element area using Gaussian quadrature
        A_e = 0.0
        for p in 1:4  # 4 Gauss points
            detJ = ShapeFunctions.get_detJ(e, p)
            w = ShapeFunctions.shape_funcs.gauss_weights[p]
            A_e += w * detJ
        end
        
        # Distribute mass equally to all 4 nodes (lumped mass)
        M_node = θ_g * A_e / 4.0
        
        # Add contribution to each node
        for i in 1:4
            node_id = nodes[i]
            M[node_id] += M_node
        end
    end
end


"""
    assemble_diffusion_matrix!(K::Matrix{Float64}, mesh, materials, gas_idx::Int)

Assemble the global diffusion stiffness matrix for a specific gas.

# Arguments
- `K::Matrix{Float64}`: Global stiffness matrix [Nnodes × Nnodes]
- `mesh`: Mesh data structure
- `materials`: Material data structure
- `gas_idx::Int`: Index of the gas species

# Formula
K_ij = ∫ (D_eff) ∇N_i · ∇N_j dΩ
where D_eff = D_g × τ is the effective diffusion coefficient.
"""
function assemble_diffusion_matrix!(K::Matrix{Float64}, mesh, materials, gas_idx::Int)
    fill!(K, 0.0)
    
    # Get gas properties
    gas_name = materials.gas_dictionary[gas_idx]
    gas = materials.gases[gas_name]
    D_g = gas.diff_coefficient
    
    # Loop over all elements
    for e in 1:mesh.num_elements
        # Get element nodes
        nodes = mesh.elements[e, :]
        
        # Get material properties for this element
        material_idx = get_element_material(mesh, e)
        if material_idx === nothing
            continue
        end
        
        soil_name = materials.soil_dictionary[material_idx]
        soil = materials.soils[soil_name]
        τ = soil.granular_tortuosity
        
        # Effective diffusion coefficient
        D_eff = D_g * τ
        
        # Element stiffness matrix [4×4]
        K_e = zeros(4, 4)
        
        # Integrate over Gauss points
        for p in 1:4
            # Get shape function derivatives in isoparametric coords
            B = ShapeFunctions.get_B(p)
            
            # Get inverse Jacobian and determinant
            invJ = ShapeFunctions.get_invJ(e, p)
            detJ = ShapeFunctions.get_detJ(e, p)
            
            # Transform derivatives to physical coordinates
            # dN/dx = invJ * dN/dξ
            dN_dx = B * invJ  # [4 nodes, 2 coords]
            
            # Gauss weight
            w = ShapeFunctions.shape_funcs.gauss_weights[p]
            
            # Compute stiffness contribution: K_e += w * detJ * D_eff * (dN/dx) * (dN/dx)^T
            for i in 1:4
                for j in 1:4
                    # Compute ∇N_i · ∇N_j = dN_i/dx * dN_j/dx + dN_i/dy * dN_j/dy
                    grad_dot = dN_dx[i, 1] * dN_dx[j, 1] + dN_dx[i, 2] * dN_dx[j, 2]
                    K_e[i, j] += w * detJ * D_eff * grad_dot
                end
            end
        end
        
        # Assemble into global matrix
        for i in 1:4
            I = nodes[i]
            for j in 1:4
                J = nodes[j]
                K[I, J] += K_e[i, j]
            end
        end
    end
end


"""
    compute_diffusion_flow!(F::Vector{Float64}, K::Matrix{Float64}, C::Vector{Float64})

Compute the diffusion flow vector F = -K × C for a specific gas.

# Arguments
- `F::Vector{Float64}`: Flow vector [Nnodes]
- `K::Matrix{Float64}`: Diffusion stiffness matrix [Nnodes × Nnodes]
- `C::Vector{Float64}`: Concentration vector [Nnodes]
"""
function compute_diffusion_flow!(F::Vector{Float64}, K::Matrix{Float64}, C::Vector{Float64})
    # F = -K * C (negative because diffusion goes from high to low concentration)
    mul!(F, K, C)
    F .*= -1.0
end


"""
    apply_boundary_conditions!(dC_dt::Vector{Float64}, C::Vector{Float64}, 
                              mesh, gas_idx::Int)

Apply boundary conditions by zeroing the rate of change at fixed concentration nodes.

# Arguments
- `dC_dt::Vector{Float64}`: Rate of concentration change [Nnodes]
- `C::Vector{Float64}`: Current concentration [Nnodes]
- `mesh`: Mesh data structure
- `gas_idx::Int`: Gas species index
"""
function apply_boundary_conditions!(dC_dt::Vector{Float64}, C::Vector{Float64}, 
                                   mesh, gas_idx::Int)
    # Apply concentration boundary conditions
    for (node_id, concentrations) in mesh.concentration_bc
        # Set rate to zero (concentration is fixed)
        dC_dt[node_id] = 0.0
        # Ensure concentration remains at BC value
        C[node_id] = concentrations[gas_idx]
    end
end


"""
    fully_explicit_diffusion_solver(mesh, materials, calc_params, time_data, log_print)

Main fully explicit solver for gas diffusion in porous media.
Solves the transient diffusion equation using forward Euler time integration.

# Arguments
- `mesh`: Mesh data structure
- `materials`: Material data structure
- `calc_params`: Calculation parameters dictionary
- `time_data`: Time stepping data structure
- `log_print`: Function for logging output

# Governing Equation
∂(θ_g C_g)/∂t = ∇ · (D_eff ∇C_g)

where:
- θ_g = gas volume fraction [-]
- C_g = gas concentration [mol/m³]
- D_eff = D_g × τ = effective diffusion coefficient [m²/s]
- D_g = gas diffusion coefficient [m²/s]
- τ = granular tortuosity [-]

# Time Integration (Forward Euler)
C_g^(n+1) = C_g^n + Δt × (1/M) × F^n

where:
- M = lumped mass vector
- F = diffusion flow vector = -K × C_g
- K = stiffness matrix from diffusion term
"""
function fully_explicit_diffusion_solver(mesh, materials, calc_params, time_data, log_print)
    log_print("\n[8/N] Starting fully explicit diffusion solver")
    
    # Get dimensions
    Nnodes = mesh.num_nodes
    NGases = length(materials.gas_dictionary)
    
    # Time stepping parameters
    dt = time_data.actual_dt
    total_time = time_data.total_time
    num_steps = time_data.num_steps
    data_saving_interval = calc_params["data_saving_interval"]
    
    # Initialize storage arrays
    M = zeros(Float64, Nnodes)  # Lumped mass vector
    K = zeros(Float64, Nnodes, Nnodes)  # Stiffness matrix (reused for each gas)
    F = zeros(Float64, Nnodes)  # Flow vector
    
    # Assemble lumped mass vector (same for all gases)
    log_print("   ✓ Assembling lumped mass vector...")
    assemble_lumped_mass_vector!(M, mesh, materials)
    
    # Check for zero or negative masses
    if any(M .<= 0.0)
        error("Lumped mass vector contains zero or negative values!")
    end
    
    # Precompute stiffness matrices for all gases
    log_print("   ✓ Assembling diffusion stiffness matrices...")
    K_matrices = Vector{Matrix{Float64}}(undef, NGases)
    @threads for gas_idx in 1:NGases
        K_matrices[gas_idx] = zeros(Float64, Nnodes, Nnodes)
        assemble_diffusion_matrix!(K_matrices[gas_idx], mesh, materials, gas_idx)
    end
    
    # Write initial state (t = 0)
    log_print("   ✓ Writing initial state to VTK...")
    write_output_vtk(mesh, materials, 0, 0.0, calc_params)
    
    # Initialize time tracking
    current_time = 0.0
    next_output_time = data_saving_interval
    output_counter = 1
    
    log_print("   ✓ Starting time integration loop...")
    log_print(@sprintf("      Total steps: %d, dt = %.4e %s", num_steps, dt, calc_params["units"]["time_unit"]))
    
    # Main time stepping loop
    for step in 1:num_steps
        # Loop over all gases (can be parallelized)
        @threads for gas_idx in 1:NGases
            # Get current concentration for this gas
            C_current = view(C_g, :, gas_idx)
            
            # Compute diffusion flow: F = -K * C
            K_gas = K_matrices[gas_idx]
            compute_diffusion_flow!(F, K_gas, C_current)
            
            # Compute rate of change: dC/dt = F / M
            for i in 1:Nnodes
                dC_g_dt[i, gas_idx] = F[i] / M[i]
            end
            
            # Apply boundary conditions
            apply_boundary_conditions!(view(dC_g_dt, :, gas_idx), C_current, mesh, gas_idx)
            
            # Update concentrations: C^(n+1) = C^n + dt * dC/dt
            for i in 1:Nnodes
                C_g[i, gas_idx] += dt * dC_g_dt[i, gas_idx]
                
                # Ensure non-negative concentrations
                if C_g[i, gas_idx] < 0.0
                    C_g[i, gas_idx] = 0.0
                end
            end
        end
        
        # Update current time
        current_time += dt
        
        # Check if we need to save output
        if current_time >= next_output_time - 1e-10 || step == num_steps
            # Write output
            write_output_vtk(mesh, materials, output_counter, current_time, calc_params)
            
            # Calculate progress percentage
            progress = 100.0 * step / num_steps
            log_print(@sprintf("      Step %d/%d (%.1f%%), Time = %.4e %s", 
                              step, num_steps, progress, current_time, 
                              calc_params["units"]["time_unit"]))
            
            # Update next output time
            next_output_time += data_saving_interval
            output_counter += 1
        end
    end
    
    log_print("   ✓ Time integration completed")
    log_print(@sprintf("   ✓ Final time: %.4e %s", current_time, calc_params["units"]["time_unit"]))
end


"""
    write_output_vtk(mesh, materials, step::Int, time::Float64, calc_params)

Write VTK output file for the current time step.

# Arguments
- `mesh`: Mesh data structure
- `materials`: Material data structure
- `step::Int`: Output file counter
- `time::Float64`: Current simulation time
- `calc_params`: Calculation parameters
"""
function write_output_vtk(mesh, materials, step::Int, time::Float64, calc_params)
    output_dir = "output"
    filename = joinpath(output_dir, "simulation")
    
    # Prepare data for VTK output
    gas_names = materials.gas_dictionary
    
    # Calculate total concentration
    total_concentration = vec(sum(C_g, dims=2))
    
    # Placeholder arrays for unused fields (filled with zeros for now)
    reaction_rates = zeros(mesh.num_nodes)
    co2_concentration = zeros(mesh.num_nodes)
    caco3_concentration = zeros(mesh.num_nodes)
    degree_of_carbonation = zeros(mesh.num_nodes)
    volumetric_binder_content = zeros(mesh.num_nodes)
    
    # Call VTK writer
    WriteVTK.write_vtk_file(
        filename,
        step,
        time,
        mesh,
        C_g,
        gas_names,
        total_concentration,
        P,
        dC_g_dt,
        reaction_rates,
        C_lime,
        co2_concentration,
        caco3_concentration,
        degree_of_carbonation,
        volumetric_binder_content,
        v,
        T,
        dT_dt
    )
end
   