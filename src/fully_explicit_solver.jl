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
        if material_idx === nothing # No material assigned
            error("Element $e has no material assigned. Check mesh material definitions.")
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
assemble_element_stiffness_matrices(mesh)

Assemble and store element geometric stiffness matrices (without material properties).

# Arguments
- `mesh`: Mesh data structure

# Returns
- `Vector{Matrix{Float64}}`: Vector of element geometric stiffness matrices [Nelements][4×4]

# Formula
K_e[i,j] = ∑_p (B · J^-1)^T · (B · J^-1) det(J) W_p
where material properties will be applied later.
"""
function assemble_element_stiffness_matrices(mesh)
    Nelements = mesh.num_elements
    
    # Preallocate vector of element matrices
    K_elements = Vector{Matrix{Float64}}(undef, Nelements)
    
    # Loop over all elements
    for e in 1:Nelements
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
            # dN/dx = B · J^-1
            dN_dx = B * invJ  # [4 nodes, 2 coords]
            
            # Gauss weight
            w = ShapeFunctions.shape_funcs.gauss_weights[p]
            
            # Compute geometric stiffness contribution: K_e += w * detJ * (dN/dx) * (dN/dx)^T
            # This is ∑_p (B · J^-1)^T · (B · J^-1) det(J) W_p
            # Using matrix multiplication: dN_dx * dN_dx' computes all ∇N_i · ∇N_j terms at once
            K_e .+= (w * detJ) .* (dN_dx * dN_dx')
        end
        
        # Store element matrix
        K_elements[e] = K_e
    end
    
    return K_elements
end


"""
    compute_diffusion_flow!(F::Vector{Float64}, K_elements::Vector{Matrix{Float64}}, 
                          C::Vector{Float64}, mesh, materials, gas_idx::Int)

Compute the diffusion flow vector F = -K × C for a specific gas using element matrices.
Material properties are applied during assembly.

# Arguments
- `F::Vector{Float64}`: Flow vector [Nnodes]
- `K_elements::Vector{Matrix{Float64}}`: Element geometric stiffness matrices [Nelements][4×4]
- `C::Vector{Float64}`: Concentration vector [Nnodes]
- `mesh`: Mesh data structure
- `materials`: Material data structure
- `gas_idx::Int`: Index of the gas species
"""
function compute_diffusion_flow!(F::Vector{Float64}, K_elements::Vector{Matrix{Float64}}, 
                                C::Vector{Float64}, mesh, materials, gas_idx::Int)
    # Zero the flow vector
    fill!(F, 0.0)
    
    # Get gas properties
    gas_name = materials.gas_dictionary[gas_idx]
    gas = materials.gases[gas_name]
    D_g = gas.diff_coefficient
    
    # Loop over all elements and assemble contributions
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
        
        # Get element geometric stiffness matrix and apply material properties
        K_e_geom = K_elements[e]
        K_e = D_eff * K_e_geom
        
        # Extract nodal concentrations for this element
        C_e = [C[nodes[i]] for i in 1:4]
        
        # Compute element flow: F_e = -K_e * C_e
        F_e = -K_e * C_e
        
        # Assemble into global flow vector
        for i in 1:4
            node_id = nodes[i]
            F[node_id] += F_e[i]
        end
    end
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
    Nelements = mesh.num_elements
    NGases = length(materials.gas_dictionary)
    
    # Time stepping parameters
    dt = time_data.actual_dt
    total_time = time_data.total_time
    load_step_time = time_data.time_per_step
    num_steps = time_data.num_steps
    data_saving_interval = calc_params["data_saving_interval"]
    
    # Initialize storage arrays
    M = zeros(Float64, Nnodes)  # Lumped mass vector
    
    # Assemble lumped mass vector (same for all gases)
    assemble_lumped_mass_vector!(M, mesh, materials)
    
    # Check for zero or negative masses
    if any(M .<= 0.0)
        error("Lumped mass vector contains zero or negative values!")
    end
    
    # Precompute geometric element stiffness matrices (same for all gases)    
    K_elements = assemble_element_stiffness_matrices(mesh)
    
    # Write initial state (t = 0)
     log_print("   ✓ Load step: 0")
    write_output_vtk(mesh, materials, 0, 0.0, calc_params)
    
    # Initialize time tracking
    current_time = 0.0
    next_output_time = load_step_time
    output_counter = 1
    
    # Main time stepping loop
    save_data = false
    for step in 1:num_steps
        #reset flow vectors
        q_diffusion = zeros(Float64, Nnodes, NGases)
        # Loop over all gases
        @threads for gas_idx in 1:NGases
            # Get gas name for this species (needed for warning messages)
            gas_name = materials.gas_dictionary[gas_idx]
            
            for e in 1:Nelements #loop elements
                # Get element nodes
                nodes = mesh.elements[e, :]
                
                # Get material properties for this element
                material_idx = get_element_material(mesh, e)
                if material_idx === nothing # No material assigned
                    error("Element $e has no material assigned. Check mesh material definitions.")
                end
                
                soil_name = materials.soil_dictionary[material_idx]
                soil = materials.soils[soil_name]
                
                # Calculate gas volume fraction θ_g = n - θ_w = n(1 - S_r)
                θ_g = soil.porosity * (1.0 - soil.saturation)

                # Get soil tortuosity
                τ = soil.granular_tortuosity

                # Get gas diffusion coefficient
                gas_name = materials.gas_dictionary[gas_idx]
                gas = materials.gases[gas_name]
                D_g = gas.diff_coefficient

                #Update diffusion flow vector ∑_p θ_g^p * D_g^p * k_elm * det(J) * W_p / τ^p 

                #nodal concentrations
                C_e = [C_g[nodes[i], gas_idx] for i in 1:4]

                q_aux= zeros(4) #local diffusion flow vector
                q_aux = (θ_g * D_g / τ) * K_elements[e] *  C_e

                for i in 1:4 #loop nodes in element       
                    node_id = nodes[i] #global node id
                    q_diffusion[node_id, gas_idx] +=  q_aux[i]
                end
            end # flux are ready for this gas

            # calculate rate of change dC/dt = q_diffusion / M
            for i in 1:Nnodes
                dC_g_dt[i, gas_idx] = ((q_boundary[i, gas_idx] - q_diffusion[i, gas_idx]) * P_boundary[i]) / M[i]
            end

            # Update concentrations: C^(n+1) = C^n + dt * dC/dt
            for i in 1:Nnodes
                C_g[i, gas_idx] += dt * dC_g_dt[i, gas_idx]
                
                # Ensure non-negative concentrations
                if C_g[i, gas_idx] < 0.0
                    C_g[i, gas_idx] = 0.0
                    # print warning in log_file
                    log_print("Warning: Negative concentration detected at node $i for gas $gas_name. Setting to zero.")
                end
            end            
        end
        
        # Update current time
        current_time += dt

        # Check if we need to save output
        if save_data || step == num_steps
            # Write output
            write_output_vtk(mesh, materials, output_counter, current_time, calc_params)
            
            # Calculate progress percentage
            progress = 100.0 * step / num_steps
            log_print(@sprintf("      Load Step %d/%d (%.1f%%), Time = %.4e %s",
                              step, num_steps, progress, current_time, 
                              calc_params["units"]["time_unit"]))
            
            # Update next output time
            next_output_time += load_step_time
            output_counter += 1
            save_data = false
        end

        #update dt to close exactly at next output time
        if current_time + dt > next_output_time
            dt = next_output_time - current_time
            #activate switch to true
            save_data = true
        else #use original dt
            dt = time_data.actual_dt
            save_data = false
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
   