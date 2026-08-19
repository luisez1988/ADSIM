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

using Statistics

"""
    heaviside(x)

Heaviside step function: returns 0 for x < 0, and 1 for x >= 0.

# Arguments
- `x`: Input value

# Returns
- 0 if x < 0, 1 if x >= 0
"""
function heaviside(x)
    return x >= 0 ? 1.0 : 0.0
end

"""
    henry_solubility(T)

Henry solubility of CO2 in water, K_H(T), following the van 't Hoff form of
Eq. (henry_vantHoff) of the manuscript.

    K_H(T) = K_H° exp[2400 (1/T - 1/T°)],   K_H° = 3.3e-4 mol m⁻³ Pa⁻¹, T° = 298.15 K

# Arguments
- `T`: Absolute temperature [K]

# Returns
- Henry solubility [mol m⁻³ Pa⁻¹]
"""
function henry_solubility(T)
    K_H_ref = 3.3e-4    # mol m⁻³ Pa⁻¹ at T_ref (Sander 2015)
    T_ref = 298.15      # K
    return K_H_ref * exp(2400.0 * (1.0 / T - 1.0 / T_ref))
end

"""
    lime_solubility(T, ρ_w)

Solubility limit of portlandite in water, A_aq,sat, from the IUPAC correlation
restated on a volumetric basis in Eq. (solubility) of the manuscript.

    ln(A_aq,sat) = 86.1534 - 3492.14/T - 13.7494 ln(T) + ln(ρ_w)

The correlation is retrograde: solubility falls as temperature rises. It returns
2.0e1 mol m⁻³ at 298.15 K.

Retained as a diagnostic only. The rate law of `extent_of_reaction_rate` no
longer caps the dissolved lime at this limit, so this function is not called
during a run; it remains available for reporting how far above saturation the
uncapped law operates.

# Arguments
- `T`: Absolute temperature [K]
- `ρ_w`: Density of water [kg/m³]

# Returns
- Solubility limit [mol per m³ of water]
"""
function lime_solubility(T, ρ_w)
    return exp(86.1534 - 3492.14 / T - 13.7494 * log(T) + log(ρ_w))
end

"""
    arrhenius_coefficient(k_o, E, T)

Arrhenius rate coefficient k_T of Eq. (arrhenius) of the manuscript.

    k_T = k_o exp(-E / R T)

# Arguments
- `k_o`: Arrhenius factor [m³ mol⁻¹ s⁻¹]
- `E`: Activation energy [J/mol]
- `T`: Absolute temperature [K]

# Returns
- Rate coefficient [m³ mol⁻¹ s⁻¹]
"""
function arrhenius_coefficient(k_o, E, T)
    R_gas = 8.3145  # J mol⁻¹ K⁻¹
    return k_o * exp(-E / (R_gas * T))
end

"""
    interfacial_area_factor(C_g_co2, β_area)

Interfacial-area factor a of the rate law, normalised at atmospheric CO2.

    a = exp[β (C_aq - C_aq,atm)]

a stands for the specific gas-liquid interfacial area available to the reaction.
Forcing CO2 into the pore network increases the area of meniscus in contact with
the gas, so the rate rises faster with pressure than the linear C_aq term alone
accounts for. Calibration against the elemental tests gives β of order 0.5
m³/mol; β = 0 disables the factor and recovers the plain second-order law.

# Reference state
C_aq,atm = K_H(T_ref) x_CO2 P_atm = 1.40e-2 mol per m³ of water is the dissolved
CO2 in equilibrium with the free atmosphere, so a = 1 describes a specimen
sitting in air and a > 1 measures the extra interface that injected CO2 buys.
Moving this reference multiplies a by a constant and divides k_o by the same
constant, so β, the residuals and the fit quality are unchanged by the choice;
only the meaning of k_o changes.

# Temperature
The factor carries NO temperature dependence: both the Henry constant and the
ideal-gas conversion are evaluated at T_ref, so a depends on the CO2 gas
concentration alone. This is deliberate. The interfacial area is set by the
capillary state of the pore network, not by how warm the specimen is, and
letting T into this factor produces a strong spurious feedback in either
direction: through K_H(T) the factor collapses as the specimen self-heats,
while through the ideal-gas term alone it would grow by an order of magnitude
over a typical adiabatic rise. The genuine temperature dependence of the rate
remains, through the Arrhenius coefficient and through C_aq(T) in
`extent_of_reaction_rate`.

# Pairing with k_o
β and k_o are fitted together and are not independent. Because a is large at
injection pressures, k_o is correspondingly small; changing one without
refitting the other rescales the rate by exp(β ΔC_aq).

# Arguments
- `C_g_co2`: CO2 concentration in the gas phase [mol per m³ of gas]
- `β_area`: Interfacial-area coefficient β [m³/mol], 0 disables the factor

# Returns
- `a`: Dimensionless area factor, a > 0
"""
function interfacial_area_factor(C_g_co2, β_area)
    if β_area == 0.0
        return 1.0
    end

    R_gas = 8.3145        # J mol⁻¹ K⁻¹
    T_ref = 298.15        # K, reference temperature of the calibration
    P_atm = 101325.0      # Pa
    x_co2_atm = 420e-6    # mol/mol, CO2 mole fraction of air

    # Dissolved CO2 the local gas concentration corresponds to at T_ref, and the
    # value in equilibrium with the atmosphere. Both per m³ of water.
    C_aq_ref = henry_solubility(T_ref) * R_gas * T_ref * C_g_co2
    C_aq_atm = henry_solubility(T_ref) * x_co2_atm * P_atm

    return exp(β_area * (C_aq_ref - C_aq_atm))
end

"""
    extent_of_reaction_rate(C_g_co2, C_lime, C_r, θ_w, T, k_o, E, β_area)

Non-negative extent-of-reaction rate r of Eq. (reaction_rate) of the manuscript,
per unit volume of pore water.

    r = k_T a C_aq (A_s - A_r) H(A_s - A_r)

with C_aq from Henry's law, k_T from Arrhenius, and a the interfacial-area factor
of [`interfacial_area_factor`](@ref). The law is second order: first order in
dissolved CO2 and first order in the lime remaining above the residual. The
Heaviside holds the rate at zero once A_s reaches A_r.

The area factor supplies pressure sensitivity that the linear C_aq term cannot:
across the calibration pressures C_aq spans only 1.14x while the fitted rate
spans an order of magnitude. Setting β_area = 0 removes it.

# Difference from the solubility-capped form
Earlier versions capped the dissolved lime at the portlandite solubility limit,
r = k_T C_aq min{A_aq,sat, A*_s}. For these mixtures A*_s exceeds A_aq,sat by
two to three orders of magnitude, so the min returned the constant A_aq,sat for
all but the last fraction of a percent of conversion and the law was effectively
zero order in lime: a constant rate followed by an abrupt arrest. Without the cap
the aqueous lime concentration tracks the remaining solid in proportion, the rate
decays as lime is consumed, and at fixed temperature the solution is an
exponential approach to A_r.

The two forms are different chemistry, not two numerical treatments of the same
assumption. The capped form assumes the pore water is saturated with portlandite
and the solid acts as a buffered reservoir; the form used here assumes
dissolution is fast relative to the reaction and never limits it.

`k_o` is not transferable between the two forms. In the capped law it multiplied
A_aq,sat, of order 2e1 mol per m³ of water; here it multiplies (A_s - A_r), of
order 1e3 mol per m³ of total volume. A value calibrated against the capped form
must be refitted. The calibration in `calibration_reactionrate/` fits this form.

# Basis note
`C_lime` and `C_r` are stored per unit total volume and enter on that basis,
while C_aq is per unit volume of water. Their product is the rate per unit volume
of water, so the caller still recovers the total-volume rate as θ_w r and the
downstream conversions are unchanged. `lime_solubility` no longer enters the rate.

# Arguments
- `C_g_co2`: CO2 concentration in the gas phase [mol per m³ of gas]
- `C_lime`: Available lime concentration [mol per m³ of total volume]
- `C_r`: Residual (shielded) lime concentration [mol per m³ of total volume]
- `θ_w`: Volumetric water content [-]
- `T`: Absolute temperature [K]
- `k_o`: Arrhenius factor [m³ mol⁻¹ s⁻¹]
- `E`: Activation energy [J/mol]
- `β_area`: Interfacial-area coefficient β [m³/mol], 0 disables the factor

# Returns
- `r`: Extent-of-reaction rate [mol per m³ of water per s], r >= 0
"""
function extent_of_reaction_rate(C_g_co2, C_lime, C_r, θ_w, T, k_o, E, β_area)
    if θ_w <= 0.0 || C_g_co2 <= 0.0
        return 0.0
    end

    R_gas = 8.3145  # J mol⁻¹ K⁻¹

    # Aqueous CO2 at the gas-liquid interface, Eq. (henry) [mol/m³ of water]
    C_aq = henry_solubility(T) * R_gas * T * C_g_co2

    # Lime remaining above the shielded residual, per unit total volume
    Δ_lime = C_lime - C_r
    A_react = heaviside(Δ_lime) > 0.0 ? Δ_lime : 0.0

    k_T = arrhenius_coefficient(k_o, E, T)
    a = interfacial_area_factor(C_g_co2, β_area)

    return k_T * a * C_aq * A_react
end

#=
IMPLEMENTATION NOTES:
=====================

This module implements a fully explicit finite element solver for transient 
gas diffusion in porous media.

The solver can simulate a combination of diffusion, advection, gravitational flow,
and reaction kinetics (e.g., lime carbonation) using a forward Euler time integration scheme.

Temperature effects from exothermic reactions are also included, but assumed adiabatic, i.e., no heat loss to surroundings including the soil matrix.
=#

using Base.Threads
using Printf
using LinearAlgebra

"""
    ElementProperties

Static, time-invariant properties of one element, resolved once before the time loop.

The solver previously re-derived all of these inside the time loop by calling
`get_element_material` and indexing `materials.soil_dictionary` / `materials.soils`,
once per element per gas per step, and again in each of the velocity and temperature
sweeps. None of it changes in time, so it is resolved once here instead.

# Fields
- `material_idx::Int`: Index of the soil material assigned to the element
- `porosity::Float64`: Porosity n [-]
- `saturation::Float64`: Degree of saturation S_r [-]
- `θ_w::Float64`: Volumetric water content n·S_r [-]
- `θ_g::Float64`: Volumetric gas content n(1-S_r) [-]
- `tortuosity::Float64`: Granular tortuosity τ [-]
- `permeability::Float64`: Intrinsic permeability K [m²]
- `residual_lime::Float64`: Residual lime concentration A_r [mol/m³ total]
- `specific_gravity::Float64`: Specific gravity of solids G_s [-]
- `specific_heat_solids::Float64`: Specific heat of solids c_s [J/(kg·K)]

The Arrhenius pair is deliberately absent: the kinetics belong to the reaction,
not to a material, so k_o and E are read once from the reactant properties rather
than carried per element.
"""
struct ElementProperties
    material_idx::Int
    porosity::Float64
    saturation::Float64
    θ_w::Float64
    θ_g::Float64
    tortuosity::Float64
    permeability::Float64
    residual_lime::Float64
    specific_gravity::Float64
    specific_heat_solids::Float64
    λ_e::Float64
end

"""
    effective_conductivity(soil, materials, n, S_r) -> Float64

Effective thermal conductivity of the mixture, Eq. (fourier) of the manuscript:

    λ_e = (1-n) λ_s + θ_w λ_w + θ_g λ_g

a volume-weighted average of the phase conductivities, assembled the same way as
the volumetric heat capacity C_mix. This is a parallel average and so an upper
bound for a granular assembly, since heat must in reality cross grain contacts;
the manuscript flags λ_e as the natural quantity to calibrate if a simulated
thermal front comes out over-diffused.

The gas conductivity is averaged over the declared species rather than weighted by
concentration: λ_g varies little between the gases involved, and the gas term is a
small fraction of λ_e, so a concentration weighting would add a per-node
dependency for no accuracy.

# Arguments
- `soil`: Soil properties of the material
- `materials`: Material data structure
- `n`: Porosity [-]
- `S_r`: Degree of saturation [-]

# Returns
- Effective thermal conductivity [W/(m·K)]
"""
function effective_conductivity(soil, materials, n, S_r)
    θ_w = n * S_r
    θ_g = n * (1.0 - S_r)

    λ_g = 0.0
    if !isempty(materials.gas_dictionary)
        for gas_name in materials.gas_dictionary
            λ_g += materials.gases[gas_name].thermal_conductivity
        end
        λ_g /= length(materials.gas_dictionary)
    end

    return (1.0 - n) * soil.thermal_conductivity_solids +
           θ_w * materials.liquid.thermal_conductivity +
           θ_g * λ_g
end

"""
    build_element_properties(mesh, materials, C_lime_residual) -> Vector{ElementProperties}

Resolve the static material properties of every element once, before time stepping.

# Arguments
- `mesh`: Mesh data structure
- `materials`: Material data structure
- `C_lime_residual::Vector{Float64}`: Residual lime concentration per soil material

# Returns
- `Vector{ElementProperties}`: One entry per element, indexed by element id
"""
function build_element_properties(mesh, materials, C_lime_residual)
    props = Vector{ElementProperties}(undef, mesh.num_elements)

    for e in 1:mesh.num_elements
        material_idx = get_element_material(mesh, e)
        if material_idx === nothing # No material assigned
            error("Element $e has no material assigned. Check mesh material definitions.")
        end

        soil_name = materials.soil_dictionary[material_idx]
        soil = materials.soils[soil_name]
        n = soil.porosity
        S_r = soil.saturation

        props[e] = ElementProperties(
            material_idx,
            n,
            S_r,
            n * S_r,             # θ_w
            n * (1.0 - S_r),     # θ_g
            soil.granular_tortuosity,
            soil.intrinsic_permeability,
            C_lime_residual[material_idx],
            soil.specific_gravity,
            soil.specific_heat_solids,
            effective_conductivity(soil, materials, n, S_r),
        )
    end

    return props
end

"""
    build_node_owner_elements(mesh) -> Vector{Int}

For each node, the highest-index element containing it.

The reaction and temperature terms are evaluated per node but depend on element-level
properties (θ_w, θ_g, the Arrhenius pair). Written as element loops that *assign* rather
than accumulate, they leave each shared node holding the value written by whichever
element came last, i.e. the highest element index. Resolving that ownership once lets
those terms be evaluated in a node loop, which is both order-independent and trivially
parallel, while reproducing the element-loop result exactly.

Where a node sits on a boundary between two materials the choice is arbitrary, but it is
arbitrary in exactly the same way it already was.

# Arguments
- `mesh`: Mesh data structure

# Returns
- `Vector{Int}`: Owning element id per node, 0 if the node belongs to no element
"""
function build_node_owner_elements(mesh)
    owner = zeros(Int, mesh.num_nodes)

    for e in 1:mesh.num_elements
        for i in 1:4
            node_id = mesh.elements[e, i]
            if e > owner[node_id]
                owner[node_id] = e
            end
        end
    end

    return owner
end

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

    Nelements = mesh.num_elements

    # A node is shared by several elements, so scatter-adding into M straight from a threaded
    # element loop is a read-modify-write race. The element range is split into fixed chunks
    # and each chunk accumulates into its own buffer, which are summed in order after the join.
    #
    # Chunks are indexed explicitly rather than by Threads.threadid(): tasks may migrate
    # between threads, so a threadid-indexed buffer can be shared by two tasks, and
    # threadid() may also exceed Threads.nthreads() when interactive threads are present.
    nchunks = max(1, min(Threads.nthreads(), Nelements))
    M_chunk = [zeros(Float64, length(M)) for _ in 1:nchunks]
    chunk_bounds = round.(Int, range(0, Nelements, length = nchunks + 1))

    # Loop over element chunks
    @threads for c in 1:nchunks
        M_local = M_chunk[c]

        for e in (chunk_bounds[c] + 1):chunk_bounds[c + 1]
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
                M_local[node_id] += M_node
            end
        end
    end

    # Reduce the chunk buffers in a fixed order, so M is reproducible run to run
    for c in 1:nchunks
        M .+= M_chunk[c]
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
    @threads for e in 1:Nelements
        # Element stiffness matrix [4×4]
        K_e = zeros(4, 4)
        
        # Integrate over Gauss points
        for p in 1:4
            # Get Jacobian determinant
            detJ = ShapeFunctions.get_detJ(e, p)

            # Precomputed physical-coordinate derivatives dN/dx = B · J^-1
            dN_dx = ShapeFunctions.get_dN_dx(e, p)  # [4 nodes, 2 coords]

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
function fully_explicit_diffusion_solver(mesh, materials, calc_params, time_data, project_name, log_print, initial_state=nothing, current_stage::Int=1)
    log_print("\n[8/8] Starting fully explicit diffusion solver")
    log_print("   Using $(Threads.nthreads()) threads for parallel execution")

    # Access global variables
    global C_g, P, T, v, P_boundary, λ_bc
    global C_lime, C_caco3, C_lime_residual, binder_content, degree_of_carbonation, Caco3_max
    global dC_g_dt, dT_dt, dC_lime_dt
    global boundary_node_influences
    global q_boundary

    # Track warnings for this step
    negative_conc_warned = Dict{Int, Bool}()  # Track warnings per gas

    # Universal gas constant [J/(mol·K)]
    R = 8.314
    M_caco3= 100.09 #g/mol
    ρ_caco3= 2.71e6 #g/m³
    
    # Reaction enthalpy for lime carbonation [J/mol CO2], read from the material file.
    # Ca(OH)2(s) + CO2(g) -> CaCO3(s) + H2O(l)
    # Negative under the IUPAC convention, i.e. exothermic. Material files that
    # predate the [reactants] table fall back to -113800 J/mol CO2, which is
    # Hess's law on the standard enthalpies of formation at 298.15 K
    # (-986.09, -393.51, -1207.6, -285.83 kJ/mol). That value assumes liquid-phase
    # product water; taking H2O(g) instead would give -69.8 kJ/mol.
    ΔH_r = materials.reactants.reaction_enthalpy  # J/mol CO2

    # Kinetics of the carbonation reaction, k_T = k_o exp(-E/RT). A property of the
    # reaction rather than of a material, so it is read once here instead of being
    # carried per element.
    k_o_reaction = materials.reactants.arrhenius_factor    # m³/(mol·s)
    E_reaction = materials.reactants.activation_energy     # J/mol
    # Interfacial-area coefficient of the same rate law. Fitted jointly with k_o,
    # so the two must be taken from the same calibration.
    β_area_reaction = materials.reactants.interfacial_area_beta  # m³/mol

    # Molar heat capacity of each gas species [J/(mol·K)], in gas_dictionary order.
    # This replaces the single hard-coded mass-basis value the solver used to assume
    # for every gas. The molar form is what Eq. (heat_advection) writes, and it lets
    # the gas contribution to the mixture be assembled straight from the nodal
    # concentrations as (ρc)_g = Σ_i C_g^i c̄_p^i, with no density conversion.
    #
    # Eq. (vol_heat_capacity) writes the gas term of C_mix on a mass basis and
    # Eq. (heat_advection) on a molar basis; they are the same quantity, so the molar
    # form is used for both. Strictly the storage term wants c_v and advection c_p,
    # but the gas carries of order 1e-3 of the mixture capacity, so the distinction is
    # far below anything else in the model.
    c_p_molar = [materials.gases[name].molar_heat_capacity
                 for name in materials.gas_dictionary]

    # Get dimensions
    Nnodes = mesh.num_nodes
    Nelements = mesh.num_elements
    NGases = length(materials.gas_dictionary)
    
    # Get solver settings to determine which fluxes to calculate
    solver_settings = calc_params["solver_settings"]
    calculate_diffusion = solver_settings["diffusion"] == 1
    calculate_advection = solver_settings["advection"] == 1
    calculate_gravity = solver_settings["gravity"] == 1
    calculate_reaction = solver_settings["reaction_kinetics"] == 1
    calculate_heat_conduction = solver_settings["heat_conduction"] == 1
    calculate_heat_advection = solver_settings["heat_advection"] == 1

    # With neither transport mechanism active the temperature field does not evolve at
    # all: the run is isothermal and the reaction enthalpy is ignored, however it is
    # set in the material file. Note this is stronger than the λ_e = 0, v_g = 0 limit
    # of the manuscript (tex:474), which still applies the local heat source; to get
    # that adiabatic limit, enable conduction and set the phase conductivities to zero.
    calculate_heat = calculate_heat_conduction || calculate_heat_advection

    # The rate law takes an Arrhenius pair (k_o, E) in place of the single
    # constant κ used before. Material files written for the old input will
    # parse with k_o = 0 and silently produce no reaction, so flag it here.
    if calculate_reaction
        log_print(@sprintf("   Reaction enthalpy ΔH_r: %.4g J/mol CO2", ΔH_r))
        log_print(@sprintf("   Arrhenius factor k_o: %.4g m³/(mol·s), activation energy E: %.4g J/mol",
                           k_o_reaction, E_reaction))
        if β_area_reaction == 0.0
            log_print("   Interfacial-area coefficient β: 0 (factor disabled, plain second-order law)")
        else
            log_print(@sprintf("   Interfacial-area coefficient β: %.4g m³/mol (a = 1 at atmospheric CO2)",
                               β_area_reaction))
        end
        if k_o_reaction <= 0.0
            log_print("Warning: reaction kinetics is enabled but lime_arrhenius_factor = 0, " *
                      "so no carbonation will occur. Check that the material file defines a " *
                      "[reactants] table with lime_arrhenius_factor and lime_activation_energy.")
        end

        # ΔH_r enters the energy equation as q̇ = ΔH_r dC_CO2/dt with dC_CO2/dt < 0,
        # so the exothermic carbonation reaction requires a negative enthalpy to
        # release heat. A magnitude entered without its sign would cool the soil,
        # which is easy to miss in the output, hence the explicit check.
        if ΔH_r > 0.0
            log_print("Warning: reaction_enthalpy = $(ΔH_r) J/mol is positive, i.e. " *
                      "endothermic, so carbonation will cool the soil. Lime carbonation " *
                      "is exothermic; enter the enthalpy as a negative number.")
        end
    end

    # Get gravity vector from calc_params
    gravity_params = calc_params["gravity"]
    g_magnitude = gravity_params["magnitude"]
    g_x = gravity_params["x_component"]
    g_y = gravity_params["y_component"]
    g_vector = [g_x, g_y] * g_magnitude  # [m/s²]    
    
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

    # Resolve static element material properties and node ownership once, so the time
    # loop never repeats a material dictionary lookup
    elem_props = build_element_properties(mesh, materials, C_lime_residual)
    node_owner = build_node_owner_elements(mesh)

    # CO2 is the only reactive species; resolve its index once
    co2_gas_idx = findfirst(name -> name == "CO2", materials.gas_dictionary)

    #Calculate total gas concentrations
    total_concentration = vec(sum(C_g, dims=2))

    #Calculate the absolute pressure
    P = total_concentration .* R .* T  # Ideal gas law: P = C_total * R * T

    # Initialize time tracking based on checkpoint or from scratch
    if initial_state !== nothing
        # Continue from checkpoint
        current_time = initial_state.current_time
        output_counter = initial_state.output_counter
        next_output_time = initial_state.next_output_time
        log_print("      Continuing from time: $(current_time) $(calc_params["units"]["time_unit"])")
        log_print("      Next output at: $(next_output_time) $(calc_params["units"]["time_unit"])")
    else
        # Start from initial conditions
        # Write initial state (t = 0)
        log_print("      Load step 0 (0.0%)")
        write_output_vtk(mesh, materials, 0, 0.0, project_name, total_concentration)
        
        current_time = 0.0
        next_output_time = load_step_time
        output_counter = 1
    end

    # Open the node probes. These follow individual nodes at data_saving_interval,
    # which is set apart from the VTK cadence, so they can be sampled much more
    # finely than the mesh is dumped. Returns nothing when no nodes were requested.
    probe_writer = WriteProbes.init_node_probes(project_name, current_stage, mesh,
                                                materials.gas_dictionary,
                                                calc_params["probing_nodes"],
                                                Float64(data_saving_interval),
                                                calc_params["units"], current_time,
                                                log_print)

    # Record the starting state, so a series begins at the initial (or restart) time
    WriteProbes.probe_sample!(probe_writer, current_time, C_g, dC_g_dt,
                              total_concentration, P, T, dT_dt, C_lime, C_caco3,
                              dC_lime_dt, degree_of_carbonation, binder_content, v;
                              force=true)

    # Flow vectors, allocated once and zeroed in place each step.
    # These are function-local rather than the like-named globals in initialize_flows.jl:
    # locals are type-stable and avoid a global lookup in the hot loop. The globals of the
    # same name are not used by this solver.
    q_diffusion = zeros(Float64, Nnodes, NGases)
    q_advection = zeros(Float64, Nnodes, NGases)
    q_gravitational = zeros(Float64, Nnodes, NGases)
    total_rate = zeros(Float64, Nnodes)   # Boundary rate of change
    q_source_sink = zeros(Float64, Nnodes) # Only for CO2 for now

    # Element-chunk decomposition for the threaded flux and velocity loops.
    #
    # nchunks is deliberately independent of Threads.nthreads(): the per-chunk buffers are
    # reduced in chunk order, so a fixed chunk count makes the result reproducible to the
    # last bit no matter how many threads run it. Tying it to the thread count would make
    # the summation order, and therefore the round-off, depend on the machine.
    #
    # Memory is nchunks * Nnodes * (3*NGases + NDim) * 8 bytes. For meshes large enough for
    # that to matter, mesh colouring would replace the buffers entirely.
    nchunks = max(1, min(Nelements, 32))
    chunk_bounds = round.(Int, range(0, Nelements, length = nchunks + 1))

    qd_chunk = [zeros(Float64, Nnodes, NGases) for _ in 1:nchunks]
    qa_chunk = [zeros(Float64, Nnodes, NGases) for _ in 1:nchunks]
    qg_chunk = [zeros(Float64, Nnodes, NGases) for _ in 1:nchunks]
    v_chunk = [zeros(Float64, Nnodes, NDim) for _ in 1:nchunks]

    # Per-element scratch, one set per chunk. Sized 4 for the nodes of a quad element.
    q_aux_chunk = [zeros(4) for _ in 1:nchunks]
    ρ_g_chunk = [zeros(4) for _ in 1:nchunks]

    # Gauss-point Darcy velocity, cached per element and Gauss point. Presently written
    # for the nodal projection only; the thermal advective flux of the energy equation
    # consumes it directly (manuscript Eq. thermal_advective_flux), which is why it is
    # retained rather than discarded once projected.
    v_gp_cache = zeros(Float64, Nelements, 4, NDim)

    # Main time stepping loop
    save_data = false
    for step in 1:num_steps

        #reset flow vectors (q_boundary is prefilled and not reset here)
        q_diffusion .= 0.0
        q_advection .= 0.0
        q_gravitational .= 0.0
        total_rate .= 0.0
        q_source_sink .= 0.0
        for c in 1:nchunks
            qd_chunk[c] .= 0.0
            qa_chunk[c] .= 0.0
            qg_chunk[c] .= 0.0
        end

        #______________________________________________________
        #Reaction kinetics calculation start here
        #
        #Evaluated per node rather than per element. The rate depends only on nodal
        #state (C_g, C_lime, T) plus element-level properties, which `node_owner`
        #resolves once, so this reproduces the element loop exactly while being
        #order-independent and parallel. It uses level-k state, so it must run before
        #the gas fluxes consume q_source_sink.
        #______________________________________________________
        if calculate_reaction && co2_gas_idx !== nothing
            @threads for node_id in 1:Nnodes
                e = node_owner[node_id]
                if e == 0
                    continue
                end
                props = elem_props[e]

                #Extent-of-reaction rate r >= 0, per unit volume of water
                r = extent_of_reaction_rate(C_g[node_id, co2_gas_idx], C_lime[node_id],
                                            props.residual_lime, props.θ_w, T[node_id],
                                            k_o_reaction, E_reaction, β_area_reaction)

                #Species rates follow from the stoichiometric coefficients.
                #Lime is stored per unit total volume, so dA/dt = -θ_w r.
                dC_lime_dt[node_id] = -props.θ_w * r

                #CO2 is stored per unit gas volume, so the sink carries θ_w/θ_g.
                if props.θ_g > 0.0
                    q_source_sink[node_id] = -M[node_id] * (props.θ_w / props.θ_g) * r
                else
                    q_source_sink[node_id] = 0.0
                end
            end
        end
        #______________________________________________________

        # Loop over element chunks, with the gas loop moved inside.
        #
        # Threading over elements rather than gas species widens the parallel section from
        # NGases (commonly 2) to the element count. Elements sharing a node would race on
        # the nodal flux arrays, so each chunk accumulates into its own buffer and the
        # buffers are reduced in chunk order afterwards. The chunk count is fixed and
        # independent of Threads.nthreads(), so the reduction order - and hence the result
        # to the last bit - does not depend on how many threads happen to be available.
        @threads for c in 1:nchunks
            qd_c = qd_chunk[c]
            qa_c = qa_chunk[c]
            qg_c = qg_chunk[c]
            q_aux = q_aux_chunk[c]
            ρ_g_buf = ρ_g_chunk[c]

            for e in (chunk_bounds[c] + 1):chunk_bounds[c + 1] #loop elements
                # Get element nodes
                nodes = mesh.elements[e, :]

                # Static element properties, resolved once before the time loop
                props = elem_props[e]

                # Gas volume fraction θ_g = n - θ_w = n(1 - S_r)
                θ_g = props.θ_g
                # Soil tortuosity
                τ = props.tortuosity
                # Intrinsic permeability
                k_intrinsic = props.permeability

            for gas_idx in 1:NGases

                # Get gas diffusion coefficient
                gas_name = materials.gas_dictionary[gas_idx]
                gas = materials.gases[gas_name]
                D_g = gas.diff_coefficient

                #Get gas dynamic viscosity
                μ_g = gas.dynamic_viscosity

                #Get current gas nodal concentrations                
                C_e = [C_g[nodes[i], gas_idx] for i in 1:4]

                #______________________________________________________    
                #Diffusion calculation start here
                #______________________________________________________
                if calculate_diffusion
                    #Update diffusion flow vector ∑_p θ_g^p * D_g^p * k_elm * det(J) * W_p / τ^p                    

                    #local diffusion flow vector (per-chunk buffer)
                    mul!(q_aux, K_elements[e], C_e, θ_g * D_g / τ, 0.0)

                    for i in 1:4 #loop nodes in element       
                        node_id = nodes[i] #global node id
                        qd_c[node_id, gas_idx] +=  q_aux[i]
                    end
                end
                #______________________________________________________              
                #Get total nodal concentrations
                C_t= [total_concentration[nodes[i]] for i in 1:4]

                #Get nodal temperatures
                T_e = [T[nodes[i]] for i in 1:4]

                #______________________________________________________    
                #Advection calculation start here
                #______________________________________________________
                if calculate_advection
                    #Zero nodal advection fluxes (per-chunk buffer)
                    fill!(q_aux, 0.0)

                    # loop Gauss points
                    for p in 1:4
                        # Get shape functions at Gauss point
                        N_p = ShapeFunctions.shape_funcs.N[p]

                        #Evaluate concentration gas species concentration at Gauss point
                        C_gp = 0.0
                        C_gp = N_p' * C_e

                        #Evaluate temperature at Gauss point
                        T_gp = 0.0
                        T_gp = N_p' * T_e

                        # Get Jacobian determinant
                        detJ = ShapeFunctions.get_detJ(e, p)

                        Wp= ShapeFunctions.shape_funcs.gauss_weights[p]

                        # Precomputed physical-coordinate derivatives dN/dx = B · J^-1
                        dN_dx = ShapeFunctions.get_dN_dx(e, p)  # [4 nodes, 2 coords]

                        #Update diffusion flow vector ∑_p K^p * T^p *C^p * k^p_elm *C_tot * det(J) * W_p / μ_g^p
                        #Contract right to left: dN_dx * (dN_dx' * C_t) avoids forming the 4×4 product
                        q_aux .+= (R * k_intrinsic * C_gp * T_gp * detJ * Wp / μ_g) .* (dN_dx * (dN_dx' * C_t))
                    end
                    for i in 1:4 #loop nodes in element       
                        node_id = nodes[i] #global node id
                        qa_c[node_id, gas_idx] +=  q_aux[i]
                    end
                end
                #______________________________________________________

                #______________________________________________________
                #Gravity calculation start here
                #______________________________________________________

                if calculate_gravity
                    #Zero nodal gravity fluxes (per-chunk buffer)
                    fill!(q_aux, 0.0)

                    # loop Gauss points
                    for p in 1:4
                        # Get shape functions at Gauss point
                        N_p = ShapeFunctions.shape_funcs.N[p]

                        #Evaluate concentration gas species concentration at Gauss point
                        C_gp = 0.0
                        C_gp = N_p' * C_e

                        # Get Jacobian determinant
                        detJ = ShapeFunctions.get_detJ(e, p)

                        Wp= ShapeFunctions.shape_funcs.gauss_weights[p]

                        # Precomputed physical-coordinate derivatives dN/dx = B · J^-1
                        dN_dx = ShapeFunctions.get_dN_dx(e, p)  # [4 nodes, 2 coords]

                        #get nodal densities
                        ρ_g = ρ_g_buf
                        fill!(ρ_g, 0.0)
                        for i in 1:4
                            for g in 1:NGases
                                gas_name = materials.gas_dictionary[g]
                                gas = materials.gases[gas_name]
                                ρ_g[i] += C_g[nodes[i], g] * gas.molar_mass
                            end
                        end

                        #Update diffusion flow vector ∑_p R * k_intrinsic * C_gp * T_gp * detJ * Wp / μ_g  * (dN_dx · g) * N_p *∑ M C_g
                        #ρ_g interpolates to a scalar at the Gauss point, so this is a scaled
                        #4-vector rather than the 4×4 outer product the previous form built
                        ρ_g_gp = N_p' * ρ_g
                        q_aux .+= (k_intrinsic * C_gp * detJ * Wp * ρ_g_gp / μ_g) .* (dN_dx * g_vector)
                    end
                    for i in 1:4 #loop nodes in element       
                        node_id = nodes[i] #global node id
                        qg_c[node_id, gas_idx] +=  q_aux[i]
                    end
                end
                #______________________________________________________

            end # all gases done for this element
            end # element chunk done
        end # end element-chunk loop

        # Reduce the per-chunk nodal fluxes in a fixed chunk order
        for c in 1:nchunks
            q_diffusion .+= qd_chunk[c]
            q_advection .+= qa_chunk[c]
            q_gravitational .+= qg_chunk[c]
        end

        # Loop over all gases to form the rate of change
        @threads for gas_idx in 1:NGases
            # Get gas name for this species (needed for the reaction branch)
            gas_name = materials.gas_dictionary[gas_idx]

            # calculate rate of change dC/dt = q_net / M
            for i in 1:Nnodes
                dC_g_dt[i, gas_idx] = ((q_boundary[i, gas_idx] - q_diffusion[i, gas_idx] - q_advection[i, gas_idx] - q_gravitational[i, gas_idx]) * P_boundary[i, gas_idx]) / M[i]

                if gas_name == "CO2" && calculate_reaction
                    #reaction contribution to the CO2 gas-phase rate [mol/(m³ gas·s)]
                    dC_g_dt_rxn = (q_source_sink[i] * P_boundary[i, gas_idx]) / M[i]

                    if C_g[i, gas_idx] + dt * (dC_g_dt[i, gas_idx] + dC_g_dt_rxn) < 0.0
                        #can't consume more CO2 than available: throttle the reaction
                        #so the gas is exactly exhausted at the end of the step, and
                        #scale the lime and heat rates by the same factor.
                        dC_g_dt_rxn_allowed = min(-C_g[i, gas_idx] / dt - dC_g_dt[i, gas_idx], 0.0)
                        scale = dC_g_dt_rxn < 0.0 ?
                                clamp(dC_g_dt_rxn_allowed / dC_g_dt_rxn, 0.0, 1.0) : 0.0
                        dC_g_dt_rxn *= scale
                        dC_lime_dt[i] *= scale
                    end

                    #include reaction source/sink term
                    dC_g_dt[i, gas_idx] += dC_g_dt_rxn
                end
            end
        end # end gas loop

        # Calculate total rate of concentration change (outside threaded loop to avoid race conditions)
        for i in 1:Nnodes
            for gas_idx in 1:NGases
                total_rate[i] += dC_g_dt[i, gas_idx]
            end
        end

        #______________________________________________________
        # Calculate Lagrangian multipliers to enforce total concentration BCs
        #______________________________________________________
        
        # Convert keys to vector for thread-safe iteration
        pressure_bc_nodes = collect(keys(mesh.absolute_pressure_bc))
        
        # Loop over all boundary nodes with prescribed pressure BC
        @threads for j in pressure_bc_nodes
            # Check if node has influence length (should be in boundary_node_influences)
            if haskey(boundary_node_influences, j)
                le = boundary_node_influences[j] #length of influence
                #Get current total concentration at node j
                C_total_j = sum(C_g[j, :])
                C_rate_imposed = 0.0 #placeholder for imposed rate will change if a transient pressure is applied

                #The correction λ/M is applied equally to every gas that is free to move
                #at this node. A gas whose concentration is prescribed here cannot absorb
                #any of it, so the constraint must be shared among the free gases only;
                #dividing by NGases would both under-correct and push the prescribed gas
                #off its boundary value.
                n_free = 0
                for g in 1:NGases
                    n_free += P_boundary[j, g]
                end

                if n_free > 0
                    λ_bc[j] = M[j] * (C_total_j + dt *total_rate[j]  - total_concentration[j]) / (dt * n_free)  # Lagrangian multiplier for node j
                else
                    λ_bc[j] = 0.0  # every gas is prescribed here, nothing left to correct
                end
            else
                λ_bc[j] = 0.0  # No influence length found, set multiplier to zero
            end
        end

        #Loop gases again to apply Lagrangian correction and update concentrations
        @threads for gas_idx in 1:NGases
            # Get gas name for warnings
            gas_name = materials.gas_dictionary[gas_idx]
            
            # Update concentrations: C^(n+1) = C^n + dt * dC/dt
            for i in 1:Nnodes
                # Apply Lagrangian correction only at nodes with boundary influence, and
                # only to gases that are actually free at this node. P_boundary already
                # zeroes dC_g_dt for a prescribed gas; without the same factor here the
                # correction would be added straight through that gate and would walk the
                # prescribed concentration away from its boundary value every step. The
                # drift only shows up once the reaction is active, because λ is built from
                # total_rate, which carries the CO2 sink term.
                lagrangian_correction = 0.0
                if haskey(boundary_node_influences, i) && haskey(mesh.absolute_pressure_bc, i)
                    lagrangian_correction = (λ_bc[i] / M[i]) * P_boundary[i, gas_idx]
                end

                C_g[i, gas_idx] += dt * (dC_g_dt[i, gas_idx] - lagrangian_correction)
                
                # Ensure non-negative and numerically stable concentrations
                C_MIN = 1e-12
                if C_g[i, gas_idx] < C_MIN
                    if C_g[i, gas_idx] < 0.0 && !get(negative_conc_warned, gas_idx, false)
                        log_print("Warning: Negative concentration detected for gas $gas_name at step $step. Setting to zero.")
                        negative_conc_warned[gas_idx] = true
                    end
                    C_g[i, gas_idx] = 0.0
                end
            end

            # Debug: Check for NaN values
            if any(isnan.(C_g[:, gas_idx]))
                nan_nodes = findall(isnan.(C_g[:, gas_idx]))
                log_print("ERROR: NaN detected in gas $gas_name at step $step, nodes: $nan_nodes")
                error("Simulation failed due to NaN values")
            end
        end
        
        # Calculate nodal gas velocities using Darcy's law
        #zero the velocity vector
        v .= 0.0
        for c in 1:nchunks
            v_chunk[c] .= 0.0
        end

        # This stays a separate pass rather than being folded into the flux loop above:
        # it consumes C_g at level k+1 (already updated) together with P at level k, and
        # merging it would silently change which state the velocity is built from.
        #loop over element chunks
        @threads for c in 1:nchunks
            v_c = v_chunk[c]
            ρ_g_vel = ρ_g_chunk[c]

        for e in (chunk_bounds[c] + 1):chunk_bounds[c + 1]
            # Get element nodes
            nodes = mesh.elements[e, :]

            # Static element properties, resolved once before the time loop
            props = elem_props[e]

            # Get intrinsic permeability
            k_intrinsic = props.permeability

            # Get nodal pressures
            P_e = [P[nodes[i]] for i in 1:4]

            # loop Gauss points
            for p in 1:4
                #Get shape functions at Gauss point
                N_p = ShapeFunctions.shape_funcs.N[p]

                # Get Jacobian determinant
                detJ = ShapeFunctions.get_detJ(e, p)

                # Gauss weight
                w = ShapeFunctions.shape_funcs.gauss_weights[p]

                # Precomputed physical-coordinate derivatives dN/dx = B · J^-1
                dN_dx = ShapeFunctions.get_dN_dx(e, p)  # [4 nodes, 2 coords]

                #Evaluate pressure gradient at Gauss point
                grad_P = dN_dx' * P_e  # [2 coords]

                #Evaluate total concentration at Gauss point
                C_total_gp = N_p' * [total_concentration[nodes[i]] for i in 1:4]

                #Calculate concentration-weighted mean dynamic viscosity
                C_TOL = 1e-12  # Numerical tolerance
                μ_g_weighted = 0.0
                if C_total_gp > C_TOL
                    for g in 1:NGases
                        C_g_gp = N_p' * [C_g[nodes[i], g] for i in 1:4]
                        gas_name = materials.gas_dictionary[g]
                        μ_g_weighted += (C_g_gp / C_total_gp) * materials.gases[gas_name].dynamic_viscosity
                    end
                else
                    # Fallback to simple mean if total concentration is negligible
                    μ_g_weighted = mean([materials.gases[materials.gas_dictionary[g]].dynamic_viscosity for g in 1:NGases])
                end

                #Calculate velocity at Gauss point using Darcy's law: v = - (k/μ) ∇P
                v_gp = - (k_intrinsic / μ_g_weighted) * grad_P

                # Calculate mass weight at this Gauss point
                mass_weight = props.θ_g * w * detJ

                #Distribute mass-weighted velocity to nodes
                for i in 1:4
                    node_id = nodes[i]
                    v_c[node_id, :] += v_gp * N_p[i] * mass_weight
                end

                #consider velocity contribution from gravity
                if calculate_gravity
                    #Calculate gravitational velocity at Gauss point (per-chunk buffer)
                    ρ_g = ρ_g_vel
                    fill!(ρ_g, 0.0)
                    for i in 1:4
                        for g in 1:NGases
                            gas_name = materials.gas_dictionary[g]
                            gas = materials.gases[gas_name]
                            ρ_g[i] += C_g[nodes[i], g] * gas.molar_mass
                        end
                    end
                    ρ_g_gp = N_p' * ρ_g

                    v_g_gp = - (k_intrinsic / μ_g_weighted) * ρ_g_gp * g_vector

                    #Distribute mass-weighted gravitational velocity to nodes
                    for i in 1:4
                        node_id = nodes[i]
                        v_c[node_id, :] += v_g_gp * N_p[i] * mass_weight
                    end

                    #Total Gauss-point velocity, retained for the energy equation
                    v_gp = v_gp + v_g_gp
                end

                #Cache the Gauss-point velocity (Eq. Darcy_Gauss). Previously this was
                #projected to nodes and discarded; the thermal advective flux needs it.
                for d in 1:NDim
                    v_gp_cache[e, p, d] = v_gp[d]
                end
            end
        end
        end # end velocity chunk loop

        # Reduce the per-chunk nodal velocities in a fixed chunk order
        for c in 1:nchunks
            v .+= v_chunk[c]
        end

        # Divide accumulated velocities by nodal mass to get average velocity
        for i in 1:Nnodes
            if M[i] > 0.0
                v[i, :] ./= M[i]
            end
        end

        # Update reaction kinetic terms for lime concentration
        if calculate_reaction
            for i in 1:Nnodes
                # Lime consumed this step. The clamp is applied to the INCREMENT, not to
                # C_lime after the fact, so that CaCO3 gains exactly what lime lost.
                # Clamping C_lime alone while crediting CaCO3 with the unclamped rate
                # would create carbonate from lime that was never there, breaking the 1:1
                # stoichiometry of Eq. (carbonation).
                ΔC_lime = dt * dC_lime_dt[i]
                if C_lime[i] + ΔC_lime < 0.0
                    if !negative_lime_warned
                        log_print("Warning: Negative lime concentration detected at step $step. Clamping to zero.")
                        negative_lime_warned = true
                    end
                    ΔC_lime = -C_lime[i]
                end
                C_lime[i] += ΔC_lime

                #Update caco3_concentration
                C_caco3[i] += - ΔC_lime
                #calculate binder content β_b= V_caco3/V_total
                binder_content[i]= C_caco3[i] * M_caco3 / ρ_caco3
                #calculate degree of carbonation DoC= C_caco3/C_caco3_max
                #Caco3_max is zero for a lime-free material, which would otherwise emit
                #NaN into the output field
                degree_of_carbonation[i] = Caco3_max[i] > 0.0 ? C_caco3[i] / Caco3_max[i] : 0.0
            end
        end

        # Calculate temperature change due to reaction (after all gas fluxes are calculated)
        #
        # Evaluated per node using the same element-ownership resolution as the reaction
        # term above: the element loop this replaces also assigned rather than accumulated,
        # so the owning element's properties reproduce it exactly.
        if calculate_heat && calculate_reaction && co2_gas_idx !== nothing
            ρ_w = materials.liquid.density
            c_w = materials.liquid.specific_heat

            @threads for node_id in 1:Nnodes
                e = node_owner[node_id]
                if e == 0
                    continue
                end
                props = elem_props[e]
                ρ_s = props.specific_gravity * ρ_w

                # Volumetric heat capacity of the gas mixture at this node,
                # (ρc)_g = Σ_i C_g^i c̄_p^i per Eq. (heat_advection), already per unit
                # volume of gas so it enters C_mix weighted by θ_g directly.
                ρc_g = 0.0
                for g in 1:NGases
                    ρc_g += C_g[node_id, g] * c_p_molar[g]
                end

                # Calculate mixture volumetric heat capacity
                # C_mix = (1-n)ρ_s*c_s + θ_w*ρ_w*c_w + θ_g*(ρc)_g
                C_mix = (1.0 - props.porosity) * ρ_s * props.specific_heat_solids +
                        props.θ_w * ρ_w * c_w + props.θ_g * ρc_g

                # Heat generation rate from reaction, Eq. (heat_rate): q̇ = ΔH_r * dC_CO2/dt
                # Both factors are negative (exothermic reaction, CO2 consumed),
                # so q̇ > 0 and represents heat released. dC_lime_dt is the molar
                # rate per unit total volume, which matches the basis of C_mix
                # and equals dC_CO2/dt on the same basis by 1:1 stoichiometry.
                if C_mix > 0.0
                    q_dot = ΔH_r * dC_lime_dt[node_id]

                    # Calculate temperature rate of change
                    dT_dt[node_id] = q_dot / C_mix
                else
                    dT_dt[node_id] = 0.0
                end
            end

            # Update temperature: T^(n+1) = T^n + dt * dT/dt
            for i in 1:Nnodes
                T[i] += dt * dT_dt[i]
                
                # Ensure physically reasonable temperatures (above absolute zero)
                if T[i] < 0.0
                    log_print("Warning: Temperature below absolute zero at node $i. Setting to 0.0 K.")
                    T[i] = 0.0
                end
            end
        end

        # Apply partial pressure boundary conditions after temperature update
        # This ensures partial pressure remains constant by updating concentrations
        # based on new temperature: C_g = P_partial / (R * T)
        for (node_id, partial_pressures) in mesh.partial_pressure_bc
            for gas_idx in 1:NGases
                # Recalculate concentration to maintain constant partial pressure
                # P_partial = C_g * R * T  =>  C_g = P_partial / (R * T)
                C_g[node_id, gas_idx] = partial_pressures[gas_idx] / (R * T[node_id])
            end
        end
        
        # Calculate total gas concentrations after all gases are updated
        total_concentration = vec(sum(C_g, dims=2))

        #Apply total pressure boundary condition at nodes
        for node_id in keys(mesh.absolute_pressure_bc)
            # Apply fixed absolute pressure BC by setting total concentration
            total_concentration[node_id] = mesh.absolute_pressure_bc[node_id] / (R * T[node_id])
        end

        #Update pressure using ideal gas law
        P= total_concentration .* R .* T  # Ideal gas law: P = C_total * R * T
        
        # Update current time
        current_time += dt

        # Sample the probed nodes. Placed after the pressure update so the row holds
        # the fully updated state of this step, and outside the VTK block so the
        # probe cadence is independent of the output cadence.
        WriteProbes.probe_sample!(probe_writer, current_time, C_g, dC_g_dt,
                                  total_concentration, P, T, dT_dt, C_lime, C_caco3,
                                  dC_lime_dt, degree_of_carbonation, binder_content, v)

        # Check if we need to save output
        if save_data || step == num_steps

            # Calculate progress percentage
            progress = 100.0 * step / num_steps
            log_print(@sprintf("      Load Step %d (%.1f%%), Time = %.4e %s",
                              output_counter, progress, current_time, 
                              calc_params["units"]["time_unit"]))

            # Write output
            write_output_vtk(mesh, materials, output_counter, current_time, project_name, total_concentration)          

            
            # Update next output time
            next_output_time += load_step_time
            output_counter += 1
            save_data = false
            negative_lime_warned = false
        end

        #update dt to close exactly at next output time
        if current_time+dt > next_output_time
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

    # Flush and close the probe files
    WriteProbes.close_probes!(probe_writer, log_print)

    # Return final time tracking values for checkpoint writing
    return (current_time=current_time, output_counter=output_counter, next_output_time=next_output_time)
end


"""
    write_output_vtk(mesh, materials, step::Int, time::Float64, project_name, total_concentration)

Write VTK output file for the current time step.

# Arguments
- `mesh`: Mesh data structure
- `materials`: Material data structure
- `step::Int`: Output file counter
- `time::Float64`: Current simulation time
- `project_name::String`: Name of the project for output files
- `total_concentration`: Total gas concentration vector
"""
function write_output_vtk(mesh, materials, step::Int, time::Float64, project_name, total_concentration)
    output_dir = "output"
    filename = joinpath(output_dir, project_name)
    
    # Prepare data for VTK output
    gas_names = materials.gas_dictionary    

    
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
        dC_lime_dt,
        C_lime,
        C_caco3,
        degree_of_carbonation,
        binder_content,
        v,
        T,
        dT_dt
    )
end
