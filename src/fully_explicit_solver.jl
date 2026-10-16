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

Return 0 for x < 0 and 1 for x >= 0.
"""
function heaviside(x)
    return x >= 0 ? 1.0 : 0.0
end

"""
    stab_tau(v, h_e, D_eff; V_TOL)

Compute the streamline-upwind stabilization parameter τ* [s], using the Brooks & Hughes (1982)
approximation of the optimal upwind value (Christie, Griffiths, Mitchell & Zienkiewicz 1976).

    τ* ≈ (h_e / 2|v|) · min(Pe/3, 1),   Pe = |v| h_e / D_eff

Limits: τ* → h_e/(2|v|) for Pe → ∞ and τ* → 0 for Pe → 0, same as the coth-based
formula, without evaluating coth at every Gauss point. Derivation in
`src/ADVECTION_STABILIZATION_NOTES.md`.

Guards (both cases occur in the verification suite):
- `|v| ≈ 0`: returns 0.
- `D_eff <= 0` (e.g. `axisym_darcy_radial`, `diffusion = 0`): returns the Pe → ∞ limit
  `h_e/(2|v|)` without forming Pe.

# Arguments
- `v`: Gauss-point velocity of the flux term [m/s]
- `h_e`: Element characteristic length [m]
- `D_eff`: Effective diffusivity of the same operator being stabilized [m²/s]
- `V_TOL`: Velocity magnitude treated as zero [m/s]
"""
function stab_tau(v::AbstractVector, h_e::Float64, D_eff::Float64; V_TOL::Float64 = 1.0e-15)
    vnorm = sqrt(v[1]^2 + v[2]^2)
    vnorm < V_TOL && return 0.0
    D_eff <= 0.0 && return h_e / (2.0 * vnorm)
    return (h_e / (2.0 * vnorm)) * min(h_e * vnorm / (6.0 * D_eff), 1.0)
end

"""
    henry_solubility(T)

Compute the Henry solubility of CO2 in water, K_H(T), Eq. (henry_vantHoff).

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

Compute the portlandite solubility limit A_aq,sat (IUPAC correlation, volumetric basis,
Eq. (solubility)).

    ln(A_aq,sat) = 86.1534 - 3492.14/T - 13.7494 ln(T) + ln(ρ_w)

Retrograde: solubility decreases with temperature. 2.0e1 mol m⁻³ at 298.15 K.

Diagnostic only. Not called during a run, since `extent_of_reaction_rate` does not cap
the dissolved lime.

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

Compute the Arrhenius rate coefficient k_T, Eq. (arrhenius).

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
    interfacial_area_factor(C_g_total, T, β_area)

Compute the interfacial-area factor a of the rate law, normalised at atmospheric pressure.

    a = exp[β (p / P_atm - 1)]

a represents the specific gas-liquid interfacial area available to the reaction. Pore
pressurisation increases the meniscus area in contact with the gas, therefore the rate
increases faster with pressure than the linear C_aq term predicts. β = 0 disables the
factor (plain second-order law).

# Driver: TOTAL gas pressure
    p = C_g,total R T,    C_g,total = Σ_α C_g,α

Every species in the mixture contributes to the pressure on the menisci. The calibration
tests also report a total chamber pressure (101.3 kPa + gauge). Using p_CO2 instead
gives a ≈ 1 wherever the pore gas is mostly air, which covers most of an injection run.

# Exponent form
Relative overpressure, dimensionless, therefore β is dimensionless. An earlier revision
used a = exp[β K_H(T_ref) (p_CO2 - x_CO2 P_atm)], with β in m³/mol. Both forms are
exponential in a pressure and differ in scale by K_H(T_ref) P_atm = 33.44 mol/m³
(β = 0.47 in the old form is 15.7 in this one). Old β values must be refitted or rescaled.

The old exponent was an absolute dissolved concentration. With a total gas concentration
as input it grew without bound (ordinary air gave a ~ 9e6). A pressure ratio is bounded by
the applied overpressure.

# Reference state
a = 1 at p = P_atm, giving the plain second-order rate at ambient pressure. Moving the
reference scales a by a constant and k_o by its inverse. β, residuals and fit quality are
unchanged; only the meaning of k_o changes.

# Below atmospheric
Exponent floored at zero, so a >= 1. The calibration covers 5 to 20 kPa GAUGE pressure,
therefore values below ambient would be extrapolation. The floor keeps the plain
second-order rate as the lower bound.

# Pairing with k_o
β and k_o are fitted together. Large a at injection pressures implies small k_o. Changing
one without refitting the other rescales the rate by exp(β Δp/P_atm).

# Arguments
- `C_g_total`: TOTAL gas concentration summed over species [mol per m³ of gas], not the CO2 component
- `T`: Absolute temperature [K], for p = C_g,total R T
- `β_area`: Interfacial-area coefficient β [-]; 0 disables the factor
"""
function interfacial_area_factor(C_g_total, T, β_area)
    if β_area == 0.0
        return 1.0
    end

    R_gas = 8.3145        # J mol⁻¹ K⁻¹
    P_atm = 101325.0      # Pa

    # Total gas pressure, ideal gas law over the whole mixture
    p = C_g_total * R_gas * T

    # Relative overpressure floored at zero (a >= 1). Calibration covers overpressure only
    overpressure = max(p / P_atm - 1.0, 0.0)

    return exp(β_area * overpressure)
end

"""
    extent_of_reaction_rate(C_g_co2, C_g_total, C_lime, C_r, θ_w, T, k_o, E, β_area)

Compute the non-negative extent-of-reaction rate r per unit volume of pore water,
Eq. (reaction_rate).

    r = k_T a C_aq (A_s - A_r) H(A_s - A_r)

C_aq from Henry's law, k_T from Arrhenius, a from [`interfacial_area_factor`](@ref).
Second order: first order in dissolved CO2 and in lime above the residual. The Heaviside
sets r = 0 once A_s reaches A_r.

The area factor provides the pressure sensitivity: across the calibration pressures C_aq
varies only 1.14x, while the fitted rate varies an order of magnitude. β_area = 0 removes it.

# Gas concentration per term
  - `C_aq` (dissolved CO2): CO2 partial concentration. Henry's law dissolves each species
    according to its own partial pressure.
  - `a` (gas-liquid interface): TOTAL gas pressure, consistent with the calibration chamber
    pressure.

C_aq on the CO2 component keeps the CO2 sink first order in CO2. Using the total instead
decouples the sink from the local CO2 amount, and the IMPES positivity limiter then drives
the step to zero at the advancing front (low CO2, high total).

# Difference from the solubility-capped form
Earlier versions used r = k_T C_aq min{A_aq,sat, A*_s}. For these mixtures A*_s exceeds
A_aq,sat by two to three orders of magnitude, therefore the min returned A_aq,sat for
nearly all of the conversion and the law was effectively zero order in lime (constant rate,
then abrupt arrest). Without the cap the rate decays with the remaining lime; at fixed
temperature the solution is an exponential approach to A_r.

Different chemistry: the capped form assumes portlandite-saturated pore water with the
solid as a buffered reservoir; this form assumes dissolution is fast relative to reaction.

`k_o` is not transferable between forms. Capped: k_o multiplied A_aq,sat (~2e1 mol per m³
of water). Here: k_o multiplies (A_s - A_r) (~1e3 mol per m³ of total volume). The
calibration in `calibration_reactionrate/` fits this form.

# Basis
`C_lime` and `C_r` per unit total volume, C_aq per unit volume of water. The product is a
rate per unit volume of water; the caller obtains the total-volume rate as θ_w r.

# Arguments
- `C_g_co2`: CO2 concentration in the gas phase [mol per m³ of gas]; drives C_aq
- `C_g_total`: Total gas concentration, summed over species [mol per m³ of gas];
  drives the interfacial-area factor through the gas pressure
- `C_lime`: Available lime concentration [mol per m³ of total volume]
- `C_r`: Residual (shielded) lime concentration [mol per m³ of total volume]
- `θ_w`: Volumetric water content [-]
- `T`: Absolute temperature [K]
- `k_o`: Arrhenius factor [m³ mol⁻¹ s⁻¹]
- `E`: Activation energy [J/mol]
- `β_area`: Interfacial-area coefficient β [-]; 0 disables the factor

# Returns
- `r` [mol per m³ of water per s], r >= 0
"""
function extent_of_reaction_rate(C_g_co2, C_g_total, C_lime, C_r, θ_w, T, k_o, E, β_area)
    if θ_w <= 0.0 || C_g_co2 <= 0.0
        return 0.0
    end

    R_gas = 8.3145  # J mol⁻¹ K⁻¹

    # Aqueous CO2 at the gas-liquid interface, Eq. (henry) [mol/m³ of water].
    # CO2 partial concentration (Henry's law per species)
    C_aq = henry_solubility(T) * R_gas * T * C_g_co2

    # Lime remaining above the shielded residual, per unit total volume
    Δ_lime = C_lime - C_r
    A_react = heaviside(Δ_lime) > 0.0 ? Δ_lime : 0.0

    k_T = arrhenius_coefficient(k_o, E, T)

    # Area factor driven by TOTAL gas concentration
    a = interfacial_area_factor(C_g_total, T, β_area)

    return k_T * a * C_aq * A_react
end


"""
    reaction_decay_constant(C_g_total, C_lime, C_r, θ_w, θ_g, T, k_o, E, β_area)

Compute the reaction decay constant λ: CO2 consumption rate per unit gas-phase CO2.

With lime and area factor fixed, the reactive term of the CO2 equation is first order:

    dC_g,CO2/dt|rxn = -(θ_w/θ_g) k_T a K_H R T (A_s - A_r) C_g,CO2 = -λ C_g,CO2

λ is the reaction eigenvalue of the nodal equation and enters the explicit stability limit.
`1/λ` is the time to deplete a node of CO2 without resupply.

λ decreases as lime is consumed and as temperature increases: Arrhenius exponent -E/RT,
Henry exponent +2400/T. For E = 7.5 kJ/mol the product scales as T exp(1498/T), and the
retrograde solubility dominates. The reaction is therefore fastest at the coldest region,
typically the face of a specimen at a convective boundary.

# Arguments
- `C_g_total`: total gas concentration [mol/m³ of gas], for the area factor
- `C_lime`: lime remaining, A_s [mol/m³ total]
- `C_r`: shielded residual lime, A_r [mol/m³ total]
- `θ_w`, `θ_g`: volumetric water and gas contents [-]
- `T`: absolute temperature [K]
- `k_o`, `E`: Arrhenius factor [m³ mol⁻¹ s⁻¹] and activation energy [J/mol]
- `β_area`: interfacial-area exponent [-]

# Returns
- λ [1/s]; zero without water, gas or reactive lime
"""
function reaction_decay_constant(C_g_total, C_lime, C_r, θ_w, θ_g, T, k_o, E, β_area)
    (θ_w > 0.0 && θ_g > 0.0 && T > 0.0) || return 0.0

    Δ_lime = C_lime - C_r
    Δ_lime > 0.0 || return 0.0

    R_gas = 8.3145  # J mol⁻¹ K⁻¹

    return (θ_w / θ_g) * arrhenius_coefficient(k_o, E, T) *
           interfacial_area_factor(C_g_total, T, β_area) *
           henry_solubility(T) * R_gas * T * Δ_lime
end

#=
IMPLEMENTATION NOTES:
=====================

Fully explicit FE solver for transient gas transport in porous media.

Diffusion, advection, gravitational flow and reaction kinetics (e.g., lime carbonation),
forward Euler time integration.

Heat from exothermic reactions included, assumed adiabatic (no heat loss to surroundings
or soil matrix).
=#

using Base.Threads
using Printf
using LinearAlgebra

"""
    ElementProperties

Time-invariant element properties, resolved once before the time loop to avoid
per-step `get_element_material` lookups.

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

k_o and E are reaction properties, read once from the reactant properties instead of per element.
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

Compute the effective thermal conductivity of the mixture, Eq. (fourier).

    λ_e = (1-n) λ_s + θ_w λ_w + θ_g λ_g

Volume-weighted (parallel) average, assembled like C_mix. Upper bound for a granular
assembly, since heat crosses grain contacts. λ_e is the calibration parameter if a
simulated thermal front is over-diffused.

λ_g is the plain average over declared species, without concentration weighting: λ_g
varies little between the gases and the gas term is a small fraction of λ_e.

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

Return, for each node, the highest-index element containing it.

Reaction and temperature terms are evaluated per node but use element-level properties
(θ_w, θ_g, Arrhenius pair). An element loop that assigns (instead of accumulating) leaves
each shared node with the value of the highest element index. Resolving ownership once
allows an order-independent, parallel node loop with the same result.

Nodes on a material interface get an arbitrary owner, consistent with the element-loop result.

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

Assemble the row-sum lumped mass vector. Mutates `M` [Nnodes].

# Formula
Eq. (lumped_mass_matrix): M_a = Σ_gp θ_g N_a(gp) w_gp W_gp, with `w_gp` the quadrature measure
of `ShapeFunctions.get_weight`: `|det J|` for plane strain, `r_gp |det J|` for axisymmetry,
Eq. (ax_lumped_mass_matrix).

The equal split `M_node = θ_g A_e / 4` matches only for a parallelogram (∫N_a dΩ = A_e/4).
Under axisymmetry `r_gp` varies within the element, therefore the node nearer the axis
receives less mass.
"""
function assemble_lumped_mass_vector!(M::Vector{Float64}, mesh, materials)
    fill!(M, 0.0)

    Nelements = mesh.num_elements

    # Shared nodes: scatter-add from a threaded element loop is a read-modify-write race.
    # Fixed element chunks, one buffer per chunk, summed in order after the join.
    #
    # Chunks indexed explicitly instead of Threads.threadid(): tasks may migrate between
    # threads, and threadid() may exceed Threads.nthreads() with interactive threads.
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

            # M_a = Σ_gp θ_g N_a(gp) w_gp W_gp. N_a weights the nodal share of the element measure
            for p in 1:4  # 4 Gauss points
                N_p = ShapeFunctions.get_N(p)
                dV = ShapeFunctions.get_weight(e, p) * ShapeFunctions.shape_funcs.gauss_weights[p]
                for i in 1:4
                    M_local[nodes[i]] += θ_g * N_p[i] * dV
                end
            end
        end
    end

    # Fixed-order reduction, reproducible M between runs
    for c in 1:nchunks
        M .+= M_chunk[c]
    end
end


"""
assemble_element_stiffness_matrices(mesh)

Assemble element geometric stiffness matrices without material properties [Nelements][4×4].

# Formula
K_e[i,j] = ∑_p (B · J^-1)^T · (B · J^-1) w_gp W_p

Material properties applied later. `w_gp` is the quadrature measure of
`ShapeFunctions.get_weight` and includes the radius for axisymmetric runs, therefore one
matrix serves both formulations: diffusive flux, Eq. (ax_diffusive_flux), and thermal
conductive flux, Eq. (ax_thermal_conductive_flux). Only the measure changes between the two.
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
            # Quadrature measure: |det J|, or r_gp |det J| when axisymmetric
            dV = ShapeFunctions.get_weight(e, p)

            # Precomputed physical-coordinate derivatives dN/dx = B · J^-1
            dN_dx = ShapeFunctions.get_dN_dx(e, p)  # [4 nodes, 2 coords]

            # Gauss weight
            w = ShapeFunctions.shape_funcs.gauss_weights[p]

            # K_e += w * w_gp * (dN/dx) * (dN/dx)^T, all ∇N_i · ∇N_j terms in one product
            K_e .+= (w * dV) .* (dN_dx * dN_dx')
        end
        
        # Store element matrix
        K_elements[e] = K_e
    end
    
    return K_elements
end


"""
    fully_explicit_diffusion_solver(mesh, materials, calc_params, time_data, log_print)

Solve transient gas transport in porous media with forward Euler time integration.

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
    global M_T, q_cond_T, q_adv_T, q_react_T, q_ext_T, q_heat
    global T_boundary, thermal_node_influences
    global boundary_node_influences
    global q_boundary

    # Track warnings for this step
    negative_conc_warned = Dict{Int, Bool}()  # Track warnings per gas

    # Universal gas constant [J/(mol·K)]
    R = 8.314
    M_caco3= 100.09 #g/mol
    ρ_caco3= 2.71e6 #g/m³
    
    # Reaction enthalpy of lime carbonation [J/mol CO2], from the material file.
    # Ca(OH)2(s) + CO2(g) -> CaCO3(s) + H2O(l)
    # Negative (IUPAC convention) for exothermic. Fallback for files without [reactants]
    # is DEFAULT_REACTION_ENTHALPY (read_materials.jl).
    ΔH_r = materials.reactants.reaction_enthalpy  # J/mol CO2

    # Carbonation kinetics, k_T = k_o exp(-E/RT). Reaction property, read once
    k_o_reaction = materials.reactants.arrhenius_factor    # m³/(mol·s)
    E_reaction = materials.reactants.activation_energy     # J/mol
    # Interfacial-area coefficient. Fitted jointly with k_o, must come from the same calibration
    β_area_reaction = materials.reactants.interfacial_area_beta  # dimensionless

    # Molar heat capacity per gas species [J/(mol·K)], gas_dictionary order.
    # Molar form of Eq. (heat_advection): (ρc)_g = Σ_i C_g^i c̄_p^i, no density conversion.
    #
    # Eq. (vol_heat_capacity) uses a mass basis and Eq. (heat_advection) a molar basis for
    # the same quantity; molar form used for both. Storage strictly requires c_v and
    # advection c_p, but the gas contributes ~1e-3 of the mixture capacity.
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
    calculate_stabilization = get(solver_settings, "advection_stabilization", 0) == 1

    # Without conduction or heat advection T is constant: isothermal run, ΔH_r ignored.
    # Stronger than the λ_e = 0, v_g = 0 limit of the manuscript (tex:474), which keeps the
    # local heat source. Adiabatic limit: enable conduction, set phase conductivities to zero.
    calculate_heat = calculate_heat_conduction || calculate_heat_advection

    # Legacy material files (single constant κ) parse with k_o = 0 and give no reaction
    if calculate_reaction
        log_print(@sprintf("   Reaction enthalpy ΔH_r: %.4g J/mol CO2", ΔH_r))
        log_print(@sprintf("   Arrhenius factor k_o: %.4g m³/(mol·s), activation energy E: %.4g J/mol",
                           k_o_reaction, E_reaction))
        if β_area_reaction == 0.0
            log_print("   Interfacial-area coefficient β: 0 (factor disabled, plain second-order law)")
        else
            log_print(@sprintf("   Interfacial-area coefficient β: %.4g [-] (a = 1 at atmospheric pressure)",
                               β_area_reaction))
        end
        if k_o_reaction <= 0.0
            log_print("Warning: reaction kinetics is enabled but lime_arrhenius_factor = 0, " *
                      "so no carbonation will occur. Check that the material file defines a " *
                      "[reactants] table with lime_arrhenius_factor and lime_activation_energy.")
        end

        # q̇ = ΔH_r dC_CO2/dt with dC_CO2/dt < 0, therefore heat release requires ΔH_r < 0.
        # A positive value cools the soil
        if ΔH_r > 0.0
            log_print("Warning: reaction_enthalpy = $(ΔH_r) J/mol is positive, i.e. " *
                      "endothermic, so carbonation will cool the soil. Lime carbonation " *
                      "is exothermic; enter the enthalpy as a negative number.")
        end
    end

    # Gravity vector from calc_params.
    #
    # Components are the physical acceleration (y up: gravity_y_component = -1.0).
    # Darcy's law, Eq. (Darcy_law):
    #
    #     v = -(K/μ_g) (∇p - ρ_g g),
    #
    # g enters the velocity with a plus sign; the gravitational flux has the leading minus
    # of Eq. (gravitational_flux). Reversing the input vector reverses the gravity drive.
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

    # Static element properties and node ownership, no dictionary lookups in the time loop
    elem_props = build_element_properties(mesh, materials, C_lime_residual)
    node_owner = build_node_owner_elements(mesh)

    # Element characteristic length for SU τ*. Static, computed once.
    # calculate_element_characteristic_length defined in time_step.jl (included earlier in kernel.jl)
    h_elem = calculate_stabilization ?
        [calculate_element_characteristic_length(mesh, e) for e in 1:Nelements] :
        Float64[]

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

    # Node probes sampled at data_saving_interval, independent of the VTK interval.
    # `nothing` when no nodes requested
    probe_writer = WriteProbes.init_node_probes(project_name, current_stage, mesh,
                                                materials.gas_dictionary,
                                                calc_params["probing_nodes"],
                                                Float64(data_saving_interval),
                                                calc_params["units"], current_time,
                                                log_print)

    # Initial (or restart) state as first probe sample
    WriteProbes.probe_sample!(probe_writer, current_time, C_g, dC_g_dt,
                              total_concentration, P, T, dT_dt, C_lime, C_caco3,
                              dC_lime_dt, degree_of_carbonation, binder_content, v;
                              force=true)

    # Flow vectors, allocated once, zeroed in place each step.
    # Function-local instead of the globals in initialize_flows.jl: type-stable, no global
    # lookup in the hot loop. Those globals are unused by this solver.
    q_diffusion = zeros(Float64, Nnodes, NGases)
    q_advection = zeros(Float64, Nnodes, NGases)
    q_gravitational = zeros(Float64, Nnodes, NGases)
    # Thermal contribution to the pressure gradient, Eq. (thermal_pressure_flux).
    # Separate from q_advection (same sign and flag) to report their relative magnitude (tex:263)
    q_thermal_press = zeros(Float64, Nnodes, NGases)
    total_rate = zeros(Float64, Nnodes)   # Boundary rate of change
    q_source_sink = zeros(Float64, Nnodes) # Only for CO2 for now

    # Composition weights for the Lagrangian pressure-BC correction. Per pressure-BC node:
    # free_conc_sum = level-n total of FREE gases (P_boundary-gated), n_free_nodes = free
    # gas count. Correction split in proportion to each free species concentration.
    # Fully recomputed each step from level-n C_g, no reset needed.
    free_conc_sum = zeros(Float64, Nnodes)
    n_free_nodes = zeros(Int, Nnodes)

    # Element-chunk decomposition for the threaded flux and velocity loops.
    #
    # nchunks independent of Threads.nthreads(): chunk-order reduction with a fixed chunk
    # count gives bitwise-reproducible results for any thread count.
    #
    # Memory: nchunks * Nnodes * (3*NGases + NDim) * 8 bytes. Mesh colouring would replace
    # the buffers for very large meshes.
    nchunks = max(1, min(Nelements, 32))
    chunk_bounds = round.(Int, range(0, Nelements, length = nchunks + 1))

    qd_chunk = [zeros(Float64, Nnodes, NGases) for _ in 1:nchunks]
    qa_chunk = [zeros(Float64, Nnodes, NGases) for _ in 1:nchunks]
    qg_chunk = [zeros(Float64, Nnodes, NGases) for _ in 1:nchunks]
    qT_chunk = [zeros(Float64, Nnodes, NGases) for _ in 1:nchunks]
    v_chunk = [zeros(Float64, Nnodes, NDim) for _ in 1:nchunks]

    # Per-element scratch, one set per chunk. Sized 4 for the nodes of a quad element.
    q_aux_chunk = [zeros(4) for _ in 1:nchunks]
    qT_aux_chunk = [zeros(4) for _ in 1:nchunks]
    ρ_g_chunk = [zeros(4) for _ in 1:nchunks]
    # SU stabilization scratch: dN_dx * v_gp and dN_dx * v_gp_T, one pair per chunk.
    # No per-Gauss-point allocation (same mul! convention as q_aux in the diffusion term)
    a_su_chunk = [zeros(4) for _ in 1:nchunks]
    aT_su_chunk = [zeros(4) for _ in 1:nchunks]

    # Energy equation: nodal heat source, per-chunk accumulators for lumped thermal
    # capacity and conductive/reactive fluxes
    Q_dot_T = zeros(Float64, Nnodes)
    MT_chunk = [zeros(Float64, Nnodes) for _ in 1:nchunks]
    qd_T_chunk = [zeros(Float64, Nnodes) for _ in 1:nchunks]
    qr_T_chunk = [zeros(Float64, Nnodes) for _ in 1:nchunks]
    qa_T_chunk = [zeros(Float64, Nnodes) for _ in 1:nchunks]
    # Heat flux projection: flux and normalising weight sum(N |detJ| W), Eq. (nodal_velocities)
    qh_chunk = [zeros(Float64, Nnodes, NDim) for _ in 1:nchunks]
    qhw_chunk = [zeros(Float64, Nnodes) for _ in 1:nchunks]
    heat_flux_weight = zeros(Float64, Nnodes)

    # Gauss-point Darcy velocity per element and Gauss point. Used by the nodal projection
    # and directly by the thermal advective flux, Eq. (thermal_advective_flux)
    v_gp_cache = zeros(Float64, Nelements, 4, NDim)

    #______________________________________________________
    # Time-dependent boundary conditions
    #
    # Flatten curve references once, only entries with a curve. Per-step cost scales with
    # the number of time-driven boundary nodes (zero for a constant model).
    # Each tuple also stores the mesh-file constant: scaled by a multiplier curve,
    # replaced by an absolute curve.
    #______________________________________________________
    time_functions = calc_params["time_functions"]

    transient_flow = [(node_id, g, ids[g], mesh.uniform_flow_bc[node_id][g])
                      for (node_id, ids) in mesh.uniform_flow_bc_tf
                      for g in 1:min(length(ids), NGases)
                      if ids[g] != 0 && haskey(mesh.uniform_flow_bc, node_id)]

    transient_conc = [(node_id, g, ids[g], mesh.concentration_bc[node_id][g])
                      for (node_id, ids) in mesh.concentration_bc_tf
                      for g in 1:min(length(ids), NGases)
                      if ids[g] != 0 && haskey(mesh.concentration_bc, node_id)]

    transient_press = [(node_id, id, mesh.absolute_pressure_bc[node_id])
                       for (node_id, id) in mesh.absolute_pressure_tf
                       if id != 0 && haskey(mesh.absolute_pressure_bc, node_id)]

    transient_temp = [(node_id, id, mesh.temperature_bc[node_id])
                      for (node_id, id) in mesh.temperature_bc_tf
                      if id != 0 && haskey(mesh.temperature_bc, node_id)]

    # No has_transient_* gate for partial_pressure_bc: re-imposed every step (follows T via
    # the ideal gas law), therefore the curve lookup is inline in that loop. Count used
    # for the log line only.
    n_transient_partial_press = sum(count(!=(0), ids) for ids in values(mesh.partial_pressure_bc_tf);
                                    init=0)

    has_transient_flow = !isempty(transient_flow)
    has_transient_conc = !isempty(transient_conc)
    has_transient_press = !isempty(transient_press)
    has_transient_temp = !isempty(transient_temp)

    if has_transient_flow || has_transient_conc || has_transient_press || has_transient_temp ||
       n_transient_partial_press > 0
        log_print("      Time-dependent boundaries: " *
                  "$(length(transient_flow)) flow, $(length(transient_conc)) concentration, " *
                  "$(length(transient_press)) pressure, $(length(transient_temp)) temperature, " *
                  "$(n_transient_partial_press) partial pressure")
    end

    # Main time stepping loop
    save_data = false
    # Reset at every output: one lime clamp warning per output interval. Initialized here
    # since the only other assignment is in the output block (undefined if lime is
    # exhausted before the first output).
    negative_lime_warned = false

    # Non-negativity clamp count per gas since the last output.
    #
    # The clamp warning appears once per gas per run, adequate for occasional round-off.
    # Thousands of clamps per step indicate an unstable mode reshaped into a bounded
    # standing oscillation (under-resolved stability limit). This went unnoticed on the
    # axisymmetric carbonation runs because a single warning at step 1 looks like round-off.
    #
    # Indexed by gas: update loop threaded over gases, one entry per task, no race, no lock.
    negative_conc_count = zeros(Int, NGases)

    # TEMPORARY DIAGNOSTIC: remove after the clamp-location investigation is resolved.
    # Per-node, per-gas clamp count to classify top offending nodes as
    # pressure-BC / concentration-BC / interior. One column per gas_idx, one thread per
    # column, no race (same as negative_conc_count).
    negative_conc_node_count = zeros(Int, Nnodes, NGases)

    # Step index of the last output, to normalise the clamp count by the number of nodal updates
    last_output_step = 0
    for step in 1:num_steps

        #______________________________________________________
        # Refresh the time-dependent boundary conditions for this step
        #
        # Forward Euler: state and rates at t_n, update at t_{n+1}. Evaluation times:
        #
        #   Neumann flux (enters the rate)        -> t_n
        #   Dirichlet value (imposed on new state) -> t_{n+1}
        #
        # dt adjusted at the end of the previous iteration to match the next output time,
        # therefore t_next is exact, no drift across output boundaries.
        #______________________________________________________
        t_next = current_time + dt

        if has_transient_flow
            # q_boundary prefilled at start-up, not zeroed here: overwrite curve-driven
            # entries only. dt passed for tables with average = true (step integral instead
            # of point sample), conserving injected mass when dt exceeds the table spacing.
            for (node_id, gas_idx, tf_id, q_const) in transient_flow
                q_boundary[node_id, gas_idx] =
                    bc_value(time_functions, tf_id, q_const, current_time, dt) *
                    flow_node_influences[node_id]
            end
        end

        if has_transient_press
            # Lagrange multiplier target set to the end-of-step pressure. T is still T_n here;
            # end-of-step re-imposition uses updated T. Difference O(dt dT/dt), below the
            # scheme error.
            for (node_id, tf_id, p_const) in transient_press
                total_concentration[node_id] =
                    bc_value(time_functions, tf_id, p_const, t_next) / (R * T[node_id])
            end
        end

        #reset flow vectors (q_boundary is prefilled and not reset here)
        q_diffusion .= 0.0
        q_advection .= 0.0
        q_gravitational .= 0.0
        q_thermal_press .= 0.0
        total_rate .= 0.0
        q_source_sink .= 0.0
        for c in 1:nchunks
            qd_chunk[c] .= 0.0
            qa_chunk[c] .= 0.0
            qg_chunk[c] .= 0.0
            qT_chunk[c] .= 0.0
        end

        #______________________________________________________
        #Reaction kinetics calculation start here
        #
        #Per node: rate depends on nodal state (C_g, C_lime, T) and element properties
        #resolved by `node_owner`. Same result as the element loop, order-independent and
        #parallel. Level-k state, must run before the gas fluxes use q_source_sink.
        #______________________________________________________
        if calculate_reaction && co2_gas_idx !== nothing
            @threads for node_id in 1:Nnodes
                e = node_owner[node_id]
                if e == 0
                    continue
                end
                props = elem_props[e]

                #Extent-of-reaction rate r >= 0, per unit volume of water
                # Total gas concentration for the area factor; C_aq uses the CO2 component
                C_g_tot = 0.0
                for g in 1:NGases
                    C_g_tot += C_g[node_id, g]
                end

                r = extent_of_reaction_rate(C_g[node_id, co2_gas_idx], C_g_tot,
                                            C_lime[node_id],
                                            props.residual_lime, props.θ_w, T[node_id],
                                            k_o_reaction, E_reaction, β_area_reaction)

                #Species rates from stoichiometric coefficients.
                #Lime per unit total volume: dA/dt = -θ_w r.
                dC_lime_dt[node_id] = -props.θ_w * r

                #CO2 per unit gas volume: sink scaled by θ_w/θ_g.
                if props.θ_g > 0.0
                    q_source_sink[node_id] = -M[node_id] * (props.θ_w / props.θ_g) * r
                else
                    q_source_sink[node_id] = 0.0
                end
            end
        end
        #______________________________________________________

        # Loop over element chunks, gas loop inside.
        #
        # Threading over elements instead of gases: parallel width = element count instead
        # of NGases (commonly 2). Per-chunk buffers avoid races on shared nodes, reduced in
        # chunk order. Fixed chunk count, bitwise-reproducible for any thread count.
        @threads for c in 1:nchunks
            qd_c = qd_chunk[c]
            qa_c = qa_chunk[c]
            qg_c = qg_chunk[c]
            qT_c = qT_chunk[c]
            q_aux = q_aux_chunk[c]
            qT_aux = qT_aux_chunk[c]
            ρ_g_buf = ρ_g_chunk[c]
            a_su = a_su_chunk[c]
            aT_su = aT_su_chunk[c]

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

                # Effective diffusivity of this species, for SU τ* only. Zero when diffusion is
                # off in solver_settings, consistent with the assembled system
                # (e.g. axisym_darcy_radial, diffusion = 0).
                D_g_eff = calculate_diffusion ? (θ_g * D_g * τ) : 0.0

                #Get current gas nodal concentrations
                C_e = [C_g[nodes[i], gas_idx] for i in 1:4]

                #______________________________________________________    
                #Diffusion calculation start here
                #______________________________________________________
                if calculate_diffusion
                    #Update diffusion flow vector ∑_p θ_g^p * D_g^p * τ^p * k_elm * det(J) * W_p

                    #local diffusion flow vector (per-chunk buffer)
                    mul!(q_aux, K_elements[e], C_e, θ_g * D_g * τ, 0.0)

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
                    fill!(qT_aux, 0.0)

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

                        # Quadrature measure: |det J|, or r_gp |det J| when axisymmetric
                        dV = ShapeFunctions.get_weight(e, p)

                        Wp= ShapeFunctions.shape_funcs.gauss_weights[p]

                        # Precomputed physical-coordinate derivatives dN/dx = B · J^-1
                        dN_dx = ShapeFunctions.get_dN_dx(e, p)  # [4 nodes, 2 coords]

                        #Update diffusion flow vector ∑_p K^p * T^p *C^p * k^p_elm *C_tot * w_gp * W_p / μ_g^p
                        #Eq. (advective_flux); Eq. (ax_advective_flux) when the measure includes r_gp
                        #Right-to-left contraction dN_dx * (dN_dx' * C_t), no 4×4 product
                        q_aux .+= (R * k_intrinsic * C_gp * T_gp * dV * Wp / μ_g) .* (dN_dx * (dN_dx' * C_t))

                        #Thermal contribution to the pressure gradient, Eq. (thermal_pressure_flux).
                        #Same operator as the advective flux with T and Σ_j C_g^j exchanged
                        #(differentiated T, Gauss-point total concentration). One extra contraction,
                        #no new element matrices (tex:368).
                        C_tot_gp = N_p' * C_t
                        qT_aux .+= (R * k_intrinsic * C_gp * C_tot_gp * dV * Wp / μ_g) .* (dN_dx * (dN_dx' * T_e))

                        #______________________________________________________
                        #Streamline-upwind (SU) stabilization of q_aux and qT_aux.
                        #Derivation in src/ADVECTION_STABILIZATION_NOTES.md. v_gp and v_gpT are
                        #the velocities of q_aux and qT_aux from the product-rule split
                        #grad(P) = R(T grad(C_tot) + C_tot grad(T)). Different vectors,
                        #therefore separate τ* values.
                        #______________________________________________________
                        if calculate_stabilization
                            v_gp = (-k_intrinsic * R * T_gp / μ_g) .* (dN_dx' * C_t)
                            v_gpT = (-k_intrinsic * R * C_tot_gp / μ_g) .* (dN_dx' * T_e)

                            τ_g = stab_tau(v_gp, h_elem[e], D_g_eff)
                            τ_gT = stab_tau(v_gpT, h_elem[e], D_g_eff)

                            mul!(a_su, dN_dx, v_gp)
                            q_aux .+= (τ_g * dV * Wp) .* (a_su .* (a_su' * C_e))

                            #Stabilizes the species concentration along the thermally driven
                            #streamline, therefore C_e contracted with v_gpT. Contracting T_e adds no
                            #diagonal weight to the C_g^i system (no discrete maximum principle) and
                            #has units of K m^3/s instead of a molar flux.
                            mul!(aT_su, dN_dx, v_gpT)
                            qT_aux .+= (τ_gT * dV * Wp) .* (aT_su .* (aT_su' * C_e))
                        end
                    end
                    for i in 1:4 #loop nodes in element
                        node_id = nodes[i] #global node id
                        qa_c[node_id, gas_idx] +=  q_aux[i]
                        qT_c[node_id, gas_idx] += qT_aux[i]
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

                        # Quadrature measure: |det J|, or r_gp |det J| when axisymmetric
                        dV = ShapeFunctions.get_weight(e, p)

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

                        #Gravitational flux, Eq. (gravitational_flux); Eq. (ax_gravitational_flux)
                        #when the measure includes r_gp. Axisymmetric gravity restricted to the axis,
                        #g = (0, -g) (Eq. ax_gravity), enforced at startup.
                        #ρ_g interpolated to a Gauss-point scalar: scaled 4-vector, no 4×4 outer product.
                        #
                        #Subtraction = leading minus of Eq. (gravitational_flux). g_vector is the
                        #user-oriented vector (downward for the usual mesh). The buffer is subtracted
                        #again in the rate equation with j_d, j_a and j_T, Eq. (FEM_compact). Both
                        #minuses give M Ċ_i = ∫ ∇N_i·(C v) with gravity velocity +(K/μ) ρ_g g, same
                        #convention as the advective term.
                        ρ_g_gp = N_p' * ρ_g
                        q_aux .-= (k_intrinsic * C_gp * dV * Wp * ρ_g_gp / μ_g) .* (dN_dx * g_vector)
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
            q_thermal_press .+= qT_chunk[c]
        end

        # Loop over all gases to form the rate of change
        @threads for gas_idx in 1:NGases
            # Get gas name for this species (needed for the reaction branch)
            gas_name = materials.gas_dictionary[gas_idx]

            # calculate rate of change dC/dt = q_net / M
            for i in 1:Nnodes
                dC_g_dt[i, gas_idx] = ((q_boundary[i, gas_idx] - q_diffusion[i, gas_idx] -
                                        q_advection[i, gas_idx] - q_thermal_press[i, gas_idx] -
                                        q_gravitational[i, gas_idx]) * P_boundary[i, gas_idx]) / M[i]

                if gas_name == "CO2" && calculate_reaction
                    #reaction contribution to the CO2 gas-phase rate [mol/(m³ gas·s)]
                    dC_g_dt_rxn = (q_source_sink[i] * P_boundary[i, gas_idx]) / M[i]

                    if C_g[i, gas_idx] + dt * (dC_g_dt[i, gas_idx] + dC_g_dt_rxn) < 0.0
                        #CO2 consumption limited to the available amount: reaction throttled
                        #to exhaust the gas at the end of the step, lime and heat rates
                        #scaled by the same factor.
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
                # total_concentration[j]: target for this step. Transient pressure: t_{n+1}
                # value set at the top of the step. Constant pressure: previous imposed value.

                #Composition-weighted correction. An equal absolute split subtracts the same
                #amount from every free species, which can make a minority species negative
                #(e.g. Air displaced by CO2 near an injection face; see advection_1d/case.toml,
                #same outlet mechanism).
                #λ_bc[j] = TOTAL correction rate (no /n_free), split at the application site in
                #proportion to each free gas level-n concentration via free_conc_sum[j].
                #Prescribed gases excluded from both the weighting sum and the split, otherwise
                #the prescribed gas departs from its boundary value.
                n_free = 0
                free_sum = 0.0
                for g in 1:NGases
                    n_free += P_boundary[j, g]
                    free_sum += P_boundary[j, g] * C_g[j, g]
                end
                free_conc_sum[j] = free_sum
                n_free_nodes[j] = n_free

                if n_free > 0
                    λ_bc[j] = M[j] * (C_total_j + dt * total_rate[j] - total_concentration[j]) / dt  # Total Lagrangian correction rate for node j
                else
                    λ_bc[j] = 0.0  # every gas is prescribed here, nothing left to correct
                end
            else
                λ_bc[j] = 0.0  # No influence length found, set multiplier to zero
                free_conc_sum[j] = 0.0
                n_free_nodes[j] = 0
            end
        end

        # Gas velocities before the concentration update: full right-hand side at level k,
        # v^k from p^k and C_g^k (same position as the manuscript algorithm listing). The
        # thermal advective flux uses v^k, without mixing old pressure and new concentrations.
        # Calculate nodal gas velocities using Darcy's law
        #zero the velocity vector
        v .= 0.0
        for c in 1:nchunks
            v_chunk[c] .= 0.0
        end

        # Separate pass from the flux loop: merging would change the state used for the velocity
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

                # Quadrature measure: |det J|, or r_gp |det J| when axisymmetric
                dV = ShapeFunctions.get_weight(e, p)

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

                #Calculate velocity at Gauss point using Darcy's law, Eq. (Darcy_law):
                #v = - (k/μ) (∇P - ρ_g g). Pressure part; gravity part added below if enabled.
                v_gp = - (k_intrinsic / μ_g_weighted) * grad_P

                # Gauss-point mass weight, used in numerator and denominator of the
                # point-to-grid transfer, Eq. (nodal_velocities). Including r_gp here satisfies
                # Eq. (ax_nodal_velocities) (r_gp in both sums) and the consistency requirement
                # of Appendix (app:velocity_transfer).
                mass_weight = props.θ_g * w * dV

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

                    #Gravity part of Eq. (Darcy_law): +(K/μ) ρ_g g, dense gas moves along g.
                    #This plus and the minus of the gravitational flux are one convention; changing
                    #one without the other makes reported velocity and transported mass inconsistent.
                    v_g_gp = (k_intrinsic / μ_g_weighted) * ρ_g_gp * g_vector

                    #Distribute mass-weighted gravitational velocity to nodes
                    for i in 1:4
                        node_id = nodes[i]
                        v_c[node_id, :] += v_g_gp * N_p[i] * mass_weight
                    end

                    #Total Gauss-point velocity, retained for the energy equation
                    v_gp = v_gp + v_g_gp
                end

                #Cache Gauss-point velocity (Eq. Darcy_Gauss) for the thermal advective flux
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

        #Loop gases again to apply Lagrangian correction and update concentrations
        @threads for gas_idx in 1:NGases
            # Get gas name for warnings
            gas_name = materials.gas_dictionary[gas_idx]
            
            # Update concentrations: C^(n+1) = C^n + dt * dC/dt
            for i in 1:Nnodes
                # Lagrangian correction only at nodes with boundary influence and only for
                # free gases. P_boundary zeroes dC_g_dt for a prescribed gas; without the same
                # factor here the prescribed concentration departs from its boundary value
                # each step. Visible only with reaction active (total_rate includes the CO2 sink).
                #
                # Composition-weighted split among free gases (see multiplier loop). C_TOL:
                # equal split when every free gas is ~0, avoiding 0/0.
                lagrangian_correction = 0.0
                if haskey(boundary_node_influences, i) && haskey(mesh.absolute_pressure_bc, i)
                    C_TOL = 1e-12
                    denom = free_conc_sum[i]
                    weight = denom > C_TOL ? (C_g[i, gas_idx] / denom) :
                             (n_free_nodes[i] > 0 ? 1.0 / n_free_nodes[i] : 0.0)
                    lagrangian_correction = (λ_bc[i] / M[i]) * weight * P_boundary[i, gas_idx]
                end

                C_g[i, gas_idx] += dt * (dC_g_dt[i, gas_idx] - lagrangian_correction)
                
                # Ensure non-negative and numerically stable concentrations
                C_MIN = 1e-12
                if C_g[i, gas_idx] < C_MIN
                    if C_g[i, gas_idx] < 0.0
                        negative_conc_count[gas_idx] += 1
                        negative_conc_node_count[i, gas_idx] += 1  # TEMPORARY DIAGNOSTIC
                        if !get(negative_conc_warned, gas_idx, false)
                            log_print("Warning: Negative concentration detected for gas $gas_name at step $step. Setting to zero.")
                            negative_conc_warned[gas_idx] = true
                        end
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
        

        # Update reaction kinetic terms for lime concentration
        if calculate_reaction
            for i in 1:Nnodes
                # Lime consumed this step. Clamp on the INCREMENT, such that CaCO3 gains exactly
                # the lime lost (1:1 stoichiometry, Eq. (carbonation)). Clamping C_lime only
                # would credit CaCO3 with nonexistent lime.
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
                #Caco3_max = 0 for lime-free material, guard against NaN
                degree_of_carbonation[i] = Caco3_max[i] > 0.0 ? C_caco3[i] / Caco3_max[i] : 0.0
            end
        end

        #______________________________________________________
        # Energy equation, Eq. (FEM_compact_energy):
        #
        #     M^L_T · dT/dt = q̃_ext - q̃_d - q̃_a + q̃_r
        #
        # Order of the manuscript algorithm listing: heat rate, lumped capacity, conductive
        # flux, reactive flux, boundary conditions, update. Runs with or without reaction;
        # reactive flux is zero without reaction.
        #______________________________________________________
        if calculate_heat
            ρ_w = materials.liquid.density
            c_w = materials.liquid.specific_heat

            # Nodal volumetric heat source, Q̇ = θ_w q̇, Eq. (heat_rate), total-volume basis:
            # q̇ = -ΔH_r r and dC_lime_dt = -θ_w r, therefore θ_w q̇ = ΔH_r dC_lime_dt.
            if calculate_reaction
                @threads for i in 1:Nnodes
                    Q_dot_T[i] = ΔH_r * dC_lime_dt[i]
                end
            else
                Q_dot_T .= 0.0
            end

            M_T .= 0.0
            q_cond_T .= 0.0
            q_adv_T .= 0.0
            q_react_T .= 0.0
            q_heat .= 0.0
            heat_flux_weight .= 0.0
            for c in 1:nchunks
                MT_chunk[c] .= 0.0
                qd_T_chunk[c] .= 0.0
                qa_T_chunk[c] .= 0.0
                qr_T_chunk[c] .= 0.0
                qh_chunk[c] .= 0.0
                qhw_chunk[c] .= 0.0
            end

            # Element assembly, same chunking as the gas fluxes
            @threads for c in 1:nchunks
                MT_c = MT_chunk[c]
                qd_c = qd_T_chunk[c]
                qa_c_T = qa_T_chunk[c]
                qr_c = qr_T_chunk[c]
                qh_c = qh_chunk[c]
                qhw_c = qhw_chunk[c]

                for e in (chunk_bounds[c] + 1):chunk_bounds[c + 1]
                    nodes = mesh.elements[e, :]
                    props = elem_props[e]
                    ρ_s = props.specific_gravity * ρ_w
                    solid_water = (1.0 - props.porosity) * ρ_s * props.specific_heat_solids +
                                  props.θ_w * ρ_w * c_w

                    T_e = [T[nodes[i]] for i in 1:4]

                    for p in 1:4
                        N_p = ShapeFunctions.shape_funcs.N[p]
                        Wp = ShapeFunctions.shape_funcs.gauss_weights[p]
                        # Quadrature measure: |det J| W, or r_gp |det J| W when axisymmetric.
                        # Used for the whole energy block, Eqs. (ax_lumped_capacity),
                        # (ax_thermal_conductive_flux), (ax_thermal_advective_flux),
                        # (ax_thermal_reactive_flux), and the heat-flux projection.
                        dV = ShapeFunctions.get_weight(e, p) * Wp

                        # (ρc)_g at the Gauss point, Eq. (heat_advection), per unit gas volume
                        # (weighted by θ_g in C_mix)
                        ρc_g_gp = 0.0
                        for g in 1:NGases
                            C_gp = 0.0
                            for i in 1:4
                                C_gp += N_p[i] * C_g[nodes[i], g]
                            end
                            ρc_g_gp += C_gp * c_p_molar[g]
                        end
                        C_mix_gp = solid_water + props.θ_g * ρc_g_gp

                        # Lumped thermal capacity, Eq. (lumped_capacity)
                        for i in 1:4
                            MT_c[nodes[i]] += C_mix_gp * N_p[i] * dV
                        end

                        # Conductive flux, Eq. (thermal_conductive_flux). Same contraction as
                        # the diffusive gas flux: cached dN_dx, right-to-left association.
                        if calculate_heat_conduction && props.λ_e > 0.0
                            dN_dx = ShapeFunctions.get_dN_dx(e, p)
                            gradT = dN_dx' * T_e
                            qd_gp = dN_dx * gradT
                            for i in 1:4
                                qd_c[nodes[i]] += props.λ_e * qd_gp[i] * dV
                            end
                        end

                        # Advective flux, Eq. (thermal_advective_flux), NON-CONSERVATIVE form
                        # (ρc)_g v·∇T instead of ∇·((ρc)_g v T).
                        #
                        # Consistent with the left-hand side C_mix ∂T/∂t (instead of
                        # ∂(C_mix T)/∂t). The two operators differ by T ∇·((ρc)_g v), linked by
                        # the gas mass balance to ∂(θ_g C)/∂t and the CO2 sink. The divergence form
                        # with this left-hand side adds a spurious source wherever gas accumulates
                        # or is consumed; on 8mm_size, 17% more energy than the reaction released.
                        #
                        # The divergence form also requires integration by parts, and the surface
                        # term ∮ N_i (ρc)_g T (v·n̂) dΓ was never assembled. Zero at interior nodes
                        # and walls (v·n̂ = 0), but at a flow-through face it acts as an enthalpy sink
                        # (ρc)_g |v| T proportional to ABSOLUTE temperature. Result on 8mm_size:
                        # inflow row near initial temperature, interior heating normally.
                        #
                        # v·∇T: no integration by parts, no surface term, and zero for a uniform
                        # temperature field at every node including boundaries.
                        #
                        # Sign: every nodal flux in Eq. (FEM_compact_energy) is subtracted, therefore
                        # positive v·∇T (temperature increasing along the flow) gives positive q̃_a
                        # and cools the node.
                        if calculate_heat_advection
                            ρc_gp = ρc_g_gp
                            dN_dx_a = ShapeFunctions.get_dN_dx(e, p)
                            gradT_a = dN_dx_a' * T_e          # [NDim]
                            vg = @view v_gp_cache[e, p, :]
                            v_gradT = 0.0
                            for d in 1:NDim
                                v_gradT += vg[d] * gradT_a[d]
                            end
                            for i in 1:4
                                qa_c_T[nodes[i]] += ρc_gp * v_gradT * N_p[i] * dV
                            end
                        end

                        # Total heat flux at this Gauss point, for output:
                        #   q_d = -lambda_e grad(T)          Eq. (fourier)
                        #   q_a = (rho c)_g T v_g            Eq. (heat_advection)
                        # Projected to the nodes below. Output diagnostic only; the solver uses
                        # the nodal flux vectors.
                        dN_dx_h = ShapeFunctions.get_dN_dx(e, p)
                        gradT_h = dN_dx_h' * T_e                # [NDim]
                        T_gp_h = 0.0
                        for i in 1:4
                            T_gp_h += N_p[i] * T_e[i]
                        end
                        for i in 1:4
                            nd = nodes[i]
                            wN = N_p[i] * dV
                            qhw_c[nd] += wN
                            for d in 1:NDim
                                q_flux = -props.λ_e * gradT_h[d]
                                if calculate_heat_advection
                                    q_flux += ρc_g_gp * T_gp_h * v_gp_cache[e, p, d]
                                end
                                qh_c[nd, d] += q_flux * wN
                            end
                        end

                        # Reactive flux, Eq. (thermal_reactive_flux)
                        if calculate_reaction
                            Q_gp = 0.0
                            for i in 1:4
                                Q_gp += N_p[i] * Q_dot_T[nodes[i]]
                            end
                            for i in 1:4
                                qr_c[nodes[i]] += Q_gp * N_p[i] * dV
                            end
                        end
                    end
                end
            end

            # Reduce in fixed chunk order
            for c in 1:nchunks
                M_T .+= MT_chunk[c]
                q_cond_T .+= qd_T_chunk[c]
                q_adv_T .+= qa_T_chunk[c]
                q_react_T .+= qr_T_chunk[c]
                q_heat .+= qh_chunk[c]
                heat_flux_weight .+= qhw_chunk[c]
            end

            # Normalise the projected flux by sum(N |detJ| W), Eq. (nodal_velocities):
            # element-uniform flux keeps its magnitude at the nodes (unlike tributary weighting).
            for i in 1:Nnodes
                if heat_flux_weight[i] > 0.0
                    for d in 1:NDim
                        q_heat[i, d] /= heat_flux_weight[i]
                    end
                end
            end

            # Convective (Robin) boundary, Eq. (conduction_flux_BC):
            #     q_d . n = h_c (T - T_inf)
            # Weak form term -∫(q_d.n) ω dS: a node above ambient loses heat. Edge integral
            # lumped with the tributary edge lengths from start-up (same as the pressure
            # boundary). h_c = 0 or unset boundary: adiabatic (natural condition).
            q_ext_T .= 0.0
            for (node_id, params) in mesh.convective_heat_bc
                h_c = params[1]
                T_inf = params[2]
                le = get(thermal_node_influences, node_id, 0.0)
                q_ext_T[node_id] -= h_c * (T[node_id] - T_inf) * le
            end

            # Update, Eq. (update_Temp)
            @threads for i in 1:Nnodes
                if M_T[i] > 0.0
                    dT_dt[i] = (q_ext_T[i] - q_cond_T[i] - q_adv_T[i] + q_react_T[i]) / M_T[i]
                else
                    dT_dt[i] = 0.0
                end
                # Prescribed node: RATE gated, therefore dT_dt = 0 there in the output
                dT_dt[i] *= T_boundary[i]
            end

            for i in 1:Nnodes
                T[i] += dt * dT_dt[i]

                # Ensure physically reasonable temperatures (above absolute zero)
                if T[i] < 0.0
                    log_print("Warning: Temperature below absolute zero at node $i. Setting to 0.0 K.")
                    T[i] = 0.0
                end
            end
        end

        # Curve-driven prescribed temperatures. The T_boundary gate zeroes the rate,
        # therefore the new value is written here. Outside calculate_heat on purpose: with
        # the energy equation off, a temperature history still drives the pressure through
        # the ideal gas law. Placed before any division by T, such that the rest of the
        # step uses the end-of-step temperature.
        if has_transient_temp
            for (node_id, tf_id, T_const) in transient_temp
                T[node_id] = bc_value(time_functions, tf_id, T_const, t_next)
            end
        end

        # Partial pressure BCs after the temperature update: C_g = P_partial / (R * T)
        #
        # Time-driven partial pressure (partial_pressure_bc_tf) evaluated at t_next, same as
        # every Dirichlet boundary. Runs every step (see n_transient_partial_press), curve
        # lookup inline.
        for (node_id, partial_pressures) in mesh.partial_pressure_bc
            tf_ids = get(mesh.partial_pressure_bc_tf, node_id, Int[])
            for gas_idx in 1:NGases
                # Recalculate concentration to maintain constant partial pressure
                # P_partial = C_g * R * T  =>  C_g = P_partial / (R * T)
                tf_id = gas_idx <= length(tf_ids) ? tf_ids[gas_idx] : 0
                p_partial = bc_value(time_functions, tf_id, partial_pressures[gas_idx], t_next)
                C_g[node_id, gas_idx] = p_partial / (R * T[node_id])
            end
        end

        # Curve-driven prescribed concentrations. Concentration Dirichlet enforced by rate
        # gating (P_boundary zeroes dC_g_dt), therefore a curve writes C_g explicitly.
        # Placed before the sum below to include the imposed values in the total.
        if has_transient_conc
            for (node_id, gas_idx, tf_id, C_const) in transient_conc
                C_g[node_id, gas_idx] = bc_value(time_functions, tf_id, C_const, t_next)
            end
        end

        # Calculate total gas concentrations after all gases are updated
        total_concentration = vec(sum(C_g, dims=2))

        #Apply total pressure boundary condition at nodes
        for node_id in keys(mesh.absolute_pressure_bc)
            # Absolute pressure BC via total concentration, at t_next (same time as the
            # Lagrange multiplier target)
            tf_id = get(mesh.absolute_pressure_tf, node_id, 0)
            P_bc = bc_value(time_functions, tf_id,
                            mesh.absolute_pressure_bc[node_id], t_next)
            total_concentration[node_id] = P_bc / (R * T[node_id])
        end

        #Update pressure using ideal gas law
        P= total_concentration .* R .* T  # Ideal gas law: P = C_total * R * T
        
        # Update current time
        current_time += dt

        # Probe sample after the pressure update (fully updated state), outside the VTK
        # block (independent sampling interval)
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

            # Clamp report. Rate normalised by nodes x gases x steps since the last output,
            # threshold independent of mesh size and output interval.
            total_clamped = sum(negative_conc_count)
            if total_clamped > 0
                steps_since = max(1, step - last_output_step)
                rate = total_clamped / (steps_since * Nnodes * NGases)
                per_gas = join((@sprintf("%s=%d", materials.gas_dictionary[g],
                                         negative_conc_count[g]) for g in 1:NGases), ", ")
                log_print(@sprintf("        Clamp: %d (%s)", total_clamped, per_gas))
                if rate > 1e-3
                    log_print(@sprintf("        ⚠ clamp %.2f%% of nodal updates - unstable mode, not round-off.", 100 * rate))
                    log_print("          Re-check the critical time step and lower the Courant number.")
                end

                # TEMPORARY DIAGNOSTIC: remove after the clamp-location investigation is
                # resolved. Hits on pressure/concentration-BC nodes indicate the Lagrangian
                # boundary correction; hits along the front do not.
                node_totals = vec(sum(negative_conc_node_count, dims=2))
                tops = String[]
                for nid in sortperm(node_totals, rev=true)[1:min(6, Nnodes)]
                    node_totals[nid] == 0 && break
                    f = String[]
                    haskey(mesh.absolute_pressure_bc, nid) && push!(f, "P")
                    haskey(mesh.concentration_bc, nid) && push!(f, "C")
                    haskey(mesh.uniform_flow_bc, nid) && push!(f, "flow")
                    isempty(f) && push!(f, "int")
                    push!(tops, @sprintf("%d:%d[%s]", nid, node_totals[nid], join(f, ",")))
                end
                isempty(tops) || log_print("          top: " * join(tops, " "))
                fill!(negative_conc_node_count, 0)
            end
            fill!(negative_conc_count, 0)
            last_output_step = step

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

Write the VTK output file for the current time step.

# Arguments
- `step`: Output file counter
- `time`: Current simulation time
- `project_name`: Project name for output files
- `total_concentration`: Total gas concentration vector
"""
function write_output_vtk(mesh, materials, step::Int, time::Float64, project_name, total_concentration)
    global C_g, P, T, v, dC_g_dt, dT_dt, dC_lime_dt
    global C_lime, C_caco3, degree_of_carbonation, binder_content, q_heat

    output_dir = "output"
    filename = joinpath(output_dir, project_name)
    
    # Prepare data for VTK output
    gas_names = materials.gas_dictionary    

    
    # Call VTK writer
    WriteVTK.write_vtk_file(filename, step, time, mesh, gas_names, (
        concentrations            = C_g,
        total_concentration       = total_concentration,
        absolute_pressure         = P,
        concentration_rates       = dC_g_dt,
        reaction_rates            = dC_lime_dt,
        lime_concentration        = C_lime,
        caco3_concentration       = C_caco3,
        degree_of_carbonation     = degree_of_carbonation,
        volumetric_binder_content = binder_content,
        gas_seepage_velocity      = v,
        temperature               = T,
        temperature_rate          = dT_dt,
        heat_flux                 = q_heat,
    ))
end
