#------------------------------------------------------------------------------
# ADSIM Materials Reader Module
# This module contains functions to read and parse .toml material files for 
# ADSIM FEM calculations
#------------------------------------------------------------------------------

using TOML

"""
GasProperties

Structure to store gas-specific properties.

# Fields
- `name::String`: Gas name
- `dynamic_viscosity::Float64`: Dynamic viscosity [Pa·s]
- `molar_mass::Float64`: Molar mass [g/mol]
- `diff_coefficient::Float64`: Diffusion coefficient [m²/s]
"""
mutable struct GasProperties
    name::String
    dynamic_viscosity::Float64
    molar_mass::Float64
    diff_coefficient::Float64
    
    function GasProperties(name::String)
        new(name, 0.0, 0.0, 0.0)
    end
end


"""
LiquidProperties

Structure to store liquid phase properties.

# Fields
- `dynamic_viscosity::Float64`: Dynamic viscosity [Pa·s]
- `density::Float64`: Density [kg/m³]
- `specific_heat::Float64`: Specific heat [J/(kg·K)]
"""
mutable struct LiquidProperties
    dynamic_viscosity::Float64
    density::Float64
    specific_heat::Float64
    
    function LiquidProperties()
        new(0.0, 0.0, 0.0)
    end
end


"""
SoilProperties

Structure to store soil material properties.

The carbonation reaction inputs (lime content, residual lime, and the Arrhenius
pair) used to live here. They now live in [`ReactionProperties`](@ref), keyed by
the same soil name, so that the soil block holds only physics that is
independent of whether reaction kinetics is switched on.

# Fields
- `name::String`: Soil name
- `specific_gravity::Float64`: Specific gravity of soil solids [-]
- `porosity::Float64`: Porosity [-]
- `saturation::Float64`: Degree of saturation [-]
- `granular_tortuosity::Float64`: Granular tortuosity [-]
- `intrinsic_permeability::Float64`: Intrinsic permeability [m²]
- `specific_heat_solids::Float64`: Specific heat of solids [J/(kg·K)]
"""
mutable struct SoilProperties
    name::String
    specific_gravity::Float64
    porosity::Float64
    saturation::Float64
    granular_tortuosity::Float64
    intrinsic_permeability::Float64
    specific_heat_solids::Float64

    function SoilProperties(name::String)
        new(name, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    end
end


"""
ReactionProperties

Structure to store the mix-design reaction inputs of one soil material.

These are the quantities of the reaction

    Ca(OH)2(s) + CO2(g) -> CaCO3(s) + H2O(l)

that depend on how much lime a given material carries, read from the
`[reaction."<soil name>"]` table of the material file. Each entry pairs with the
soil of the same name, so different layers can carry different lime contents.

The kinetics of the reaction (k_o, E) are a property of the chemistry rather
than of a soil layer and live in [`ReactantProperties`](@ref) instead.

# Fields
- `name::String`: Soil name this reaction data belongs to
- `lime_content::Float64`: Lime content β_l, as a fraction of the weight of solids [-]
- `residual_lime::Float64`: Residual (non-reactive) lime, as a fraction of the lime content [-]
"""
mutable struct ReactionProperties
    name::String
    lime_content::Float64
    residual_lime::Float64

    function ReactionProperties(name::String)
        new(name, 0.0, 0.0)
    end
end


# Reaction enthalpy of lime carbonation used when the material file predates the
# `[reactants]` table, so that older inputs reproduce their previous results.
#
# Ca(OH)2(s) + CO2(g) -> CaCO3(s) + H2O(l). Hess's law on the standard enthalpies
# of formation at 298.15 K (-986.09, -393.51, -1207.6, -285.83 kJ/mol) gives
# -113.8 kJ/mol CO2. Negative under the IUPAC convention, i.e. exothermic. The
# liquid-phase product water applies here; taking H2O(g) instead would give
# -69.8 kJ/mol.
const DEFAULT_REACTION_ENTHALPY = -113800.0  # J/mol CO2


"""
ReactantProperties

Structure to store the reactant (chemistry) properties shared by every material.

These describe the carbonation reaction itself, so they do not vary from one
soil layer to the next. How much lime each layer carries is a mix-design input
and lives in [`ReactionProperties`](@ref) instead.

# Fields
- `reaction_enthalpy::Float64`: Enthalpy of the carbonation reaction ΔH_r [J/mol CO2].
  Negative is exothermic, following the IUPAC convention.
- `arrhenius_factor::Float64`: Arrhenius factor k_o of the rate coefficient
  k_T = k_o exp(-E/RT) [m³/(mol·s)]
- `activation_energy::Float64`: Activation energy E of the carbonation reaction [J/mol]
"""
mutable struct ReactantProperties
    reaction_enthalpy::Float64
    arrhenius_factor::Float64
    activation_energy::Float64

    function ReactantProperties()
        new(DEFAULT_REACTION_ENTHALPY, 0.0, 0.0)
    end
end


"""
MaterialData

Structure to store all material properties including gases, liquids, and soils.

# Fields
- `gas_dictionary::Vector{String}`: List of gas names
- `gases::Dict{String, GasProperties}`: Gas properties indexed by name
- `liquid::LiquidProperties`: Liquid phase properties
- `soil_dictionary::Vector{String}`: List of soil names
- `soils::Dict{String, SoilProperties}`: Soil properties indexed by name
- `reactions::Dict{String, ReactionProperties}`: Reaction inputs indexed by soil name
- `reactants::ReactantProperties`: Reactant properties shared by all materials
"""
mutable struct MaterialData
    gas_dictionary::Vector{String}
    gases::Dict{String, GasProperties}
    liquid::LiquidProperties
    soil_dictionary::Vector{String}
    soils::Dict{String, SoilProperties}
    reactions::Dict{String, ReactionProperties}
    reactants::ReactantProperties

    function MaterialData()
        new(String[],
            Dict{String, GasProperties}(),
            LiquidProperties(),
            String[],
            Dict{String, SoilProperties}(),
            Dict{String, ReactionProperties}(),
            ReactantProperties())
    end
end


"""
read_materials_file(filename::String) -> MaterialData

Read an ADSIM materials .toml file and return a MaterialData structure containing
all material properties.

# Arguments
- `filename::String`: Path to the .toml materials file

# Returns
- `MaterialData`: Structure containing all material properties

# Example
```julia
materials = read_materials_file("materials.toml")
println("Number of gases: ", length(materials.gas_dictionary))
println("Number of soils: ", length(materials.soil_dictionary))
```
"""
function read_materials_file(filename::String)
    materials = MaterialData()
    
    # Read TOML file
    data = TOML.parsefile(filename)
    
    # Parse gas dictionary
    materials.gas_dictionary = data["gas_dictionary_"]
    
    # Parse gas properties
    parse_gas_properties!(materials, data["gas"])
    
    # Parse liquid properties
    parse_liquid_properties!(materials, data["liquid"])
    
    # Parse soil dictionary
    materials.soil_dictionary = data["soil_dictionary_"]

    # Parse soil properties
    parse_soil_properties!(materials, data["soil"])

    # Parse the reaction tables. Both fall back to the soil tables for material
    # files written before the reaction inputs were split out of the soil block.
    reaction_data = get(data, "reaction", Dict{String, Any}())
    reactant_data = get(data, "reactants", Dict{String, Any}())

    parse_reaction_properties!(materials, reaction_data, data["soil"], filename)
    parse_reactant_properties!(materials, reactant_data, reaction_data, data["soil"], filename)

    return materials
end


"""
parse_gas_properties!(materials::MaterialData, gas_data::Dict)

Parse gas properties from TOML data and store in MaterialData structure.

# Arguments
- `materials::MaterialData`: Material data structure to populate
- `gas_data::Dict`: Dictionary containing gas property data from TOML
"""
function parse_gas_properties!(materials::MaterialData, gas_data::Dict)
    for gas_name in materials.gas_dictionary
        gas_props = GasProperties(gas_name)
        gas_info = gas_data[gas_name]
        
        # Read gas properties
        gas_props.dynamic_viscosity = Float64(gas_info["dynamic_viscosity"])
        gas_props.molar_mass = Float64(gas_info["molar_mass"])
        gas_props.diff_coefficient = Float64(gas_info["diff_coefficient"])
        
        materials.gases[gas_name] = gas_props
    end
end


"""
parse_liquid_properties!(materials::MaterialData, liquid_data::Dict)

Parse liquid properties from TOML data and store in MaterialData structure.

# Arguments
- `materials::MaterialData`: Material data structure to populate
- `liquid_data::Dict`: Dictionary containing liquid property data from TOML
"""
function parse_liquid_properties!(materials::MaterialData, liquid_data::Dict)
    # Read liquid properties
    materials.liquid.dynamic_viscosity = Float64(liquid_data["dynamic_viscosity"])
    materials.liquid.density = Float64(liquid_data["density"])
    materials.liquid.specific_heat = Float64(liquid_data["specific_heat"])
end


"""
parse_soil_properties!(materials::MaterialData, soil_data::Dict)

Parse soil properties from TOML data and store in MaterialData structure.

# Arguments
- `materials::MaterialData`: Material data structure to populate
- `soil_data::Dict`: Dictionary containing soil property data from TOML
"""
function parse_soil_properties!(materials::MaterialData, soil_data::Dict)
    for soil_name in materials.soil_dictionary
        soil_props = SoilProperties(soil_name)
        soil_info = soil_data[soil_name]
        
        # Read physical properties
        soil_props.specific_gravity = Float64(soil_info["specific_gravity"])
        soil_props.porosity = Float64(soil_info["porosity"])
        soil_props.saturation = Float64(soil_info["saturation"])
        soil_props.granular_tortuosity = Float64(soil_info["granular_tortuosity"])
        soil_props.intrinsic_permeability = Float64(soil_info["intrinsic_permeability"])

        # Read thermal properties
        soil_props.specific_heat_solids = Float64(soil_info["specific_heat_solids"])

        materials.soils[soil_name] = soil_props
    end
end


"""
parse_reactant_properties!(materials::MaterialData, reactant_data::Dict,
                           reaction_data::Dict, soil_data::Dict, filename::String)

Parse the reactant properties shared by all materials from TOML data.

The Arrhenius pair describes the reaction rather than a soil layer, so it is read
from `[reactants]`. Material files that predate that grouping carry it per
material, in either `[reaction."<soil name>"]` or `[soil."<soil name>"]`; those
are still honoured. Because one reaction now serves every material, a legacy file
whose materials disagree on the kinetics cannot be represented exactly: the first
material's values are taken and the conflict is reported.

# Arguments
- `materials::MaterialData`: Material data structure to populate
- `reactant_data::Dict`: Contents of the `[reactants]` table, empty if absent
- `reaction_data::Dict`: Contents of the `[reaction]` tables, used as a fallback
- `soil_data::Dict`: Contents of the `[soil]` tables, used as the oldest fallback
- `filename::String`: Material file path, quoted in the warnings
"""
function parse_reactant_properties!(materials::MaterialData, reactant_data::Dict,
                                    reaction_data::Dict, soil_data::Dict,
                                    filename::String)
    materials.reactants.reaction_enthalpy =
        Float64(get(reactant_data, "reaction_enthalpy", DEFAULT_REACTION_ENTHALPY))

    if haskey(reactant_data, "lime_arrhenius_factor") ||
       haskey(reactant_data, "lime_activation_energy")
        materials.reactants.arrhenius_factor =
            Float64(get(reactant_data, "lime_arrhenius_factor", 0.0))
        materials.reactants.activation_energy =
            Float64(get(reactant_data, "lime_activation_energy", 0.0))
        return
    end

    # Legacy layout: gather the kinetics from wherever each material declared them.
    found = Tuple{String, Float64, Float64}[]
    for soil_name in materials.soil_dictionary
        for source in (get(reaction_data, soil_name, nothing),
                       get(soil_data, soil_name, nothing))
            source === nothing && continue
            if haskey(source, "lime_arrhenius_factor") ||
               haskey(source, "lime_activation_energy")
                push!(found, (soil_name,
                              Float64(get(source, "lime_arrhenius_factor", 0.0)),
                              Float64(get(source, "lime_activation_energy", 0.0))))
                break
            end
        end
    end

    isempty(found) && return

    _, k_o, E = first(found)
    materials.reactants.arrhenius_factor = k_o
    materials.reactants.activation_energy = E

    conflicting = [f for f in found if f[2] != k_o || f[3] != E]
    if !isempty(conflicting)
        @warn """
        Materials in $(filename) declare different carbonation kinetics, but the
        Arrhenius pair is now a single property of the reaction. The values of
        '$(first(found)[1])' are used: lime_arrhenius_factor = $(k_o),
        lime_activation_energy = $(E).

        Ignored: $(join(["$(c[1]) (k_o = $(c[2]), E = $(c[3]))" for c in conflicting], ", "))

        Set the pair explicitly in a [reactants] table to silence this.
        """
    end
end


"""
parse_reaction_properties!(materials::MaterialData, reaction_data::Dict,
                           soil_data::Dict, filename::String)

Parse the per-material mix-design reaction inputs from TOML data.

Reads the `[reaction."<soil name>"]` table of each soil. Material files written
before the reaction inputs were split out of the soil block carry these keys
inside `[soil."<soil name>"]` instead; those are still honoured, with a warning,
so that existing input files keep running unchanged.

The kinetics are read separately by [`parse_reactant_properties!`](@ref), since
they belong to the reaction rather than to a material.

# Arguments
- `materials::MaterialData`: Material data structure to populate
- `reaction_data::Dict`: Dictionary containing reaction data from TOML, empty if absent
- `soil_data::Dict`: Dictionary containing soil data from TOML, used as the legacy fallback
- `filename::String`: Material file path, quoted in the deprecation warning
"""
function parse_reaction_properties!(materials::MaterialData, reaction_data::Dict,
                                    soil_data::Dict, filename::String)
    legacy_soils = String[]

    for soil_name in materials.soil_dictionary
        reaction_props = ReactionProperties(soil_name)

        if haskey(reaction_data, soil_name)
            reaction_info = reaction_data[soil_name]
        else
            # Legacy layout: the reaction keys sit in the soil table.
            reaction_info = get(soil_data, soil_name, Dict{String, Any}())
            if any(haskey(reaction_info, k) for k in ("lime_content", "residual_lime",
                                                      "lime_arrhenius_factor",
                                                      "lime_activation_energy"))
                push!(legacy_soils, soil_name)
            end
        end

        reaction_props.lime_content = Float64(get(reaction_info, "lime_content", 0.0))
        reaction_props.residual_lime = Float64(get(reaction_info, "residual_lime", 0.0))

        materials.reactions[soil_name] = reaction_props
    end

    if !isempty(legacy_soils)
        @warn """
        Reaction inputs were read from the [soil] tables of $(filename).
        These now belong in their own tables and should be moved:

            [reactants]
            reaction_enthalpy = ...
            lime_arrhenius_factor = ...
            lime_activation_energy = ...

            [reaction."$(first(legacy_soils))"]
            lime_content = ...
            residual_lime = ...

        Without a [reactants] table the built-in enthalpy of
        $(DEFAULT_REACTION_ENTHALPY) J/mol CO2 is used.
        Affected materials: $(join(legacy_soils, ", "))
        """
    end
end


"""
get_gas_properties(materials::MaterialData, gas_name::String) -> Union{GasProperties, Nothing}

Get the properties for a specific gas by name.

# Arguments
- `materials::MaterialData`: Material data structure
- `gas_name::String`: Name of the gas

# Returns
- `Union{GasProperties, Nothing}`: Gas properties or nothing if not found

# Example
```julia
co2 = get_gas_properties(materials, "CO2")
if co2 !== nothing
    println("CO2 viscosity: ", co2.dynamic_viscosity)
end
```
"""
function get_gas_properties(materials::MaterialData, gas_name::String)
    return get(materials.gases, gas_name, nothing)
end


"""
get_soil_properties(materials::MaterialData, soil_name::String) -> Union{SoilProperties, Nothing}

Get the properties for a specific soil by name.

# Arguments
- `materials::MaterialData`: Material data structure
- `soil_name::String`: Name of the soil

# Returns
- `Union{SoilProperties, Nothing}`: Soil properties or nothing if not found

# Example
```julia
soil = get_soil_properties(materials, "Soil 1")
if soil !== nothing
    println("Porosity: ", soil.porosity)
end
```
"""
function get_soil_properties(materials::MaterialData, soil_name::String)
    return get(materials.soils, soil_name, nothing)
end


"""
get_reaction_properties(materials::MaterialData, soil_name::String) -> Union{ReactionProperties, Nothing}

Get the carbonation reaction inputs for a specific soil material by name.

# Arguments
- `materials::MaterialData`: Material data structure
- `soil_name::String`: Name of the soil

# Returns
- `Union{ReactionProperties, Nothing}`: Reaction properties or nothing if not found

# Example
```julia
reaction = get_reaction_properties(materials, "Soil 1")
if reaction !== nothing
    println("Lime content: ", reaction.lime_content)
end
```
"""
function get_reaction_properties(materials::MaterialData, soil_name::String)
    return get(materials.reactions, soil_name, nothing)
end


"""
get_reactant_properties(materials::MaterialData) -> ReactantProperties

Get the reactant properties shared by all materials.

# Arguments
- `materials::MaterialData`: Material data structure

# Returns
- `ReactantProperties`: Reactant properties

# Example
```julia
reactants = get_reactant_properties(materials)
println("Reaction enthalpy: ", reactants.reaction_enthalpy)
```
"""
function get_reactant_properties(materials::MaterialData)
    return materials.reactants
end


"""
get_num_gases(materials::MaterialData) -> Int

Get the total number of gases defined in the material data.

# Arguments
- `materials::MaterialData`: Material data structure

# Returns
- `Int`: Number of gases
"""
function get_num_gases(materials::MaterialData)
    return length(materials.gas_dictionary)
end


"""
get_num_soils(materials::MaterialData) -> Int

Get the total number of soils defined in the material data.

# Arguments
- `materials::MaterialData`: Material data structure

# Returns
- `Int`: Number of soils
"""
function get_num_soils(materials::MaterialData)
    return length(materials.soil_dictionary)
end


"""
get_liquid_properties(materials::MaterialData) -> LiquidProperties

Get the liquid phase properties.

# Arguments
- `materials::MaterialData`: Material data structure

# Returns
- `LiquidProperties`: Liquid properties

# Example
```julia
liquid = get_liquid_properties(materials)
println("Liquid density: ", liquid.density)
```
"""
function get_liquid_properties(materials::MaterialData)
    return materials.liquid
end


# Export all public functions and types
export MaterialData, GasProperties, LiquidProperties, SoilProperties
export ReactionProperties, ReactantProperties
export read_materials_file, get_gas_properties, get_soil_properties
export get_num_gases, get_num_soils, get_liquid_properties
export get_reaction_properties, get_reactant_properties
