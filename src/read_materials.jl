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

# Fields
- `name::String`: Soil name
- `specific_gravity::Float64`: Specific gravity of soil solids [-]
- `porosity::Float64`: Porosity [-]
- `saturation::Float64`: Degree of saturation [-]
- `granular_tortuosity::Float64`: Granular tortuosity [-]
- `intrinsic_permeability::Float64`: Intrinsic permeability [m²]
- `lime_content::Float64`: Lime content [-]
- `residual_lime::Float64`: Residual lime [-]
- `reaction_rate::Float64`: Chemical reaction rate [1/s]
- `specific_heat_solids::Float64`: Specific heat of solids [J/(kg·K)]
"""
mutable struct SoilProperties
    name::String
    specific_gravity::Float64
    porosity::Float64
    saturation::Float64
    granular_tortuosity::Float64
    intrinsic_permeability::Float64
    lime_content::Float64
    residual_lime::Float64
    reaction_rate::Float64
    specific_heat_solids::Float64
    
    function SoilProperties(name::String)
        new(name, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
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
"""
mutable struct MaterialData
    gas_dictionary::Vector{String}
    gases::Dict{String, GasProperties}
    liquid::LiquidProperties
    soil_dictionary::Vector{String}
    soils::Dict{String, SoilProperties}
    
    function MaterialData()
        new(String[], 
            Dict{String, GasProperties}(),
            LiquidProperties(),
            String[],
            Dict{String, SoilProperties}())
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
    materials.soil_dictionary = data["soil_dictionary"]
    
    # Parse soil properties
    parse_soil_properties!(materials, data["soil"])
    
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
        soil_props.lime_content = Float64(soil_info["lime_content"])
        soil_props.residual_lime = Float64(soil_info["residual_lime"])
        soil_props.reaction_rate = Float64(get(soil_info, "reaction_rate", 0.0))
        
        # Read thermal properties
        soil_props.specific_heat_solids = Float64(soil_info["specific_heat_solids"])
        
        materials.soils[soil_name] = soil_props
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
export read_materials_file, get_gas_properties, get_soil_properties
export get_num_gases, get_num_soils, get_liquid_properties
