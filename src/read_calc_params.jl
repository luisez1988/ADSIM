# ADSIM - Read Calculation Parameters from TOML
# This module contains functions to read and parse calculation parameters

using TOML

"""
    read_calc_params(filename::String)

Read calculation parameters from a TOML file and return a dictionary with all parameters.

# Arguments
- `filename::String`: Path to the calculation parameters TOML file

# Returns
- Dictionary containing all calculation parameters organized by sections
"""
function read_calc_params(filename::String)
    calc_data = TOML.parsefile(filename)
    return calc_data
end

"""
    get_units(calc_data::Dict)

Extract unit system from calculation parameters.

# Arguments
- `calc_data::Dict`: Dictionary containing calculation parameters

# Returns
- Dictionary with unit definitions
"""
function get_units(calc_data::Dict)
    units = calc_data["units"]
    return Dict(
        "geometry_unit" => units["geometry_unit"],
        "mass_unit" => units["mass_unit"],
        "pressure_unit" => units["pressure_unit"],
        "temperature_unit" => units["temperature_unit"],
        "time_unit" => units["time_unit"]
    )
end

"""
    get_gravity(calc_data::Dict)

Extract gravity parameters from calculation data.

# Arguments
- `calc_data::Dict`: Dictionary containing calculation parameters

# Returns
- Dictionary with gravity magnitude and components
"""
function get_gravity(calc_data::Dict)
    gravity = calc_data["gravity"]
    return Dict(
        "magnitude" => gravity["gravity_magnitude"],
        "x_component" => gravity["gravity_x_component"],
        "y_component" => gravity["gravity_y_component"]
    )
end

"""
    get_solver_type(calc_data::Dict)

Extract solver type from calculation data.

# Arguments
- `calc_data::Dict`: Dictionary containing calculation parameters

# Returns
- String with solver type
"""
function get_solver_type(calc_data::Dict)
    return calc_data["solver"]["solver_type"]
end

"""
    get_time_stepping(calc_data::Dict)

Extract time stepping parameters from calculation data.

# Arguments
- `calc_data::Dict`: Dictionary containing calculation parameters

# Returns
- Dictionary with time stepping parameters
"""
function get_time_stepping(calc_data::Dict)
    time_data = calc_data["time_stepping"]
    return Dict(
        "total_simulation_time" => time_data["total_simulation_time"],
        "time_per_step" => time_data["time_per_step"],
        "courant_number" => time_data["courant_number"]
    )
end

"""
    get_data_saving_interval(calc_data::Dict)

Extract data saving interval from calculation data.

# Arguments
- `calc_data::Dict`: Dictionary containing calculation parameters

# Returns
- Float with data saving interval
"""
function get_data_saving_interval(calc_data::Dict)
    return calc_data["data_saving"]["data_saving_interval"]
end

"""
    get_probing_nodes(calc_data::Dict)

Extract probing node information from calculation data.

# Arguments
- `calc_data::Dict`: Dictionary containing calculation parameters

# Returns
- Dictionary with number of nodes and array of node IDs to probe
"""
function get_probing_nodes(calc_data::Dict)
    probing = calc_data["probing"]
    return Dict(
        "number_of_nodes" => probing["number_of_nodes"],
        "nodes_to_probe" => probing["nodes_to_probe"]
    )
end

"""
    get_probing_elements(calc_data::Dict)

Extract probing element information from calculation data.

# Arguments
- `calc_data::Dict`: Dictionary containing calculation parameters

# Returns
- Dictionary with number of elements and array of element IDs to probe
"""
function get_probing_elements(calc_data::Dict)
    probing = calc_data["probing"]
    return Dict(
        "number_of_elements" => probing["number_of_elements"],
        "elements_to_probe" => probing["elements_to_probe"]
    )
end

"""
    get_all_calc_params(filename::String)

Read calculation parameters from a TOML file and return a structured dictionary
with all parsed parameters organized by category.

# Arguments
- `filename::String`: Path to the calculation parameters TOML file

# Returns
- Dictionary containing all calculation parameters in structured format
"""
function get_all_calc_params(filename::String)
    calc_data = read_calc_params(filename)
    
    return Dict(
        "units" => get_units(calc_data),
        "gravity" => get_gravity(calc_data),
        "solver_type" => get_solver_type(calc_data),
        "time_stepping" => get_time_stepping(calc_data),
        "data_saving_interval" => get_data_saving_interval(calc_data),
        "probing_nodes" => get_probing_nodes(calc_data),
        "probing_elements" => get_probing_elements(calc_data)
    )
end
