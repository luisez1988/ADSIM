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
    get_solver_settings(calc_data::Dict)

Extract solver settings from calculation data including dimension and calculation mode flags.

# Arguments
- `calc_data::Dict`: Dictionary containing calculation parameters

# Returns
- Dictionary with solver dimension and calculation mode components
"""
function get_solver_settings(calc_data::Dict)
    solver = calc_data["solver"]

    # The dimension string used to be carried through and only printed, so a typo would run
    # silently as plane strain. Validate it: the two formulations differ in physics, not in
    # presentation, and selecting the wrong one must not be a quiet outcome.
    dimension = solver["solver_type"]
    if !(dimension in ("2D-Plane", "2D-Axisymmetric"))
        error("Unknown solver_type \"$dimension\" in [solver]. " *
              "Expected \"2D-Plane\" or \"2D-Axisymmetric\".")
    end

    return Dict(
        "dimension" => dimension,
        # Axisymmetric analysis per Appendix (app:axisymmetric): the first coordinate is the
        # radius r, the second the axial coordinate z, and every quadrature measure carries
        # the factor r_gp.
        "axisymmetric" => dimension == "2D-Axisymmetric",
        "diffusion" => solver["diffusion"],
        "advection" => solver["advection"],
        "gravity" => solver["gravity"],
        "reaction_kinetics" => solver["reaction_kinetics"],
        # Heat transport. Optional so calculation files written before the energy
        # equation still load, and off by default. With neither mechanism active the
        # temperature field does not evolve at all and the run is isothermal, whatever
        # reaction enthalpy the material file declares.
        "heat_conduction" => get(solver, "heat_conduction", 0),
        "heat_advection" => get(solver, "heat_advection", 0)
    )
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
    # Probing is optional: a calculation file written without a [probing] section,
    # or with the section but no node list, simply runs without probes
    probing = get(calc_data, "probing", Dict{String, Any}())
    nodes = Int[Int(n) for n in get(probing, "nodes_to_probe", Int[])]
    return Dict(
        "number_of_nodes" => Int(get(probing, "number_of_nodes", length(nodes))),
        "nodes_to_probe" => nodes
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
    # Optional for the same reason as the node list. Element probing is parsed but
    # not yet consumed by the solver; only nodes are written as time series.
    probing = get(calc_data, "probing", Dict{String, Any}())
    elements = Int[Int(e) for e in get(probing, "elements_to_probe", Int[])]
    return Dict(
        "number_of_elements" => Int(get(probing, "number_of_elements", length(elements))),
        "elements_to_probe" => elements
    )
end

"""
    get_time_functions_section(calc_data::Dict, data_dir::AbstractString)

Build the time function table from the `[[time_function]]` section.

Optional for the same reason as the probing section: a calculation file written
without it yields an empty table and every boundary condition stays constant.
`data_dir` is the folder holding the calculation file, and any CSV referenced by
a table curve is resolved against it so the data travels with the input set.

# Returns
- `Dict{Int, TimeFunction}` keyed by the id referenced from the mesh file
"""
function get_time_functions_section(calc_data::Dict, data_dir::AbstractString)
    return read_time_functions(calc_data, data_dir)
end

"""
    log_analysis_type(solver_settings::Dict)

Generate analysis type log message based on enabled solver components.

# Arguments
- `solver_settings::Dict`: Dictionary containing solver settings

# Returns
- String with formatted analysis type information
"""
function log_analysis_type(solver_settings::Dict)
    msg = "   ✓ Dimension: $(solver_settings["dimension"])\n"
    if get(solver_settings, "axisymmetric", false)
        msg *= "   ✓ Geometry: axisymmetric - x is the radius r, y the axis of revolution z;\n" *
               "               every quadrature measure carries the factor r_gp\n"
    end

    # Build component list (values are 0 or 1)
    components = String[]
    solver_settings["diffusion"] == 1 && push!(components, "Diffusion")
    solver_settings["advection"] == 1 && push!(components, "Advection")
    solver_settings["gravity"] == 1 && push!(components, "Gravity")
    solver_settings["reaction_kinetics"] == 1 && push!(components, "Reaction Kinetics")
    solver_settings["heat_conduction"] == 1 && push!(components, "Heat Conduction")
    solver_settings["heat_advection"] == 1 && push!(components, "Heat Advection")

    if isempty(components)
        msg *= "   ✓ Solver: WARNING - No components selected!"
    else
        msg *= "   ✓ Solver: $(join(components, " + "))"
    end

    # Say so explicitly rather than leaving the reader to infer it from two absent
    # entries in the list above: with no heat transport the temperature is frozen,
    # so a nonzero reaction enthalpy is deliberately ignored.
    if solver_settings["reaction_kinetics"] == 1 &&
       solver_settings["heat_conduction"] == 0 && solver_settings["heat_advection"] == 0
        msg *= "\n   ✓ Thermal: isothermal (no heat transport active, reaction enthalpy ignored)"
    end

    return msg
end

"""
    validate_reaction_kinetics_requirements(solver_settings::Dict, materials)

Validate that CO2 gas is present in materials when reaction_kinetics is enabled.

# Arguments
- `solver_settings::Dict`: Dictionary containing solver settings
- `materials`: MaterialData structure containing gas definitions

# Throws
- `ErrorException`: If reaction_kinetics is enabled but CO2 gas is not defined

# Returns
- `true` if validation passes
"""
function validate_reaction_kinetics_requirements(solver_settings::Dict, materials)
    # Check if reaction kinetics is enabled
    if solver_settings["reaction_kinetics"] == 1
        # Check if CO2 exists in the gas dictionary
        if !("CO2" in materials.gas_dictionary)
            error("""
            Validation Error: Reaction kinetics is enabled but CO2 gas is not defined.
            
            When using reaction_kinetics = 1 in the solver settings, you must include 
            CO2 in the gas_dictionary_ and define its properties in the materials file.
            
            Example materials file entry:
            gas_dictionary_ = ["CO2"]
            
            [gas."CO2"]
            dynamic_viscosity = 1.48e-5
            molar_mass = 0.04401
            diff_coefficient = 1.6e-5            
            """)
        end
    end
    return true
end

"""
    get_all_calc_params(filename::String, data_dir::AbstractString = dirname(filename))

Read calculation parameters from a TOML file and return a structured dictionary
with all parsed parameters organized by category.

# Arguments
- `filename::String`: Path to the calculation parameters TOML file
- `data_dir::AbstractString`: Folder any relative CSV path of a time function is
  resolved against, defaulting to the folder holding the calculation file

# Returns
- Dictionary containing all calculation parameters in structured format
"""
function get_all_calc_params(filename::String, data_dir::AbstractString = dirname(filename))
    calc_data = read_calc_params(filename)

    return Dict(
        "units" => get_units(calc_data),
        "gravity" => get_gravity(calc_data),
        "solver_settings" => get_solver_settings(calc_data),
        "time_stepping" => get_time_stepping(calc_data),
        "data_saving_interval" => get_data_saving_interval(calc_data),
        "probing_nodes" => get_probing_nodes(calc_data),
        "probing_elements" => get_probing_elements(calc_data),
        "time_functions" => get_time_functions_section(calc_data, data_dir)
    )
end
