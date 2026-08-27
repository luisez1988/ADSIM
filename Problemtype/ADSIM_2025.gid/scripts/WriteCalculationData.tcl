#===============================================================================
# ADSIM Calculation Data Writer (TOML format)
# This script writes calculation parameters, units, and simulation settings
# to a .toml file for use in the FEM code
#===============================================================================

proc ADSIM::WriteCalculationData { filename } {
    # Initialize the file
    GiD_WriteCalculationFile init $filename
    
    # Get root of XML tree
    set root [$::gid_groups_conds::doc documentElement]
    
    # Write header
    ADSIM::WriteCalculationHeader
    
    # Write units
    ADSIM::WriteUnits $root
    
    # Write gravity vector
    ADSIM::WriteGravity $root
    
    # Write solver settings
    ADSIM::WriteSolver $root
    
    # Write time stepping data
    ADSIM::WriteTimeStepping $root
    
    # Write data saving and probing settings
    ADSIM::WriteDataSaving $root

    # Write time functions (curves that drive boundary conditions in time)
    ADSIM::WriteTimeFunctions $root

    # Close the file
    GiD_WriteCalculationFile end
}

#===============================================================================
# Write time functions
#
# One [[time_function]] table per curve defined in the Time functions container.
# The id is the 1 based position of the curve in the tree, which is exactly what
# ADSIM::GetTimeFunctionMap in WriteMeshFile.tcl assigns when it turns the name
# picked in a boundary condition into the integer written to the mesh file. The
# two must stay in step: they are both derived from document order.
#
# Only the keys that matter for the chosen type are written, so a ramp does not
# emit period and phase and the calculation file stays readable.
#===============================================================================
proc ADSIM::WriteTimeFunctions { root } {
    set xp {//container[@n="time_functions"]/blockdata}
    set blocks [$root selectNodes $xp]

    if {[llength $blocks] == 0} {
        return
    }

    GiD_WriteCalculationFile puts "# Curves f(t) driving the time dependent boundary conditions."
    GiD_WriteCalculationFile puts "# Referenced by id from the *_tf blocks of the mesh file."
    GiD_WriteCalculationFile puts ""

    set id 0
    foreach block $blocks {
        incr id

        set name [$block @name]
        set type [$block selectNodes {string(value[@n="type"]/@v)}]
        set mode [$block selectNodes {string(value[@n="mode"]/@v)}]

        if {$type == ""} { set type "constant" }
        if {$mode == ""} { set mode "absolute" }

        GiD_WriteCalculationFile puts "\[\[time_function\]\]"
        GiD_WriteCalculationFile puts "id = $id"
        GiD_WriteCalculationFile puts "name = \"$name\""
        GiD_WriteCalculationFile puts "type = \"$type\""
        GiD_WriteCalculationFile puts "mode = \"$mode\""

        switch -- $type {
            "ramp" {
                foreach key {t_start t_end v_start v_end} {
                    set v [$block selectNodes [format {string(value[@n="%s"]/@v)} $key]]
                    if {$v == ""} { set v 0.0 }
                    GiD_WriteCalculationFile puts "$key = $v"
                }
            }
            "periodic" {
                set wave [$block selectNodes {string(value[@n="wave"]/@v)}]
                if {$wave == ""} { set wave "sine" }
                GiD_WriteCalculationFile puts "wave = \"$wave\""
                foreach key {mean amplitude period phase t_start} {
                    set v [$block selectNodes [format {string(value[@n="%s"]/@v)} $key]]
                    if {$v == ""} { set v 0.0 }
                    GiD_WriteCalculationFile puts "$key = $v"
                }
            }
            "table" {
                # Only the file name is written: the CSV itself is copied into the
                # generated folder by ADSIM::CopyTimeSeriesFiles, and the solver
                # resolves a relative name against the folder holding this file.
                set file [$block selectNodes {string(value[@n="file"]/@v)}]
                GiD_WriteCalculationFile puts "file = \"[file tail $file]\""

                foreach key {interpolation before after} {
                    set v [$block selectNodes [format {string(value[@n="%s"]/@v)} $key]]
                    if {$v == ""} { set v "hold" }
                    GiD_WriteCalculationFile puts "$key = \"$v\""
                }
                foreach key {time_offset time_scale} {
                    set v [$block selectNodes [format {string(value[@n="%s"]/@v)} $key]]
                    if {$v == ""} { set v [expr {$key eq "time_scale" ? 1.0 : 0.0}] }
                    GiD_WriteCalculationFile puts "$key = $v"
                }

                set avg [$block selectNodes {string(value[@n="average"]/@v)}]
                GiD_WriteCalculationFile puts "average = [expr {$avg == 1 ? "true" : "false"}]"
            }
            default {
                set v [$block selectNodes {string(value[@n="value"]/@v)}]
                if {$v == ""} { set v 0.0 }
                GiD_WriteCalculationFile puts "value = $v"
            }
        }

        GiD_WriteCalculationFile puts ""
    }
}

#===============================================================================
# Write calculation file header
#===============================================================================
proc ADSIM::WriteCalculationHeader { } {
    GiD_WriteCalculationFile puts "# ADSIM calculation parameters file header (need to add disclaimer and license)"
    GiD_WriteCalculationFile puts "# ADSIM_version = \"2025 v0.1.0\""
    GiD_WriteCalculationFile puts "# File_created = \"[clock format [clock seconds] -format %Y-%m-%d]\""
    GiD_WriteCalculationFile puts ""
}

#===============================================================================
# Write units section
#===============================================================================
proc ADSIM::WriteUnits { root } {
    # Get units container
    set xp_units {//container[@n="Units_Dimensions"]/container[@n="units"]}
    set units_container [$root selectNodes $xp_units]
    
    if {$units_container == ""} {
        return
    }
    
    GiD_WriteCalculationFile puts "\[units\]"
    
    # Geometry unit
    set geom_unit [$units_container selectNodes {string(value[@n="units_mesh"]/@v)}]
    GiD_WriteCalculationFile puts "geometry_unit = \"$geom_unit\""
    
    # Mass unit
    set mass_unit [$units_container selectNodes {string(value[@n="units_mass"]/@v)}]
    GiD_WriteCalculationFile puts "mass_unit = \"$mass_unit\""
    
    # Pressure unit
    set pressure_unit [$units_container selectNodes {string(value[@n="units_pressure"]/@v)}]
    GiD_WriteCalculationFile puts "pressure_unit = \"$pressure_unit\""
    
    # Temperature unit
    set temp_unit [$units_container selectNodes {string(value[@n="units_temperature"]/@v)}]
    GiD_WriteCalculationFile puts "temperature_unit = \"$temp_unit\""
    
    # Time unit
    set time_unit [$units_container selectNodes {string(value[@n="units_time"]/@v)}]
    GiD_WriteCalculationFile puts "time_unit = \"$time_unit\""
    
    GiD_WriteCalculationFile puts ""
}

#===============================================================================
# Write gravity vector
#===============================================================================
proc ADSIM::WriteGravity { root } {
    # Get gravity container
    set xp_gravity {//container[@n="Units_Dimensions"]/container[@n="gravity"]}
    set gravity_container [$root selectNodes $xp_gravity]
    
    if {$gravity_container == ""} {
        return
    }
    
    GiD_WriteCalculationFile puts "\[gravity\]"
    
    # Gravity magnitude
    set grav_mag [$gravity_container selectNodes {string(value[@n="gravity_acceleration"]/@v)}]
    GiD_WriteCalculationFile puts "gravity_magnitude = $grav_mag"
    
    # Gravity x component
    set grav_x [$gravity_container selectNodes {string(value[@n="gravity_x"]/@v)}]
    GiD_WriteCalculationFile puts "gravity_x_component = $grav_x"
    
    # Gravity y component
    set grav_y [$gravity_container selectNodes {string(value[@n="gravity_y"]/@v)}]
    GiD_WriteCalculationFile puts "gravity_y_component = $grav_y"
    
    GiD_WriteCalculationFile puts ""
}

#===============================================================================
# Write solver settings
#===============================================================================
proc ADSIM::WriteSolver { root } {
    # Get solver container
    set xp_solver {//container[@n="Calculation_Data"]/container[@n="solver_type"]}
    set solver_container [$root selectNodes $xp_solver]
    
    if {$solver_container == ""} {
        return
    }
    
    GiD_WriteCalculationFile puts "\[solver\]"
    
    # Dimension selection
    set dimension [$solver_container selectNodes {string(value[@n="dimension_selection"]/@v)}]
    GiD_WriteCalculationFile puts "solver_type = \"$dimension\""
    
    # Calculation mode components
    set calc_mode_container [$solver_container selectNodes {container[@n="calculation_mode"]}]
    
    if {$calc_mode_container != ""} {
        # Diffusion
        set diffusion [$calc_mode_container selectNodes {string(value[@n="diffusion"]/@v)}]
        GiD_WriteCalculationFile puts "diffusion = $diffusion"
        
        # Advection
        set advection [$calc_mode_container selectNodes {string(value[@n="advection"]/@v)}]
        GiD_WriteCalculationFile puts "advection = $advection"
        
        # Gravity
        set gravity [$calc_mode_container selectNodes {string(value[@n="gravity"]/@v)}]
        GiD_WriteCalculationFile puts "gravity = $gravity"
        
        # Reaction kinetics
        set reaction_kinetics [$calc_mode_container selectNodes {string(value[@n="reaction_kinetics"]/@v)}]
        GiD_WriteCalculationFile puts "reaction_kinetics = $reaction_kinetics"

        # Heat transport. With both off the temperature is held fixed and the
        # reaction enthalpy is ignored.
        set heat_conduction [$calc_mode_container selectNodes {string(value[@n="heat_conduction"]/@v)}]
        if {$heat_conduction == ""} { set heat_conduction 0 }
        GiD_WriteCalculationFile puts "heat_conduction = $heat_conduction"

        set heat_advection [$calc_mode_container selectNodes {string(value[@n="heat_advection"]/@v)}]
        if {$heat_advection == ""} { set heat_advection 0 }
        GiD_WriteCalculationFile puts "heat_advection = $heat_advection"
    }
    
    GiD_WriteCalculationFile puts ""
}

#===============================================================================
# Write time stepping data
#===============================================================================
proc ADSIM::WriteTimeStepping { root } {
    # Get time stepping container
    set xp_time {//container[@n="Calculation_Data"]/container[@n="calculation_stepping_data"]}
    set time_container [$root selectNodes $xp_time]
    
    if {$time_container == ""} {
        return
    }
    
    GiD_WriteCalculationFile puts "\[time_stepping\]"
    
    # Total simulation time
    set sim_time [$time_container selectNodes {string(value[@n="simulation_time"]/@v)}]
    GiD_WriteCalculationFile puts "total_simulation_time = $sim_time"
    
    # Time per step
    set time_step [$time_container selectNodes {string(value[@n="time_step"]/@v)}]
    GiD_WriteCalculationFile puts "time_per_step = $time_step"
    
    # Courant number
    set courant [$time_container selectNodes {string(value[@n="courant_number"]/@v)}]
    GiD_WriteCalculationFile puts "courant_number = $courant"
    
    GiD_WriteCalculationFile puts ""
}

#===============================================================================
# Write data saving and probing settings
#===============================================================================
proc ADSIM::WriteDataSaving { root } {
    # Get data saving container
    set xp_saving {//container[@n="Calculation_Data"]/container[@n="data_saving"]}
    set saving_container [$root selectNodes $xp_saving]
    
    if {$saving_container == ""} {
        return
    }
    
    GiD_WriteCalculationFile puts "\[data_saving\]"
    
    # Data saving interval
    set save_interval [$saving_container selectNodes {string(value[@n="data_saving_interval"]/@v)}]
    GiD_WriteCalculationFile puts "data_saving_interval = $save_interval"
    
    GiD_WriteCalculationFile puts ""
    GiD_WriteCalculationFile puts "\[probing\]"
    
    # Nodes to follow
    set num_nodes [$saving_container selectNodes {string(.//container[@n="nodes_to_follow"]/value[@n="number_nodes"]/@v)}]
    GiD_WriteCalculationFile puts "# Nodes to follow during simulation"
    GiD_WriteCalculationFile puts "number_of_nodes = $num_nodes"
    
    # Build nodes array
    set node_list [list]
    for {set i 1} {$i <= $num_nodes} {incr i} {
        set node_xpath [format {string(.//container[@n="nodes_to_follow"]/value[@n="node_%d"]/@v)} $i]
        set node_id [$saving_container selectNodes $node_xpath]
        if {$node_id != "" && $node_id > 0} {
            lappend node_list $node_id
        }
    }
    if {[llength $node_list] > 0} {
        GiD_WriteCalculationFile puts "nodes_to_probe = \[[join $node_list ", "]\]"
    } else {
        GiD_WriteCalculationFile puts "nodes_to_probe = \[\]"
    }
    
    GiD_WriteCalculationFile puts ""
    
    # Elements to follow
    set num_elements [$saving_container selectNodes {string(.//container[@n="elements_to_follow"]/value[@n="number_elements"]/@v)}]
    GiD_WriteCalculationFile puts "# Elements to follow during simulation (data at Gaussian points)"
    GiD_WriteCalculationFile puts "number_of_elements = $num_elements"
    
    # Build elements array
    set elem_list [list]
    for {set i 1} {$i <= $num_elements} {incr i} {
        set elem_xpath [format {string(.//container[@n="elements_to_follow"]/value[@n="element_%d"]/@v)} $i]
        set elem_id [$saving_container selectNodes $elem_xpath]
        if {$elem_id != "" && $elem_id > 0} {
            lappend elem_list $elem_id
        }
    }
    if {[llength $elem_list] > 0} {
        GiD_WriteCalculationFile puts "elements_to_probe = \[[join $elem_list ", "]\]"
    } else {
        GiD_WriteCalculationFile puts "elements_to_probe = \[\]"
    }
    
    GiD_WriteCalculationFile puts ""
}
