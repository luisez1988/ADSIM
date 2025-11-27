#===============================================================================
# ADSIM Mesh File Writer
# This script writes the geometrical mesh data and boundary/initial conditions
# to a .mesh file for use in the FEM code
#===============================================================================

proc ADSIM::WriteMeshFile { filename } {
    variable current_xml_root
    
    # Initialize the mesh file
    GiD_WriteCalculationFile init $filename
    
    # Get project information
    set project_path [GiD_Info Project ModelName]
    set model_name [file tail $project_path]
    set root [$::gid_groups_conds::doc documentElement] ;# xml document to get some tree data

    set current_xml_root $root
    
    # Write header
    ADSIM::WriteMeshHeader
    
    # Write counters
    ADSIM::WriteMeshCounters
    
    # Write nodal coordinates
    ADSIM::WriteMeshCoordinates
    
    # Write elements
    ADSIM::WriteMeshElements
    
    # Write concentration boundary conditions
    ADSIM::WriteMeshConcentrationBC $root
    
    # Write uniform flow boundary conditions
    ADSIM::WriteMeshFlowBC $root
    
    # Write absolute pressure boundary conditions
    ADSIM::WriteMeshPressureBC $root
    
    # Write initial gas concentrations
    ADSIM::WriteMeshInitialConcentrations $root
    
    # Write initial temperature
    ADSIM::WriteMeshInitialTemperature $root
    
    # Write material assignation
    ADSIM::WriteMeshMaterials $root
    
    # Close the file
    GiD_WriteCalculationFile end
}

#===============================================================================
# Write mesh file header
#===============================================================================
proc ADSIM::WriteMeshHeader { } {
    GiD_WriteCalculationFile puts "### ADSIM MESH FILE ###"
    GiD_WriteCalculationFile puts "### Version: ADSIM_2025 ###"
    GiD_WriteCalculationFile puts ""
}

#===============================================================================
# Write mesh counters (number of nodes and elements)
#===============================================================================
proc ADSIM::WriteMeshCounters { } {
    # Get number of nodes
    set numNodes [GiD_Info Mesh NumNodes]
    
    # Get number of elements
    set numElements [GiD_Info Mesh NumElements]
    
    # Write counters
    GiD_WriteCalculationFile puts "MESH $numNodes $numElements"
    GiD_WriteCalculationFile puts ""
}

#===============================================================================
# Write nodal coordinates
#===============================================================================
proc ADSIM::WriteMeshCoordinates { } {
    set Nodes [GiD_Info Mesh nodes -sublist]
    GiD_WriteCalculationFile puts "coordinates"
    
    # Loop through all nodes
    for {set i 0} {$i < [llength $Nodes]} {incr i} {
        set xcoor [lindex $Nodes  $i 1]
        set ycoor [lindex $Nodes  $i 2]
        GiD_WriteCalculationFile puts "$xcoor $ycoor"
    }
    
    GiD_WriteCalculationFile puts "end coordinates"
    GiD_WriteCalculationFile puts ""
}


#===============================================================================
# Write elements (element connectivity)
#===============================================================================
proc ADSIM::WriteMeshElements { } {
    GiD_WriteCalculationFile puts "elements"
    set ElementList [GiD_Info Mesh Elements Quadrilateral -sublist]
    
    # Loop through all elements (quadrilaterals)
    for {set i 0} {$i < [llength $ElementList]} {incr i} {
        set node1 [lindex $ElementList  $i 1]
        set node2 [lindex $ElementList  $i 2]
        set node3 [lindex $ElementList  $i 3]
        set node4 [lindex $ElementList  $i 4]
        GiD_WriteCalculationFile puts "$node1 $node2 $node3 $node4"
    }

    GiD_WriteCalculationFile puts "end elements"
    GiD_WriteCalculationFile puts ""
}

#===============================================================================
# Write concentration boundary conditions
#===============================================================================
proc ADSIM::WriteMeshConcentrationBC { root } {
    GiD_WriteCalculationFile puts "concentration_bc"
     

    # Get number of gases defined
    set xp_gases {//container[@n="materials"]/container[@n="m_gas"]/blockdata}
    set gas_blocks [$root selectNodes $xp_gases]
    set num_gases [llength $gas_blocks]

    # Line conditions
    set ov_type "line"
    set xp [format_xpath {container[@n="BC"]/condition[@n="gas_fixities"]/group[@ov=%s]} $ov_type]
    set formats ""

    foreach gNode [$root selectNodes $xp] {
    set aux ""
    for {set i 1} {$i <= $num_gases} {incr i} {
        # set string as : concentration_gas_1_
        set format_xpath [format {string(value[@n="concentration_gas_%d_"]/@v)} $i]
        set v1 [$gNode selectNodes $format_xpath]
        append aux "$v1 " 
        
    }
    
    #write nodes with values
    dict set formats [$gNode @n] "%d $aux\n"
    }

    #point conditions
    set ov_type "point" 
    set xp [format_xpath {container[@n="BC"]/condition[@n="gas_fixities"]/group[@ov=%s]} $ov_type]

    foreach gNode [$root selectNodes $xp] {
    set aux ""
    for {set i 1} {$i <= $num_gases} {incr i} {
        # set string as : concentration_gas_1_
        set format_xpath [format {string(value[@n="concentration_gas_%d_"]/@v)} $i]
        set v1 [$gNode selectNodes $format_xpath]
        append aux "$v1 "
        
    }
    #write nodes with values
    dict set formats [$gNode @n] "%d $aux\n"
    }

    # Add a counter 
    set counter [GiD_WriteCalculationFile nodes -count $formats]
    GiD_WriteCalculationFile puts $counter
    GiD_WriteCalculationFile nodes $formats
    GiD_WriteCalculationFile puts "end concentration_bc"
    GiD_WriteCalculationFile puts ""
}

#===============================================================================
# Write uniform flow boundary conditions
#===============================================================================
proc ADSIM::WriteMeshFlowBC { root } {
    GiD_WriteCalculationFile puts "uniform_flow_bc"
    
    # Get number of gases defined
    set xp_gases {//container[@n="materials"]/container[@n="m_gas"]/blockdata}
    set gas_blocks [$root selectNodes $xp_gases]
    set num_gases [llength $gas_blocks]

    # Line conditions
    set ov_type "line"
    set xp [format_xpath {container[@n="BC"]/condition[@n="Gas_Flow"]/group[@ov=%s]} $ov_type]
    set formats ""

    foreach gNode [$root selectNodes $xp] {
    set aux ""
    for {set i 1} {$i <= $num_gases} {incr i} {
        # set string as : flow_gas_1
        set format_xpath [format {string(value[@n="flow_gas_%d"]/@v)} $i]
        set v1 [$gNode selectNodes $format_xpath]
        append aux "$v1 "
        
    }
    
    #write nodes with values
    dict set formats [$gNode @n] "%d $aux\n"
    }

    #point conditions
    set ov_type "point"
    set xp [format_xpath {container[@n="BC"]/condition[@n="Gas_Flow"]/group[@ov=%s]} $ov_type]
    foreach gNode [$root selectNodes $xp] {
    set aux ""
    for {set i 1} {$i <= $num_gases} {incr i} {
        # set string as : flow_gas_1
        set format_xpath [format {string(value[@n="flow_gas_%d"]/@v)} $i]
        set v1 [$gNode selectNodes $format_xpath]
        append aux "$v1 "
        
    }
    #write nodes with values
    dict set formats [$gNode @n] "%d $aux\n"
    }

    # Add a counter 
    set counter [GiD_WriteCalculationFile nodes -count $formats]
    GiD_WriteCalculationFile puts $counter
    GiD_WriteCalculationFile nodes $formats
    GiD_WriteCalculationFile puts "end uniform_flow_bc"
    GiD_WriteCalculationFile puts ""
}

#===============================================================================
# Write absolute pressure boundary conditions
#===============================================================================
proc ADSIM::WriteMeshPressureBC { root } {
    GiD_WriteCalculationFile puts "absolute_pressure"

    # Line conditions
    set ov_type "line"
    set xp [format_xpath {container[@n="BC"]/condition[@n="gas_absolute_boundary"]/group[@ov=%s]} $ov_type]
    set formats ""

    foreach gNode [$root selectNodes $xp] {
        set v1 [$gNode selectNodes {string(value[@n="absolute_pressure"]/@v)}]
        # set gas_idx [$gNode selectNodes {string(value[@n="vacating_gas_index"]/@v)}]
        #write nodes with pressure and gas index
        dict set formats [$gNode @n] "%d $v1\n"
    }

    #point conditions
    set ov_type "point"
    set xp [format_xpath {container[@n="BC"]/condition[@n="gas_absolute_boundary"]/group[@ov=%s]} $ov_type]
    foreach gNode [$root selectNodes $xp] {
        set v1 [$gNode selectNodes {string(value[@n="absolute_pressure"]/@v)}]
        # set gas_idx [$gNode selectNodes {string(value[@n="vacating_gas_index"]/@v)}]
        #write nodes with pressure and gas index
        dict set formats [$gNode @n] "%d $v1\n"
    }

    # Add a counter 
    set counter [GiD_WriteCalculationFile nodes -count $formats]
    GiD_WriteCalculationFile puts $counter
    GiD_WriteCalculationFile nodes $formats
    GiD_WriteCalculationFile puts "end absolute_pressure"
    GiD_WriteCalculationFile puts ""
}

#===============================================================================
# Write initial gas concentrations
#===============================================================================
proc ADSIM::WriteMeshInitialConcentrations { root } {
    GiD_WriteCalculationFile puts "initial_concentrations"
    
    # Get number of gases defined
    set xp_gases {//container[@n="materials"]/container[@n="m_gas"]/blockdata}
    set gas_blocks [$root selectNodes $xp_gases]
    set num_gases [llength $gas_blocks]

    # Surface conditions
    set ov_type "surface"
    set xp [format_xpath {container[@n="initial_conditions"]/condition[@n="Initial_Concentrations"]/group[@ov=%s]} $ov_type]
    set formats ""

    foreach gNode [$root selectNodes $xp] {
        set aux ""
        for {set i 1} {$i <= $num_gases} {incr i} {
            # set string as : concentration_gas_1_
            set format_xpath [format {string(value[@n="concentration_gas_%d_"]/@v)} $i]
            set v1 [$gNode selectNodes $format_xpath]
            append aux "$v1 "
        }
        
        #write elements with values
        dict set formats [$gNode @n] "%d $aux\n"
    }

    # Add a counter 
    set counter [GiD_WriteCalculationFile elements -count $formats]
    GiD_WriteCalculationFile puts $counter
    GiD_WriteCalculationFile elements $formats
    GiD_WriteCalculationFile puts "end initial_concentrations"
    GiD_WriteCalculationFile puts ""
}

#===============================================================================
# Write initial temperature
#===============================================================================
proc ADSIM::WriteMeshInitialTemperature { root } {
    GiD_WriteCalculationFile puts "initial_temperature"

    # Surface conditions (applied to elements)
    set ov_type "surface"
    set xp [format_xpath {//container[@n="initial_conditions"]/condition[@n="Initial_Temperature"]/group[@ov=%s]} $ov_type]
    set formats ""

    foreach gNode [$root selectNodes $xp] {
        set v1 [$gNode selectNodes {string(value[@n="Temperature"]/@v)}]
        #write elements with values
        dict set formats [$gNode @n] "%d $v1\n"
    }

    # Add a counter 
    set counter [GiD_WriteCalculationFile elements -count $formats]
    GiD_WriteCalculationFile puts $counter
    GiD_WriteCalculationFile elements $formats
    GiD_WriteCalculationFile puts "end initial_temperature"
    GiD_WriteCalculationFile puts ""
}

#===============================================================================
# Write material assignation
#===============================================================================
proc ADSIM::WriteMeshMaterials { root } {
    GiD_WriteCalculationFile puts "materials"
    
    # Get all material definitions from the problem type
    set xp_materials {//container[@n="materials"]/container[@n="m_soil"]/blockdata}
    set material_blocks [$root selectNodes $xp_materials]
    
    # Create a mapping of material names to indices
    set material_index 1
    set material_map [dict create]
    foreach mat_block $material_blocks {
        set mat_name [$mat_block @name]
        dict set material_map $mat_name $material_index
        incr material_index
    }
    
    # Surface conditions
    set ov_type "surface"
    set xp [format_xpath {//container[@n="m_soil"]/blockdata/condition[@n="soil_basic"]/group[@ov=%s]} $ov_type]
    set formats ""

    foreach gNode [$root selectNodes $xp] {
        # Find the parent blockdata to get the material name
        set parent [$gNode selectNodes {ancestor::blockdata}]
        if {$parent != ""} {
            set mat_name [$parent @name]
            
            if {[dict exists $material_map $mat_name]} {
                set mat_index [dict get $material_map $mat_name]
                #write elements with material index
                dict set formats [$gNode @n] "%d $mat_index\n"
            }
        }
    }

    # Add a counter 
    set counter [GiD_WriteCalculationFile elements -count $formats]
    GiD_WriteCalculationFile puts $counter
    GiD_WriteCalculationFile elements $formats
    GiD_WriteCalculationFile puts "end materials"    
}
