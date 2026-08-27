proc InitGIDProject { dir } {

    global _dir
    set _dir $dir
    set adsim_version "ADSIM v0.1.0"
    set GiDVersionRequired "14.0"

    ADSIM::SetDir $dir ;#store to use it later
    ADSIM::LoadScripts
    GidUtils::OpenWindow CUSTOMLIB

	# disclaimer at start-up
	set self_close 0
    GidUtils::Splash [file join $_dir images Splash_screen_ADSIM.png] .splash $self_close [list $adsim_version 445 10]

	# check minimal GiD version
	set message_required "This $adsim_version interface is developed for GiD $GiDVersionRequired or later. \n \n It is advised to update your GiD software."
	set title_required "ADSIM - Required GiD Version"
    if { [GidUtils::VersionCmp $GiDVersionRequired] < 0 } { tk_messageBox -title $title_required -message $message_required -icon warning -type ok }

	# add ADSIM menu
	GiDMenu::Create "ADSIM" "PRE" 5 =
    GiDMenu::InsertOption "ADSIM" [list "Generate ADSIM Files"] 0 PRE "ADSIM::Calculate" "Control-B" "[file join $dir images generate.png]" replace =
    GiDMenu::InsertOption "ADSIM" [list "Import time series CSV..."] 1 PRE "ADSIM::ImportTimeSeries" "" "[file join $dir images icon_scientific.png]" replace =
    GiDMenu::InsertOption "ADSIM" [list "---"] 2 PRE "" "" "" replace =
	GiDMenu::InsertOption "ADSIM" [list "GitHub Repository..."] 3 PRE "ADSIM::Github" "" "[file join $dir images icon_Github.png]" replace =
	GiDMenu::InsertOption "ADSIM" [list "Scientific Manuscripts..."] 4 PRE "ADSIM::Scientific" "" "[file join $dir images icon_scientific.png]" replace =
	GiDMenu::InsertOption "ADSIM" [list "Sample projects..."] 5 PRE "ADSIM::Sample" "" "[file join $dir images icon_samples.png]" replace =
    GiDMenu::InsertOption "ADSIM" [list "Disclaimer..."] 6 PRE "ADSIM::Disclaimer" "" "[file join $dir images icon_ADSIM.png]" replace =
	GiDMenu::InsertOption "ADSIM" [list "About..."] 7 PRE "ADSIM::About" "" "[file join $dir images icon_ADSIM.png]" replace =

	# simplify and adaptd HELP menu
	#GiDMenu::RemoveOption "Help" [list "Customization help"] "PRE" _
	GiDMenu::RemoveOption "Help" [list "Tutorials"] "PRE" _
    GiDMenu::RemoveOption "Help" [list "What is new"] "PRE" _
    GiDMenu::RemoveOption "Help" [list "FAQ"] "PRE" _
    GiDMenu::RemoveOption "Help" [list "Register problem type"] "PRE" _
    GiDMenu::RemoveOption "Help" [list "Register from file"] "PRE" _

    # remove options in DATA menu
	GidChangeDataLabel "Interval" ""
	GidChangeDataLabel "Local axes" ""
	GidChangeDataLabel "Materials" ""
	GidChangeDataLabel "Conditions" ""
	GidChangeDataLabel "Problem data" ""

	# remove CALCULATE menu
	GiDMenu::Delete "Calculate" PRE

    # apply menu update
    GiDMenu::UpdateMenus

    # open window tree
    GidUtils::OpenWindow CUSTOMLIB

}

proc EndGIDProject { } {

}

namespace eval ADSIM {
}

namespace eval ADSIM_2025 {
    variable problemtype_dir
}

proc ADSIM::SetDir { dir } {
    variable problemtype_dir
    set problemtype_dir $dir
}

proc ADSIM::GetDir { } {
  variable problemtype_dir
    return $problemtype_dir
}

proc ADSIM::Github { } {

    global _dir
    set url "https://github.com/EspinosaAndres/GRAPE"
    eval exec [auto_execok start] \"\" [list $url]

}

proc ADSIM::Scientific { } {

    global _dir
    tk_messageBox -title "ADSIM - Scientific" -message "in construction" -icon info -type ok

}

proc ADSIM::Sample { } {

    global _dir
    tk_messageBox -title "ADSIM - Sample" -message "in construction" -icon info -type ok
    # set TestDoc [file join $_dir doc "SampleManual_2022.pdf"]
    #eval exec [auto_execok start] \"\" [list $TestDoc]

}

proc ADSIM::About { } {

    tk_messageBox -title "ADSIM - About" -message "ADSIM version 0.x.x: Project in construction.
    Stay tunned!" -icon info

}

proc ADSIM::Disclaimer { } {

set answer [tk_messageBox -title "ADSIM - Disclaimer" -message "Copyright (C) 2020  Luis Zambrano-Cruzatty, Paula Sarmiento, Andrés Espinosa

This program is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>

ADSIM Development Team
E-mail:	lezambra@ncsu.edu
Web: 	go.ncsu.edu/zgem" -icon info -type ok]

}

#proc ADSIM::CreateBatch { channel } {

#    set project_path [GiD_Info Project ModelName]
#    set model_name [file tail $project_path]
#    set exe_name [GiD_Info Project ProblemType]
#    GiD_File fprintf -nonewline $channel "\""
#    GiD_File fprintf -nonewline $channel $project_path
#    GiD_File fprintf -nonewline $channel ".GP\\"
#    GiD_File fprintf -nonewline $channel $exe_name
#    GiD_File fprintf -nonewline $channel ".exe\" "
#    GiD_File fprintf -nonewline $channel "\""
#    GiD_File fprintf -nonewline $channel $project_path
#    GiD_File fprintf -nonewline $channel ".GP\\"
#    GiD_File fprintf -nonewline $channel $model_name
#    GiD_File fprintf $channel "\""
#    GiD_File fprintf $channel "PAUSE"

#}

proc ADSIM::Calculate { } {
	set project_path [GiD_Info Project ModelName]
    set model_name [file tail $project_path]
	
	# Check if project has been saved
	if {$project_path == "" || $model_name == ""} {
		tk_messageBox -title "ADSIM - Error" -message "Please save your project before generating ADSIM files!" -icon error -type ok
		return
	}
	
	# Create output folder with .ADSIM extension
	set output_folder "${project_path}.ADSIM"
	
    # Recreate output folder: delete if it exists, then create fresh
    if {[file exists $output_folder]} {
        file delete -force $output_folder
    }
    file mkdir $output_folder
	
	# Set file paths to the .ADSIM folder
	set calc_params_file [file join $output_folder "${model_name}_calc.toml"]
	set material_file [file join $output_folder "${model_name}_mat.toml"]
	set mesh_file [file join $output_folder "${model_name}.mesh"]

	
    set answer [tk_messageBox -title "ADSIM - Generate ADSIM Files" -message "Before generating the ADSIM input files, please confirm:\n\n - All materials are assigned\n - Boundary conditions are defined\n - Calculation parameters are set\n - The mesh has been (re)generated.\n\nSave the project, and ensure the project name contains no spaces!" -icon warning -type okcancel]

    switch -- $answer {
      cancel return
      ok {
	  #GiD_Process Mescape Utilities Calculate	  
    #   Generates the calculation parameters, material data, and mesh file
	  ADSIM::WriteCalculationData $calc_params_file
	  ADSIM::WriteMaterialData $material_file
	  ADSIM::WriteMeshFile $mesh_file

	  # Copy any CSV a time function reads. This has to happen after the writers,
	  # because the output folder is deleted and recreated above.
	  set missing [ADSIM::CopyTimeSeriesFiles $output_folder]

	  # Show success message
	  if {[llength $missing] > 0} {
	      tk_messageBox -title "ADSIM - Missing time series files" -message "The input files were generated in:\n\n$output_folder\n\nBut these time series files could not be found and were NOT copied:\n\n[join $missing \n]\n\nThe calculation will fail until they are placed next to the project. Use the ADSIM menu, Import time series CSV, to bring them in." -icon warning -type ok
	  } else {
	      tk_messageBox -title "ADSIM - Files Generated" -message "ADSIM input files have been successfully generated in:\n\n$output_folder" -icon info -type ok
	  }
	  }
    }
}


#===============================================================================
# Copy the CSV files read by table time functions into the generated folder.
#
# The calculation file names each table only by its basename, and the solver
# resolves that against the folder holding the calculation file, so the data has
# to sit beside it. The generated folder is wiped on every run, which is why this
# is a copy at generation time rather than something the user does once.
#
# A name is looked for beside the project first, then as given, so both a bare
# file name and a full path typed into the tree work.
#
# Returns the list of files that could not be found, for the caller to report.
#===============================================================================
proc ADSIM::CopyTimeSeriesFiles { output_folder } {
    set missing [list]

    set root [$::gid_groups_conds::doc documentElement]
    set project_dir [file dirname [GiD_Info Project ModelName]]

    foreach block [$root selectNodes {//container[@n="time_functions"]/blockdata}] {
        set type [$block selectNodes {string(value[@n="type"]/@v)}]
        if {$type ne "table"} {
            continue
        }

        set file [$block selectNodes {string(value[@n="file"]/@v)}]
        if {$file eq ""} {
            lappend missing "[$block @name] (no file given)"
            continue
        }

        set candidates [list [file join $project_dir $file] $file]
        set source ""
        foreach candidate $candidates {
            if {[file exists $candidate] && [file isfile $candidate]} {
                set source $candidate
                break
            }
        }

        if {$source eq ""} {
            lappend missing "[$block @name]: $file"
            continue
        }

        set target [file join $output_folder [file tail $file]]
        if {[catch {file copy -force $source $target} err]} {
            lappend missing "[$block @name]: $file ($err)"
        }
    }

    return $missing
}

proc ADSIM::LoadScripts { } {
    variable problemtype_dir
    # Common scripts
    #set script_files [list calculate.tcl createGOM.tcl createCPS.tcl createOPD.tcl]
    set script_files [list WriteMaterialData.tcl WriteCalculationData.tcl WriteMeshFile.tcl]
    foreach filename $script_files {
        uplevel #0 [list source [file join $problemtype_dir scripts $filename]]
    }
}



### auxiliary procedures to be called from the .bas templates to write some of its parts

proc ADSIM::WriteDamping { condition_name } {
    set result ""
    foreach item [GiD_Info conditions $condition_name mesh] {
        set data_by_element([lindex $item 1]) [lindex $item 3]
    }
    foreach element_id [GiD_Mesh list element] {
        if { [info exists data_by_element($element_id)] } {
            append result "$data_by_element($element_id) \n"
        } else {
            append result "0 \n"
        }
    }
    return $result
}

proc ADSIM::WriteNumberOfMaterialPointsSP { condition_name } {
    set result ""
    foreach item [GiD_Info conditions $condition_name mesh] {
        set data_by_element([lindex $item 1]) [lindex $item 3]
    }
    foreach element_id [GiD_Mesh list element] {
        if { [info exists data_by_element($element_id)] } {
            append result "$data_by_element($element_id) 0 \n"
        } else {
            append result "0 0 \n"
        }
    }
    return $result
}

proc ADSIM::WriteNumberOfMaterialPointsDP { condition_name } {
    set result ""
    foreach item [GiD_Info conditions $condition_name mesh] {
        set data1_by_element([lindex $item 1]) [lindex $item 3]
		set data2_by_element([lindex $item 1]) [lindex $item 4]
    }
    foreach element_id [GiD_Mesh list element] {
        if { [info exists data1_by_element($element_id)] || [info exists data2_by_element($element_id)] } {
            append result "$data1_by_element($element_id) $data2_by_element($element_id) \n"
        } else {
            append result "0 0 \n"
        }
    }
    return $result
}

#to be used by TKWIDGET to adda layer thickness
#e.g.
#QUESTION: your_question
#VALUE: soil_thickness
#TKWIDGET: GidUtils::TkwidgetPickSoilLayer

proc GidUtils::TkwidgetPickSoilLayer { event args } {
    global tkwidgedprivpicknodebuttons
    switch $event {
        INIT {
            lassign $args PARENT current_row_variable GDN STRUCT QUESTION
            upvar $current_row_variable ROW
            set entry ""
            set entry_gridded 0
            foreach item [grid slaves $PARENT -row [expr $ROW-1]] {
                if { [winfo class $item] == "Entry"  || [winfo class $item] == "TEntry" } {
                    #assumed that it is the only entry of this row
                    set entry $item
                    set entry_gridded 1
                    break
                }
            }
            if { $entry == "" } {
                set entry [GidUtils::DarkTrickTryFindUngriddedEntry $PARENT $QUESTION]
            }
            if { $entry != "" } {
                set tkwidgedprivpicknodebuttons($QUESTION) [ttk::button $PARENT.bpicknode$QUESTION \
                        -image [gid_themes::GetImage "dimension_dist.png" small_icons] \
                        -command [list GetThickness $entry]]
                grid $tkwidgedprivpicknodebuttons($QUESTION) -row [expr $ROW-1] -column 2 -sticky w
                grid configure $entry -sticky ew
                if { !$entry_gridded } {
                    grid remove $entry
                    grid remove $tkwidgedprivpicknodebuttons($QUESTION)
                }
            }
            return ""
        }
        SYNC {
            #lassign $args GDN STRUCT QUESTION
            #DWLocalSetValue $GDN $STRUCT $QUESTION $value
        }
        DEPEND {
            lassign $args GDN STRUCT QUESTION ACTION VALUE
            if { [info exists tkwidgedprivpicknodebuttons($QUESTION)] && \
                     [winfo exists $tkwidgedprivpicknodebuttons($QUESTION)] } {
                if { $ACTION == "HIDE" } {
                    grid remove $tkwidgedprivpicknodebuttons($QUESTION)
                } else {
                    #RESTORE
                    grid $tkwidgedprivpicknodebuttons($QUESTION)
                }
            }
        }
        CLOSE {
            array unset tkwidgedprivpicknodebuttons
        }
        default {
            return [list ERROR [_ "Unexpected tkwidget event"]]
        }
    }
    #a tkwidget procedure must return "" if Ok or [list ERROR $description] or [list WARNING $description]
    return ""
}

#Required for the previous function
proc GetThickness {entry} {
    set p1 [GidUtils::GetCoordinates [_ "Enter first point (ESC to leave)"]]
    if { $p1=="" } return ""
    set p2 [GidUtils::GetCoordinates [_ "Enter second point (ESC to leave)"]]
    if { $p2=="" } return ""

    set d [MathUtils::VectorDistance $p1 $p2]

    $entry delete 0 end
    $entry insert end $d

}


#to be used by TKWIDGET to pick soil surface
#e.g.
#QUESTION: your_question
#VALUE: soil_surface (Y coordinate)
#TKWIDGET: GidUtils::TkwidgetSoilElevation

proc GidUtils::TkwidgetSoilElevation { event args } {
    global tkwidgedprivpicknodebuttons
    switch $event {
        INIT {
            lassign $args PARENT current_row_variable GDN STRUCT QUESTION
            upvar $current_row_variable ROW
            set entry ""
            set entry_gridded 0
            foreach item [grid slaves $PARENT -row [expr $ROW-1]] {
                if { [winfo class $item] == "Entry"  || [winfo class $item] == "TEntry" } {
                    #assumed that it is the only entry of this row
                    set entry $item
                    set entry_gridded 1
                    break
                }
            }
            if { $entry == "" } {
                set entry [GidUtils::DarkTrickTryFindUngriddedEntry $PARENT $QUESTION]
            }
            if { $entry != "" } {
                set tkwidgedprivpicknodebuttons($QUESTION) [ttk::button $PARENT.bpicknode$QUESTION \
                        -image [gid_themes::GetImage "point.png" small_icons] \
                        -command [list GetElevation $entry]]
                grid $tkwidgedprivpicknodebuttons($QUESTION) -row [expr $ROW-1] -column 2 -sticky w
                grid configure $entry -sticky ew
                if { !$entry_gridded } {
                    grid remove $entry
                    grid remove $tkwidgedprivpicknodebuttons($QUESTION)
                }
            }
            return ""
        }
        SYNC {
            #lassign $args GDN STRUCT QUESTION
            #DWLocalSetValue $GDN $STRUCT $QUESTION $value
        }
        DEPEND {
            lassign $args GDN STRUCT QUESTION ACTION VALUE
            if { [info exists tkwidgedprivpicknodebuttons($QUESTION)] && \
                     [winfo exists $tkwidgedprivpicknodebuttons($QUESTION)] } {
                if { $ACTION == "HIDE" } {
                    grid remove $tkwidgedprivpicknodebuttons($QUESTION)
                } else {
                    #RESTORE
                    grid $tkwidgedprivpicknodebuttons($QUESTION)
                }
            }
        }
        CLOSE {
            array unset tkwidgedprivpicknodebuttons
        }
        default {
            return [list ERROR [_ "Unexpected tkwidget event"]]
        }
    }
    #a tkwidget procedure must return "" if Ok or [list ERROR $description] or [list WARNING $description]
    return ""
}

#Required for the previous function
proc GetElevation {entry} {
    set p1 [GidUtils::GetCoordinates [_ "Enter first point (ESC to leave)"]]
    if { $p1=="" } return ""
    set Elev [lindex $p1 1]
    $entry delete 0 end
    $entry insert end $Elev
}

#Required for the previous function
proc PickNodeId {entry} {
    set nodeId [GiD_Info incolor pick nodes]
    if { $nodeId == "" || $nodeId <= 0 } {
        return ""
    }
    $entry delete 0 end
    $entry insert end $nodeId
}


######################################################################
#  auxiliary procs invoked from the tree (see .spd xml description)
proc ADSIM_2025::GetMaterialsList { domNode } {
    set x_path {//container[@n="m_soils"]}
    set dom_materials [$domNode selectNodes $x_path]
    if { $dom_materials == "" } {
        error [= "xpath '%s' not found in the spd file" $x_path]
    }
    set result [list]
    foreach dom_material [$dom_materials childNodes] {
        set name [$dom_material @name]
        lappend result $name
    }
    return [join $result ,]
}

######################################################################
#  auxiliary procs invoked from the tree (see .spd xml description)
proc ADSIM_2025::GetGasDistributionsList { domNode } {
    set x_path {//container[@n="m_concentrations"]}
    set dom_materials [$domNode selectNodes $x_path]
    if { $dom_materials == "" } {
        error [= "xpath '%s' not found in the spd file" $x_path]
    }
    set result [list]
    foreach dom_material [$dom_materials childNodes] {
        set name [$dom_material @name]
        lappend result $name
    }
    return [join $result ,]
}


proc check_dim_points {dim1 dim2 npoints domNode} {
      set dim_path {string(//container[@n="Units_Dimensions"]/value[@n="NDIM"]/@v)}
	  set point_path {string(//container[@n="Units_Dimensions"]/value[@n="NLAYERS"]/@v)}
	  set problem_type [$domNode selectNodes $dim_path]
	  set num_points [$domNode selectNodes $point_path]
	  if {[string equal $problem_type $dim1]} {
	  if {[string equal $num_points $npoints]} {
	  return normal}
	  }
	  if {[string equal $problem_type $dim2]} {
	  if {[string equal $num_points $npoints]} {return normal}
	  }
      return hidden
}

proc find_material_id {material_name root} {
    # Search for soils
    set xp_soil {//container[@n="materials"]/container[@n="m_soil"]/blockdata}
    set soil_list [$root selectNodes $xp_soil]
    set int 1
    foreach gNode $soil_list {	
	    set type [$gNode getAttribute "name"]	
	    if {$type eq $material_name} {
	      return $int
	    }	    
	    set int [expr $int + 1]    
    }
    
    # If no matching blockdata element is found, return an empty string.
    return ""
}

######################################################################
# Check if a gas with the given index exists in Materials
# Returns "normal" if gas exists, "hidden" if it doesn't
#===============================================================================
# Time functions
#===============================================================================

#-------------------------------------------------------------------------------
# Fill the "Time variation" combo of a boundary condition.
#
# Returns "Constant" followed by the name of every curve defined in the Time
# functions container, in tree order. "Constant" is always first and is the
# default, so a boundary left alone behaves exactly as it did before time
# functions existed.
#-------------------------------------------------------------------------------
proc ADSIM_2025::GetTimeFunctionsList { } {
    set result [list "Constant"]

    set root [customlib::GetBaseRoot]
    if {$root == ""} {
        return [join $result ,]
    }

    set xp {//container[@n="time_functions"]/blockdata}
    foreach node [$root selectNodes $xp] {
        set name [$node @name]
        if {$name ne ""} {
            lappend result $name
        }
    }

    return [join $result ,]
}

#-------------------------------------------------------------------------------
# Show or hide a field of a time function blockdata according to its type.
#
# `types` is a comma separated list of the types the field belongs to, so
# state="[ADSIM_2025::TimeFunctionFieldState ramp,periodic %W]" shows the field
# for a ramp and for a periodic curve and hides it otherwise. Same contract as
# CheckGasExists above: return the literal "normal" or "hidden".
#
# %W is substituted by GiD with the node the state is being evaluated for, so
# the xpath is relative to it and finds the type of *this* curve. domNode comes
# last, the convention every other state callback in the GiD tree follows.
#-------------------------------------------------------------------------------
proc ADSIM_2025::TimeFunctionFieldState { types domNode } {
    if {$domNode == ""} {
        return "normal"
    }

    # Sibling value inside the same blockdata
    set current [$domNode selectNodes {string(../value[@n="type"]/@v)}]
    if {$current eq ""} {
        return "normal"
    }

    if {[lsearch -exact [split $types ,] $current] >= 0} {
        return "normal"
    }
    return "hidden"
}

#-------------------------------------------------------------------------------
# Browse for a time series CSV and copy it next to the project.
#
# The generated input folder is deleted and rebuilt on every run, so the file
# has to live with the project and be copied in at generation time; that copy is
# done by ADSIM::CopyTimeSeriesFiles. This command only puts the file somewhere
# stable and tells the user the name to type into the CSV file field.
#-------------------------------------------------------------------------------
proc ADSIM::ImportTimeSeries { } {
    set project_path [GiD_Info Project ModelName]
    if {$project_path == "" || $project_path == "UNNAMED"} {
        tk_messageBox -title "ADSIM - Import time series" \
            -message "Please save your project before importing a time series file." \
            -icon error -type ok
        return
    }

    set source [tk_getOpenFile -title "Select a time series CSV file" \
                    -filetypes {{"CSV files" {.csv}} {"Text files" {.txt}} {"All files" {*}}}]
    if {$source == ""} {
        return
    }

    set target_dir [file dirname $project_path]
    set basename [file tail $source]
    set target [file join $target_dir $basename]

    if {[file normalize $source] ne [file normalize $target]} {
        if {[catch {file copy -force $source $target} err]} {
            tk_messageBox -title "ADSIM - Import time series" \
                -message "Could not copy the file:\n\n$err" -icon error -type ok
            return
        }
    }

    tk_messageBox -title "ADSIM - Import time series" \
        -message "Copied to:\n\n$target\n\nNow create a Time function of type 'table' and enter this name in its CSV file field:\n\n$basename" \
        -icon info -type ok
}

proc ADSIM_2025::CheckGasExists { gas_index } {
    set root [customlib::GetBaseRoot]
    if {$root == ""} {
        # If we can't get root, show first gas by default
        if {$gas_index == 1} {
            return "normal"
        } else {
            return "hidden"
        }
    }
    
    set xp {//container[@n="materials"]/container[@n="m_gas"]/blockdata}
    set gas_blocks [$root selectNodes $xp]
    set num_gases [llength $gas_blocks]
    
    if {$gas_index <= $num_gases} {
        return "normal"
    } else {
        return "hidden"
    }
}

######################################################################
# Get list of gas indices for dropdown menu
# Returns comma-separated list like "1,2,3" based on number of gases defined
proc ADSIM_2025::GetGasIndices {} {
    set root [customlib::GetBaseRoot]
    if {$root == ""} {
        # Default to single gas if root not available
        return "1"
    }
    
    set xp {//container[@n="materials"]/container[@n="m_gas"]/blockdata}
    set gas_blocks [$root selectNodes $xp]
    set num_gases [llength $gas_blocks]
    
    if {$num_gases == 0} {
        return "1"
    }
    
    # Build comma-separated list of indices
    set indices ""
    for {set i 1} {$i <= $num_gases} {incr i} {
        if {$indices eq ""} {
            set indices "$i"
        } else {
            set indices "$indices,$i"
        }
    }
    
    return $indices
}

proc ADSIM_2025::selectNode {domNode} {
    # Invoke GiD's node selection
    set nodeId [GiD_Info incolor pick nodes]
    
    if {$nodeId ne "" && $nodeId > 0} {
        # Return the node ID to be set as the value
        return $nodeId
    } else {
        WarnWin "No node selected!"
        return ""
    }
}
######################################################################
# Check if a the number of nodes to track and changes the state of the values
# Returns "normal" if node exists, "hidden" if it doesn't
proc ADSIM_2025::CheckNumberofNodes { node_index } {
    set root [customlib::GetBaseRoot]
    if {$root == ""} {
        # If we can't get root, show first node by default
        if {$node_index == 1} {
            return "normal"
        } else {
            return "hidden"
        }
    }
    
    set xp {string(//container[@n="Calculation_Data"]/container[@n="data_saving"]/container[@n="nodes_to_follow"]/value[@n="number_nodes"]/@v)}
    set num_nodes [$root selectNodes $xp]
    
    if {$num_nodes eq ""} {
        set num_nodes 0
    }
    
    if {$node_index <= $num_nodes} {
        return "normal"
    } else {
        return "hidden"
    }
}

######################################################################
# Check if a the number of elements to track and changes the state of the values
# Returns "normal" if node exists, "hidden" if it doesn't
proc ADSIM_2025::CheckNumberofElements { element_index } {
    set root [customlib::GetBaseRoot]
    if {$root == ""} {
        # If we can't get root, show first node by default
        if {$element_index == 1} {
            return "normal"
        } else {
            return "hidden"
        }
    }
    
    set xp {string(//container[@n="Calculation_Data"]/container[@n="data_saving"]/container[@n="elements_to_follow"]/value[@n="number_elements"]/@v)}
    set num_elements [$root selectNodes $xp]
    
    if {$num_elements eq ""} {
        set num_elements 0
    }
    
    if {$element_index <= $num_elements} {
        return "normal"
    } else {
        return "hidden"
    }
}