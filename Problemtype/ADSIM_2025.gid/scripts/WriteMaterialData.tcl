#===============================================================================
# ADSIM Material Data Writer (TOML format)
# This script writes material properties for gases, liquids, and soils
# to a .toml file for use in the FEM code
#===============================================================================

proc ADSIM::WriteMaterialData { filename } {
    # Initialize the file
    GiD_WriteCalculationFile init $filename
    
    # Get root of XML tree
    set root [$::gid_groups_conds::doc documentElement]
    
    # Write header
    ADSIM::WriteMaterialHeader
    
    # Write all dictionaries first (TOML requirement: all top-level keys before tables)
    ADSIM::WriteGasDictionary $root
    ADSIM::WriteSoilDictionary $root
    GiD_WriteCalculationFile puts ""
    
    # Write gas properties
    ADSIM::WriteGasProperties $root
    
    # Write liquid properties
    ADSIM::WriteLiquidProperties $root
    
    # Write soil properties
    ADSIM::WriteSoilProperties $root
    
    # Close the file
    GiD_WriteCalculationFile end
}

#===============================================================================
# Write material file header
#===============================================================================
proc ADSIM::WriteMaterialHeader { } {
    GiD_WriteCalculationFile puts "# ADSIM material file header (need to add disclaimer and license)"
    GiD_WriteCalculationFile puts "# ADSIM_version = \"2025 v0.x.x\""
    GiD_WriteCalculationFile puts "# File_created = \"[clock format [clock seconds] -format %Y-%m-%d]\""
    GiD_WriteCalculationFile puts ""
}

#===============================================================================
# Write gas dictionary (must be at top level before any tables)
#===============================================================================
proc ADSIM::WriteGasDictionary { root } {
    # Get all gas materials
    set xp_gases {//container[@n="materials"]/container[@n="m_gas"]/blockdata}
    set gas_blocks [$root selectNodes $xp_gases]
    
    if {[llength $gas_blocks] == 0} {
        return
    }
    
    # Build gas dictionary array
    set gas_names [list]
    foreach gas_block $gas_blocks {
        set gas_name [$gas_block @name]
        lappend gas_names "\"$gas_name\""
    }
    
    # Write gas dictionary as array
    GiD_WriteCalculationFile puts "gas_dictionary_ = \[[join $gas_names ", "]\]"
}

#===============================================================================
# Write gas properties
#===============================================================================
proc ADSIM::WriteGasProperties { root } {
    # Get all gas materials
    set xp_gases {//container[@n="materials"]/container[@n="m_gas"]/blockdata}
    set gas_blocks [$root selectNodes $xp_gases]
    
    if {[llength $gas_blocks] == 0} {
        return
    }
    
    # Write gas properties
    GiD_WriteCalculationFile puts "# Gas properties"
    foreach gas_block $gas_blocks {
        set gas_name [$gas_block @name]
        GiD_WriteCalculationFile puts "\[gas.\"$gas_name\"\]"
        
        # Dynamic viscosity
        set dyn_visc [$gas_block selectNodes {string(value[@n="Dynamic_viscosity_"]/@v)}]
        GiD_WriteCalculationFile puts "dynamic_viscosity = $dyn_visc"
        
        # Molar mass
        set molar_mass [$gas_block selectNodes {string(value[@n="molar_mass_"]/@v)}]
        GiD_WriteCalculationFile puts "molar_mass = $molar_mass"
        
        # Diffusion coefficient
        set diff_coef [$gas_block selectNodes {string(value[@n="diffusion_coefficient_1"]/@v)}]
        GiD_WriteCalculationFile puts "diff_coefficient = $diff_coef"
        
        GiD_WriteCalculationFile puts ""
    }
}

#===============================================================================
# Write liquid properties
#===============================================================================
proc ADSIM::WriteLiquidProperties { root } {
    # Get liquid properties container
    set xp_liquid {//container[@n="materials"]/container[@n="m_liquid"]}
    set liquid_container [$root selectNodes $xp_liquid]
    
    if {$liquid_container == ""} {
        return
    }
    
    GiD_WriteCalculationFile puts "# Liquid properties"
    GiD_WriteCalculationFile puts "\[liquid\]"
    
    # Dynamic viscosity
    set dyn_visc [$liquid_container selectNodes {string(value[@n="Dynamic_viscosity_"]/@v)}]
    GiD_WriteCalculationFile puts "dynamic_viscosity = $dyn_visc"
    
    # Density
    set density [$liquid_container selectNodes {string(value[@n="density_"]/@v)}]
    GiD_WriteCalculationFile puts "density = $density"
    
    # Specific heat
    set spec_heat [$liquid_container selectNodes {string(value[@n="specific_heat_water_"]/@v)}]
    GiD_WriteCalculationFile puts "specific_heat = $spec_heat"
    
    GiD_WriteCalculationFile puts ""
}

#===============================================================================
# Write soil dictionary (must be at top level before any tables)
#===============================================================================
proc ADSIM::WriteSoilDictionary { root } {
    # Get all soil materials
    set xp_soils {//container[@n="materials"]/container[@n="m_soil"]/blockdata}
    set soil_blocks [$root selectNodes $xp_soils]
    
    if {[llength $soil_blocks] == 0} {
        return
    }
    
    # Build soil dictionary array
    set soil_names [list]
    foreach soil_block $soil_blocks {
        set soil_name [$soil_block @name]
        lappend soil_names "\"$soil_name\""
    }
    
    # Write soil dictionary as array
    GiD_WriteCalculationFile puts "soil_dictionary_ = \[[join $soil_names ", "]\]"
}

#===============================================================================
# Write soil properties
#===============================================================================
proc ADSIM::WriteSoilProperties { root } {
    # Get all soil materials
    set xp_soils {//container[@n="materials"]/container[@n="m_soil"]/blockdata}
    set soil_blocks [$root selectNodes $xp_soils]
    
    if {[llength $soil_blocks] == 0} {
        return
    }
    
    # Write soil properties
    GiD_WriteCalculationFile puts "# Soil properties"
    foreach soil_block $soil_blocks {
        set soil_name [$soil_block @name]
        GiD_WriteCalculationFile puts "\[soil.\"$soil_name\"\]"
        
        # Get the first group for this soil (assumes all groups have same values)
        set groups [$soil_block selectNodes {condition[@n="soil_basic"]/group}]
        
        if {[llength $groups] > 0} {
            set first_group [lindex $groups 0]
            
            # Physical properties
            GiD_WriteCalculationFile puts "# Physical properties"
            
            set spec_grav [$first_group selectNodes {string(.//value[@n="specific_gravity"]/@v)}]
            GiD_WriteCalculationFile puts "specific_gravity = $spec_grav"
            
            set porosity [$first_group selectNodes {string(.//value[@n="initial_porosity_"]/@v)}]
            GiD_WriteCalculationFile puts "porosity = $porosity"
            
            set saturation [$first_group selectNodes {string(.//value[@n="saturation_"]/@v)}]
            GiD_WriteCalculationFile puts "saturation = $saturation"
            
            set tort [$first_group selectNodes {string(.//value[@n="granular_tortuosity_"]/@v)}]
            GiD_WriteCalculationFile puts "granular_tortuosity = $tort"
            
            set perm [$first_group selectNodes {string(.//value[@n="intrinsic_permeabiliy_"]/@v)}]
            GiD_WriteCalculationFile puts "intrinsic_permeability = $perm"
            
            set lime [$first_group selectNodes {string(.//value[@n="lime_content_"]/@v)}]
            GiD_WriteCalculationFile puts "lime_content = $lime"

            set reaction_rate [$first_group selectNodes {string(.//value[@n="reaction_rate_"]/@v)}]
            GiD_WriteCalculationFile puts "lime_reaction_rate = $reaction_rate"
            
            set res_lime [$first_group selectNodes {string(.//value[@n="lime_impure_"]/@v)}]
            GiD_WriteCalculationFile puts "residual_lime = $res_lime"
            
            # Thermal properties
            GiD_WriteCalculationFile puts "# Thermal properties"
            
            set spec_heat_solid [$first_group selectNodes {string(.//value[@n="specific_heat_solid_"]/@v)}]
            GiD_WriteCalculationFile puts "specific_heat_solids = $spec_heat_solid"
        }
        
        GiD_WriteCalculationFile puts ""
    }
}
