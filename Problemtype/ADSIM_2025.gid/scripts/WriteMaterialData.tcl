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

    # Write reactant properties (chemistry shared by all materials)
    ADSIM::WriteReactantProperties $root

    # Write soil properties
    ADSIM::WriteSoilProperties $root

    # Write per-material reaction properties
    ADSIM::WriteReactionProperties $root

    # Close the file
    GiD_WriteCalculationFile end
}

#===============================================================================
# Write material file header
#===============================================================================
proc ADSIM::WriteMaterialHeader { } {
    GiD_WriteCalculationFile puts "# ADSIM material file header (need to add disclaimer and license)"
    GiD_WriteCalculationFile puts "# ADSIM_version = \"2025 v0.1.0\""
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

        # Thermal properties
        set gas_conductivity [$gas_block selectNodes {string(value[@n="thermal_conductivity_gas_"]/@v)}]
        # Models saved before this field existed return an empty string.
        if {$gas_conductivity == ""} { set gas_conductivity 0.0 }
        GiD_WriteCalculationFile puts "thermal_conductivity = $gas_conductivity"

        set molar_cp [$gas_block selectNodes {string(value[@n="molar_heat_capacity_"]/@v)}]
        if {$molar_cp == ""} { set molar_cp 0.0 }
        GiD_WriteCalculationFile puts "molar_heat_capacity = $molar_cp"
        
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

    # Thermal conductivity
    set water_conductivity [$liquid_container selectNodes {string(value[@n="thermal_conductivity_water_"]/@v)}]
    # Models saved before this field existed return an empty string.
    if {$water_conductivity == ""} { set water_conductivity 0.0 }
    GiD_WriteCalculationFile puts "thermal_conductivity = $water_conductivity"
    
    GiD_WriteCalculationFile puts ""
}

#===============================================================================
# Write reactant properties
#
# Chemistry and kinetics of the carbonation reaction, common to every soil
# material. The per-material lime quantities are written separately by
# WriteReactionProperties.
#===============================================================================
proc ADSIM::WriteReactantProperties { root } {
    # Get reactant properties container
    set xp_reactants {//container[@n="materials"]/container[@n="m_reactants"]}
    set reactants_container [$root selectNodes $xp_reactants]

    if {$reactants_container == ""} {
        return
    }

    GiD_WriteCalculationFile puts "# Reactant properties"
    GiD_WriteCalculationFile puts "\[reactants\]"

    # Reaction enthalpy
    set reaction_enthalpy [$reactants_container selectNodes {string(value[@n="reaction_enthalpy_"]/@v)}]
    GiD_WriteCalculationFile puts "reaction_enthalpy = $reaction_enthalpy"

    # Reaction kinetics
    set arrhenius_factor [$reactants_container selectNodes {string(value[@n="arrhenius_factor_"]/@v)}]
    GiD_WriteCalculationFile puts "lime_arrhenius_factor = $arrhenius_factor"

    set activation_energy [$reactants_container selectNodes {string(value[@n="activation_energy_"]/@v)}]
    GiD_WriteCalculationFile puts "lime_activation_energy = $activation_energy"

    # Interfacial-area coefficient of the same rate law. Zero disables the factor.
    # Models saved before this field existed return an empty string, so fall back to
    # zero rather than writing an unparsable line.
    set interfacial_area_beta [$reactants_container selectNodes {string(value[@n="interfacial_area_beta_"]/@v)}]
    if {$interfacial_area_beta == ""} {
        set interfacial_area_beta 0.0
    }
    GiD_WriteCalculationFile puts "lime_interfacial_area_beta = $interfacial_area_beta"

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

            # Thermal properties
            GiD_WriteCalculationFile puts "# Thermal properties"
            
            set spec_heat_solid [$first_group selectNodes {string(.//value[@n="specific_heat_solid_"]/@v)}]
            GiD_WriteCalculationFile puts "specific_heat_solids = $spec_heat_solid"

            set solid_conductivity [$first_group selectNodes {string(.//value[@n="thermal_conductivity_solid_"]/@v)}]
            # Models saved before this field existed return an empty string.
            if {$solid_conductivity == ""} { set solid_conductivity 0.0 }
            GiD_WriteCalculationFile puts "thermal_conductivity_solids = $solid_conductivity"
        }

        GiD_WriteCalculationFile puts ""
    }
}

#===============================================================================
# Write reaction properties
#
# One table per soil material, holding how much lime that material carries.
# Keyed by the same soil name, so the solver pairs each reaction block with its
# material. The kinetics are common to all materials and are written by
# WriteReactantProperties instead.
#===============================================================================
proc ADSIM::WriteReactionProperties { root } {
    # Get all soil materials
    set xp_soils {//container[@n="materials"]/container[@n="m_soil"]/blockdata}
    set soil_blocks [$root selectNodes $xp_soils]

    if {[llength $soil_blocks] == 0} {
        return
    }

    # Write reaction properties
    GiD_WriteCalculationFile puts "# Reaction properties"
    foreach soil_block $soil_blocks {
        set soil_name [$soil_block @name]
        GiD_WriteCalculationFile puts "\[reaction.\"$soil_name\"\]"

        # Get the first group for this soil (assumes all groups have same values)
        set groups [$soil_block selectNodes {condition[@n="soil_basic"]/group}]

        if {[llength $groups] > 0} {
            set first_group [lindex $groups 0]

            set lime [$first_group selectNodes {string(.//value[@n="lime_content_"]/@v)}]
            GiD_WriteCalculationFile puts "lime_content = $lime"

            set res_lime [$first_group selectNodes {string(.//value[@n="lime_impure_"]/@v)}]
            GiD_WriteCalculationFile puts "residual_lime = $res_lime"
        }

        GiD_WriteCalculationFile puts ""
    }
}
