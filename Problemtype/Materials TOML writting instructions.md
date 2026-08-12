# TOML writting instructions for material data

This document shows the structure of the TOML file for the material properties.

## Header

The file contains "commented" info about the version of ADSIM used.
```toml
# ADSIM material file header (need to add disclaimer and license)
# ADSIM_version = "2025 v0.x.x"
# File_created = "2025-01-01"
```

## material properties
The following is the structure of the gas properties.
```toml
gas_dictionary_=[gas_1, gas_2, ... gas_n]

#gas properties
[gas.gas_1]
dynamic_viscosity= val
molar_mass= val
diff_coeffiecient= val

[gas.gas_2]
dynamic_viscosity= val
molar_mass= val
diff_coeffiecient= val

.
.
.

[gas.gas_n]
dynamic_viscosity= val
molar_mass= val
diff_coeffiecient= val

#liquid properties
[liquid]
dynamic_viscosity= val
density= val
specific_heat= val

#reactant properties (chemistry and kinetics of the carbonation reaction,
#shared by all soils)
[reactants]
reaction_enthalpy= val
lime_arrhenius_factor= val
lime_activation_energy= val

#soil dictionary
soil_dictionary= [soil_1, soil_2, ... soil_n]

#soil properties
[soil.soil_1]
#physical properties
specific_gravity= val
porosity= val
saturation= val
granular_tortuosity= val
intrinsic_permeability= val
#thermal properties
specific_heat_solids= val

[soil.soil_2]
#physical properties
specific_gravity= val
porosity= val
saturation= val
granular_tortuosity= val
intrinsic_permeability= val
#thermal properties
specific_heat_solids= val
.
.
.
[soil.soil_n]
#physical properties
specific_gravity= val
porosity= val
saturation= val
granular_tortuosity= val
intrinsic_permeability= val
#thermal properties
specific_heat_solids= val

#reaction properties, one table per soil material, keyed by the same name
[reaction.soil_1]
lime_content= val
residual_lime= val

[reaction.soil_2]
lime_content= val
residual_lime= val
.
.
.
[reaction.soil_n]
lime_content= val
residual_lime= val

```

## Reaction inputs

The carbonation reaction

```
Ca(OH)2(s) + CO2(g) -> CaCO3(s) + H2O(l)
```

is described by two groups of inputs, kept apart from the soil block so that the
soil tables hold only physics that applies whether or not reaction kinetics is
switched on.

`[reactants]` describes the reaction itself, so it is the same for every material:

| key                      | meaning                                                    |
| ------------------------ | ---------------------------------------------------------- |
| `reaction_enthalpy`      | ΔH_r per mole of CO2 [J/mol]. **Negative**, since carbonation is exothermic. |
| `lime_arrhenius_factor`  | k_o in k_T = k_o exp(-E/RT) [m³/(mol·s)]                    |
| `lime_activation_energy` | E, activation energy [J/mol]                               |

`[reaction."<soil name>"]` holds the mix-design inputs of one soil material,
matched to the soil of the same name, so different layers can carry different
amounts of lime:

| key             | meaning                                                  |
| --------------- | -------------------------------------------------------- |
| `lime_content`  | β_l, lime as a fraction of the weight of solids [-]      |
| `residual_lime` | Non-reactive lime, as a fraction of the lime content [-] |

Material files written before this split carry all four keys inside their
`[soil]` tables. Those are still read, with a warning, and `[reactants]` defaults
to -113800 J/mol when absent, so older files reproduce their previous results
unchanged. Because one Arrhenius pair now serves every material, a legacy file
whose soils declare *different* kinetics cannot be represented exactly: the first
material's pair is used and the conflict is reported.


