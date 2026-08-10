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
lime_content= val
residual_lime= val
lime_arrhenius_factor= val
lime_activation_energy= val
#thermal properties
specific_heat_solids= val

[soil.soil_2]
#physical properties
specific_gravity= val
porosity= val
saturation= val
granular_tortuosity= val
intrinsic_permeability= val
lime_content= val
residual_lime= val
lime_arrhenius_factor= val
lime_activation_energy= val
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
lime_content= val
residual_lime= val
lime_arrhenius_factor= val
lime_activation_energy= val
#thermal properties
specific_heat_solids= val

```


