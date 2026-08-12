# TOML writing instructions for calculation parameters

This document shows the structure of the TOML file for the calculation parameters, units, and simulation settings.

## Header

The file contains "commented" info about the version of ADSIM used.
```toml
# ADSIM calculation parameters file header (need to add disclaimer and license)
# ADSIM_version = "2025 v0.x.x"
# File_created = "2025-01-01"
```

## Units and Dimensions

The following section defines the unit system used in the simulation.
```toml
[units]
geometry_unit = "m"
mass_unit = "kg"
pressure_unit = "kPa"
temperature_unit = "K"
time_unit = "s"
```

## Gravity Vector

Define gravity acceleration magnitude and components.
```toml
[gravity]
gravity_magnitude = 9.81
gravity_x_component = 0.0
gravity_y_component = -1.0
```

## Solver Settings

Define the solver type for the problem.
```toml
[solver]
solver_type = "2D-Plane"
```

## Time Stepping Data

Parameters for time-dependent calculations.
```toml
[time_stepping]
total_simulation_time = 1.0
time_per_step = 0.01
courant_number = 0.98
```

## Data Saving and Probing

Settings for data output and monitoring specific nodes/elements.
```toml
[data_saving]
data_saving_interval = 1.0

[probing]
# Nodes to follow
number_of_nodes = 3
nodes_to_probe = [1, 2, 3]

# Elements to follow
number_of_elements = 2
elements_to_probe = [1, 2]
```

### Node probes

Every node listed in `nodes_to_probe` is followed through the run and written as
its own time series:

```
output/<project>_node<ID>_stage<N>.csv
```

The file is a plain-text CSV with a commented preamble (project, stage, node id,
node coordinates, sampling interval and unit labels) followed by one header row
and one row per sample. Columns carry the same names as the VTK point data, so a
probe series and a VTK snapshot of the same node agree field for field:

| column                              | meaning                                  |
| ----------------------------------- | ---------------------------------------- |
| `time`                              | Simulation time of the sample            |
| `Total_Concentration`               | Summed gas concentration                 |
| `Absolute_Pressure`                 | Pressure from the ideal gas law          |
| `Reaction_Rate`                     | dC_lime/dt                               |
| `Lime_Concentration`                | C_lime                                   |
| `CaCO3_Concentration`               | C_caco3                                  |
| `Degree_of_Carbonation`             | C_caco3 / C_caco3,max                    |
| `Volumetric_Binder_Content`         | Binder content                           |
| `Temperature`, `Temperature_Rate`   | T and dT/dt                              |
| `<Gas>_Concentration`               | One column per gas species               |
| `<Gas>_Concentration_Rate`          | One rate column per gas species          |
| `Gas_Seepage_Velocity_X/Y[/Z]`      | Darcy velocity components                |

`data_saving_interval` is the **probe sampling interval**. It is set in GiD under
Calculation Data → Simulation probing, in the same container as the node list,
and is independent of the VTK output cadence (which follows `time_per_step`), so
nodes can be sampled far more finely than the mesh is written.

Sampling does not affect the solution. A sample is taken at the first time step
at or past each scheduled probe time and stamped with the true current time,
rather than shortening `dt` to land on it exactly; the recorded times therefore
sit at or just after the exact multiples of the interval. The schedule advances
in whole multiples of the interval, so it does not drift over long runs.

Each calculation stage writes its own file. A restarted run opens
`..._stage2.csv` beginning at the restart time, whose first row repeats the final
state of `..._stage1.csv`.

Node ids outside the mesh, and repeated ids, are reported in the log and skipped
rather than aborting the run. If `number_of_nodes` disagrees with the length of
`nodes_to_probe`, the listed ids win and the mismatch is reported.

Probing is optional: with `nodes_to_probe = []`, or without a `[probing]` section
at all, no probe files are written.

> `elements_to_probe` is parsed but not yet consumed by the solver. Only node
> time series are written at present.

## Complete Example

```toml
# ADSIM calculation parameters file header (need to add disclaimer and license)
# ADSIM_version = "2025 v0.x.x"
# File_created = "2025-01-01"

[units]
geometry_unit = "m"
mass_unit = "kg"
pressure_unit = "kPa"
temperature_unit = "K"
time_unit = "s"

[gravity]
gravity_magnitude = 9.81
gravity_x_component = 0.0
gravity_y_component = -1.0

[solver]
solver_type = "2D-Plane"

[time_stepping]
total_simulation_time = 100.0
time_per_step = 0.1
courant_number = 0.98

[data_saving]
data_saving_interval = 1.0

[probing]
# Nodes to follow during simulation
number_of_nodes = 3
nodes_to_probe = [10, 25, 50]

# Elements to follow during simulation (data at Gaussian points)
number_of_elements = 2
elements_to_probe = [5, 15]
```
