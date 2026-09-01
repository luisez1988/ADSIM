# Mesh file writing instructions      

The mesh in ADSIM is written in plain text files with extension `.mesh`. The structure follows the following order:

## Header
This indicates the version of the ADSIM used.
## Counters
The counters are written in the following format:

``` 
MESH numNodes Ntotal numElements NelmTotal
```
`Ntotal` is the total number of nodes

`NelmTotal` is the total number of elements

## Nodal coordinates

The file must be structured in the following manner:
```
coordinates
0.0 0.0
1.0 0.0
1.0 1.0
0.0 1.0
end coordinates
```
The coordinates must be wrapped on the `coordinates` and `end coordinates` text. The coordinates follow and the number of coordinates must be exactly `Ntotal`

## Elements

The file must be written in the following manner:
```
elements
1 2 3 4
4 5 6 7
2 5 6 3
end elements
```

The nodes ids corresponding that define the element conductivities are wrapped around `elements` and `end elements`. The order of the information implicitly defines the element id.

## Concentration boundary conditions

Use the following structure:
```
concentration_bc
counter
.
.
.
node_id value_gas1 value_gas2 ... value_gasn
.
.
.
end concentration_bc
```

The structure is wrapped as shown. Each line contains the node id and gas concentrations assigned to that node id. Eg. `2 0.5 0.25` indicates that node 2 has a concentration of 0.5 $\text{mol}/\text{L}^3$ of gas 1 and 0.25 $\text{mol}/\text{L}^3$ of gas 2. The `counter` is needed to initialize the arrays that will contain the data.

## Uniform flow boundary conditions
Use the following structure:
```
uniform_flow_bc
counter
.
.
.
node_id value_gas1 value_gas2 ... value_gasn
.
.
.
end uniform_flow_bc
```

Each line contains the node id within a line where the condition was applied. After this, the value of normal flow is set, but it will be recalculated in real time to consider length of elements. Eg. `2 1.0 0.0` indicates that node 2 has a flow of 0.5 $\text{mol}/\text{L}^2 \text{T}$ of gas 1 and 0.25 $\text{mol}/\text{L}^2 \text{T}$ of gas 2.

## Absolute pressure boundary condition

Use the following structure:
```
absolute_pressure
counter
.
.
.
NodeId value
.
.
.
ebd absolute_pressure
```
Here the value of pressure is assigned to the `NodeId`. Eg, `1 250` indicates node 1 has an absolute pressure of 250  $\text{F} /\text{L}^2$

## Initial gas concentrations

Use the following structure:
```
initial_concentrations
counter
.
.
.
ElmId value_gas1 value_gas2 ... value_nasn
.
.
.
end initial_concentrations
```
E.g., `5 0.5 0.2` indicates element 5 has an initial concentration on its nodes of 0.5 and 0.2 $\text{mol}/\text{L}^3$ for gas 1 and 2 respectively.

## Initial temperature

use the following structure:
```
initial_temperature
counter
.
.
.
ElmId value
.
.
.
end initial_temperature
```
E.g., `5 295` indicates element 5 has an initial temperature on its nodes of 295 $\text{T}_p$. Later, a subroutine will assign temperatures to nodes.

## Material assignation
Use the following structure:
```
material
.
.
.
ElmID MatIndex
.
.
.
end material
```
The structure consists of a list of elements with the material assigned to each. E.g., `12 1` indicates material index 1 is assigned to element 12.

## Time-dependent boundary conditions

Any of the five boundaries below can be driven by a *time function*: a named
curve `f(t)` defined once in the calculation TOML and referenced from the mesh
file by integer id. Id `0` means "no curve": the boundary holds the constant
value from its own block, which is what a mesh file without any of these extra
blocks gives.

The ids live in blocks that mirror the value block they belong to, one line per
node, in the same node order and covering the same nodes:

| block | payload | drives |
|---|---|---|
| `concentration_bc_tf` | `NodeID tf_1 ... tf_nGas` | `concentration_bc` |
| `uniform_flow_bc_tf` | `NodeID tf_1 ... tf_nGas` | `uniform_flow_bc` |
| `absolute_pressure_tf` | `NodeID tf` | `absolute_pressure` |
| `temperature_bc_tf` | `NodeID tf` | `temperature_bc` |
| `partial_pressure_bc_tf` | `NodeID tf_1 ... tf_nGas` | `partial_pressure_bc` |

Each block is optional, and so is any node inside it. For example:

```
temperature_bc
2
1 350
4 350
end temperature_bc

temperature_bc_tf
2
1 2
4 0
end temperature_bc_tf
```

drives node 1 with curve 2 and leaves node 4 fixed at 350.

The ids are 1-based positions of the `[[time_function]]` entries in
`<project>_calc.toml`, in file order. Both files are written from the same GiD
project in one pass, so they always agree; the solver stops with an error if a
mesh file references an id the calculation file does not define.

### Ramping a partial pressure boundary (`partial_pressure_bc_tf`)

Like `concentration_bc_tf`/`uniform_flow_bc_tf`, this block is per-gas: one
time-function id per gas, in `gas_dictionary` order. To ramp CO2 up from the
domain's initial condition instead of snapping it on at full strength from
step 1 (which is hard on the smallest elements nearest the boundary — see
`src/ADVECTION_STABILIZATION_NOTES.md`'s discussion of the top-corner clamp
cluster for why this matters), use a `mode = "absolute"` ramp curve whose
`v_start` matches the domain's initial partial pressure for that gas and whose
`v_end` matches the constant value already in `partial_pressure_bc`:

```
partial_pressure_bc
7
225 113000 0.0
...
end partial_pressure_bc

partial_pressure_bc_tf
7
225 1 0
...
end partial_pressure_bc_tf
```

with, in `<project>_calc.toml`:

```toml
[[time_function]]
id = 1
name = "CO2 partial pressure ramp"
type = "ramp"
mode = "absolute"
t_start = 0.0
t_end = 60.0          # tune to the mesh - long enough to blunt the initial shock
v_start = 0.0          # matches the domain's initial CO2 partial pressure
v_end = 113000.0       # matches the constant already in partial_pressure_bc
```

Curve id `0` for Air (the second entry on the `partial_pressure_bc_tf` line)
leaves Air's boundary value constant, since it starts at 0 and has nothing to
ramp from.

### Curve definitions

Curves are `[[time_function]]` tables in `<project>_calc.toml`. Every curve has
an `id`, a `name`, a `type` and a `mode`:

- `mode = "absolute"` - the curve **is** the boundary value, in the units of the
  boundary condition, and the number in the mesh file is ignored.
- `mode = "multiplier"` - the applied boundary is (mesh file value) x f(t), with
  f dimensionless.

```toml
# Linear ramp, held flat outside [t_start, t_end]
[[time_function]]
id = 1
name = "co2_ramp"
type = "ramp"
mode = "absolute"
t_start = 0.0
t_end = 3600.0
v_start = 0.0
v_end = 1.2

# mean + amplitude * wave(2*pi*(t - t_start)/period + phase)
# wave = sine | square | triangle | sawtooth, all peaking at +/- amplitude
[[time_function]]
id = 2
name = "daily_ambient"
type = "periodic"
mode = "absolute"
wave = "sine"
mean = 293.15
amplitude = 10.0
period = 86400.0
phase = 0.0
t_start = 0.0

# Lookup in a two column CSV of time and value
[[time_function]]
id = 3
name = "measured_flux"
type = "table"
mode = "absolute"
file = "measured_flux.csv"
interpolation = "linear"   # linear | step (zero order hold)
before = "hold"            # hold | zero | repeat | error
after = "hold"
time_offset = 0.0          # added to every time in the file, after scaling
time_scale = 1.0           # multiplies every time; use 3600 for a file in hours
average = false            # see below
```

### CSV files

Two columns, time first. Blank lines and lines beginning with `#` are skipped,
and a header row is detected automatically. Comma, semicolon, tab and plain
whitespace all separate. Times must be strictly increasing; a file that breaks
that is rejected at load with the offending line number.

```
# measured chamber pressure
time,pressure
0.0,101325.0
3600.0,102000.0
```

`file` is resolved relative to the folder holding the calculation file, so the
CSV travels with the input set. The problemtype copies it into the generated
`.ADSIM` folder automatically.

### Time step and file spacing

A table is a continuous function of `t`, not a sequence of steps to march
through, and the lookup is a stateless binary search. Nothing couples the file's
sample spacing to the solver's time step, and a checkpoint restart resumes at
the right point of every curve.

- **Time step finer than the file** (the usual case): the interval between two
  samples is sub-sampled exactly. Nothing to do.
- **Time step coarser than the file**: samples in between are never visited and
  the boundary is an aliased version of the data. The solver warns once at
  start-up, naming the curve. For a *flow* boundary set `average = true`, which
  integrates the curve over each step instead of sampling it at one instant and
  keeps the total injected mass exactly right. For a prescribed value, reduce
  the time step instead.

Note that the critical time step is computed once from the initial state, so a
curve that drives a pressure well above its starting value may need a lower
Courant number to stay stable for the whole run.
