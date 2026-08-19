# ADSIM Verification Suite

Automated regression/verification tests. Each case runs the ADSIM solver on a
fixed input set, then compares the numerical result against a **closed-form
analytical solution** and reports **PASS/FAIL** on the relative L2 error.

Run it after making changes to the solver:

```bash
julia verification/run_verification.jl              # run every case
julia verification/run_verification.jl diffusion_1d # run selected case(s)
julia verification/run_verification.jl --verbose    # also show solver output
julia verification/run_verification.jl --list        # list available cases
```

The command exits `0` if all cases pass and `1` otherwise, so it can be dropped
into CI or a pre-commit hook.

## No more manual ParaView exports

The old workflow required opening each result in ParaView, applying *Plot Over
Line*, and exporting a CSV by hand before anything could be compared. That step
is gone. ADSIM already writes every field into its `.vtk` snapshots, so the
suite parses those files directly ([`src/vtk_reader.jl`](src/vtk_reader.jl)) and
reconstructs the 1D profile itself (`line_profile`). Nothing manual, fully
repeatable.

## How a run works

For each case the runner ([`src/verify.jl`](src/verify.jl)):

1. copies the case's input files into `src/data/`;
2. runs the solver exactly as a user would (`buildscripts/run_cli.jl <project>`);
3. reads the resulting `src/output/<project>_NNNNNN.vtk` snapshots;
4. for every saved step (t = 0 skipped), extracts the field along the 1D axis
   and evaluates the analytical solution at the same node positions;
5. computes the **relative L2 error** `‖num − exact‖ / ‖exact‖`
   ([`src/metrics.jl`](src/metrics.jl)) and also reports RMSE and max-abs;
6. the case **passes** if the worst relative L2 across all checked steps is
   within the case's `tolerance`.

## Layout

Everything for a case lives in one folder — `verification/cases/<case>/` — so
there is a single site per case, no mirrored trees.

```
verification/
  run_verification.jl     # the command you run
  src/
    vtk_reader.jl         # parses ADSIM .vtk (replaces ParaView CSV export)
    analytical.jl         # closed-form solutions + self-contained erfc
    metrics.jl            # relative L2 / rmse / max-abs
    verify.jl             # orchestration
  cases/
    <case>/
      case.toml           # this case's config
      <project>.mesh      # ADSIM inputs (committed)
      <project>_calc.toml
      <project>_mat.toml
      gid/                # GiD project (provenance / to regenerate the mesh)
```

The inputs and `case.toml` are `*.toml`/`*.mesh` files that the repo's global
`.gitignore` would normally hide, so `.gitignore` re-includes everything under
`verification/cases/**` — do not remove those rules or the suite stops being
reproducible. A case folder is only picked up if it contains a `case.toml`
(e.g. `reaction/` currently holds inputs but no config yet, so it is skipped).

## Adding a case

1. Make a folder `verification/cases/<case>/` and drop the ADSIM inputs in it
   (`<project>.mesh`, `<project>_calc.toml`, `<project>_mat.toml`).
2. If it needs a new analytical model, add a branch to
   `Analytical.evaluate` in [`src/analytical.jl`](src/analytical.jl).
3. Add `verification/cases/<case>/case.toml`:

```toml
project        = "My_test"              # base name of the input files in this folder
field          = "CO2_Concentration"    # VTK scalar to compare
axis           = "y"                     # 1D axis of the mesh
boundary_coord = 1.0                     # axis location of the driving boundary
tolerance      = 0.02                    # PASS if worst relative L2 <= 2%
check_steps    = []                      # [] = all saved steps except t=0

[analytical]
type = "diffusion_series"                # dispatch key in Analytical.evaluate
# ... model parameters ...
```

### Scoring a node time series instead of a profile

Some cases are a history at one point rather than a profile at each instant — a
0-D kinetics test has no spatial variation to plot. Those are scored from the
node probe CSV instead of the VTK snapshots:

```toml
source     = "probe"                 # default is "vtk"
probe_node = 1                       # which node; must appear in [probing]
field      = "Degree_of_Carbonation" # a probe CSV column, not a VTK field
```

The calculation file has to actually request that node, or there is nothing to
read:

```toml
[probing]
number_of_nodes = 1
nodes_to_probe = [1]
```

Probe rows are written at `data_saving_interval`, a separate cadence from the VTK
`time_per_step`. The comparison drops `t = 0` and evaluates the analytical
solution at each sample time with `x = 0`, so the analytical type must ignore its
position argument (`carbonation_isothermal` and `darcy_velocity` both do). Column
names come from the probe header — `time`, `Degree_of_Carbonation`,
`Lime_Concentration`, `Temperature`, `<gas>_Concentration` and so on; a wrong name
reports the available ones.

## Analytical references

- **diffusion_1d** — 1D pure diffusion, sine-series (consolidation) solution.
- **advection_1d** — 1D advection–diffusion, Cleary & Adrian (1973) finite-column
  solution with a first-type inlet boundary.
- **darcy_1d** — steady Darcy flow. A single gas fills a 1 m column with fixed
  concentrations at both ends; since the gas is ideal (`P = C R T`), the
  concentration difference `dC` sets a pressure gradient and the code should
  return the Darcy flux `v = (K/mu) R T dC / L`. The test compares the magnitude
  of the computed `Gas_Seepage_Velocity` at steady state against that value.
- **reaction_0d** — carbonation kinetics with no transport. One element with every
  node on a fixed-concentration boundary, so the CO2 field never moves, and
  `reaction_enthalpy = 0` so the run is isothermal and the rate coefficient is
  constant. The rate law then integrates in closed form to
  `DoC(t) = (1 - p_r)(1 - exp(-lambda t))` with
  `lambda = k_T a theta_w K_H R T C_co2` — Eq. (5.4) of
  `calibration_reactionrate/integrator_mathematics.md`, exact to round-off. The
  residual is therefore ADSIM's forward-Euler lime update against the true
  exponential, a first-order error of order `lambda*dt/2`. Halving
  `courant_number` halves it, which is the real check: 0.056 % → 0.028 % →
  0.014 % at Courant 0.98 → 0.49 → 0.245.

  Scored on `Degree_of_Carbonation` from the node probe. Do not give this case a
  nonzero reaction enthalpy without also giving it a thermal boundary condition —
  the reference stops being valid the moment `T` moves.

### Comparing a vector field (e.g. velocity)

`field` may name a VTK **vector** (like `Gas_Seepage_Velocity`) instead of a
scalar. Add `component = "magnitude"` (or `"x"`/`"y"`/`"z"`) to pick what to
compare. Use `snapshot = "last"` to score only the final, steady-state snapshot
instead of every step — handy when the quantity is meaningful only once the run
has settled.

> The analytical `D` (and `v`, `R`) in each case config must match the solver's
> **effective** transport coefficients for that input set, not just the raw
> material properties. If a case starts failing after a physics change, check
> whether the effective coefficient changed before assuming a regression.
