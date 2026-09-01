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

### Scoring an axisymmetric solution on a 2D mesh

`axis = "radius"` scores a 2D case against a solution written in `r`. Each node's
position becomes its distance from `center_x`, `center_y` (both default 0):

```toml
axis     = "radius"
center_x = 0.0
center_y = 0.0
```

Unlike the Cartesian path, this one keeps **every node as its own sample** rather
than averaging nodes that share a position. That is deliberate. On a mesh that is
not itself axisymmetric, the spread between nodes at equal radius *is* the error
an axisymmetric reference is meant to expose; averaging each ring would quietly
discard exactly the signal the case exists to measure.

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

- **conduction_disk_2d** — transient conduction in a solid disk of radius 1 m,
  uniform at 373.15 K, rim held at 273.15 K. Scored against the Fourier–Bessel
  series (Carslaw & Jaeger, §7.6)
  `(T-Ts)/(T0-Ts) = Σ 2/(jₙJ₁(jₙ)) J₀(jₙr/R) exp(-α jₙ² t/R²)`,
  with `J₀`, `J₁` and the zeros of `J₀` implemented in `analytical.jl` (A&S
  9.4.1–9.4.6 plus a Newton polish) so the suite keeps its no-dependency rule.

  **This is the suite's only genuinely 2D case, and it earns its place.** Every
  other case runs on an axis-aligned strip, where the isoparametric Jacobian is
  diagonal and a whole class of mapping errors cancels identically. Here the
  elements are skewed general quads and the reference is axisymmetric, so a
  conduction operator that has acquired a direction dependence cannot hide.

  It was added after finding that `shape_functions.jl` formed `dN/dx` as
  `B * inv(J)` where the chain rule wants `B * transpose(inv(J))`. With the bug
  present, all seven 1D cases returned *bit-identical* error numbers — the suite
  was completely blind to it — while this case went from 0.079 % to 4.007 %
  relative L2, a worst-node error of 28.5 K, and 20 K of spread between nodes at
  equal radius that should agree to 0.05 K. Any change that makes the operator
  direction-dependent will show up here and nowhere else.

### Axisymmetric cases

Four cases cover the axisymmetric formulation of the manuscript appendix, where every
quadrature measure carries a factor `r_gp`. They are deliberately of two kinds, and both
kinds are needed.

- **axisym_diffusion_1d**, **axisym_advection_1d** — the *same* mesh, material file and
  reference as `diffusion_1d` and `advection_1d`, with only `solver_type` changed. Read
  axisymmetrically, those strips are solid cylinders of radius 0.1 m transporting along
  their axis. For a field varying only with `z` the radius factor enters the lumped mass
  and the stiffness row of a node in the same proportion and cancels from the rate, so
  these must reproduce the plane results **exactly**: measured agreement is 5e-15 on
  concentration and identically zero on temperature, and the scored errors match to three
  decimals (0.051 % and 1.707 %).

  They fail if the radius is applied to the mass but not the stiffness, to one flux term
  but not another, or with the wrong coordinate taken as the radius.

- **axisym_conduction_radial**, **axisym_darcy_radial** — genuine radial transport, where
  the radius does *not* cancel. A solid cylinder cooled at its rim, scored against the
  same Fourier–Bessel series as `conduction_disk_2d`; and a well bore feeding gas outward
  through an annulus, whose steady profile has `C²` linear in `ln r`.

  **These are what prove the flag is not being ignored.** A solver that silently ran plane
  strain would pass the two column cases and cannot pass these. Run in plane mode on the
  identical mesh they score 6.11 % and 13.38 % against tolerances of 0.5 % and 2 % —
  46× and 557× worse than the axisymmetric run.

  `axisym_conduction_radial` also covers the axis of symmetry: the mesh reaches `r = 0`
  and carries no boundary condition there, because a surface generated by an edge on the
  axis has zero area. Its axis nodes evolve smoothly (300 → 331.1 → 369.3 K) and remain
  the coldest point, with nothing non-finite anywhere.

### IMPES cases

Four cases — **impes_advection_1d**, **impes_axisym_advection_1d**, **impes_darcy_1d**,
**impes_axisym_darcy_radial** — are paired copies of the four that exercise the
advective/pressure path, with `[solver] time_integration = "IMPES"` and
`advection_stabilization = 1` set and *nothing else changed*: same mesh, same material
file, same analytical reference, same tolerances. A disagreement with the explicit sibling
is therefore attributable to the time integration and to nothing else.

They are pairs rather than replacements because IMPES makes only the advective term
implicit. Every other term, and the whole reaction and energy path, stays explicit and is
already covered by the existing cases.

Two things read differently in an IMPES case:

- **Δt is adaptive.** `courant_number` still scales the result, but the solver recomputes
  the step every step from the velocity field, starting at the explicit-equivalent step and
  growing by at most 1.1× per step. The step count in the *run log* is the number to
  compare, not the one in the start-up report, which reports a limit rather than a step.
- **The gain is in steps, not in tolerance.** `impes_advection_1d` takes 136 steps where
  `advection_1d` takes 29 136, at 1.936 % worst rel-L2 against the explicit 1.707 %. The
  IMPES error is larger only at the first output, where its step is still ramping, and
  smaller at every later one (0.307 % vs 0.516 % at t = 2000 s).

`darcy_1d` / `impes_darcy_1d` are the sharpest of the four: they score the Darcy velocity
directly and so isolate the pressure solve from species transport. Their errors agree to
all printed digits (0.086 %), which is what says the implicit pressure equation reproduces
the explicit one rather than merely resembling it.

See `src/IMPES_FORMULATION_NOTES.md` for the formulation.

`axisym_darcy_radial` scores the **concentration**, not the velocity, and with a
factor-of-two contrast. Both choices are about measuring the formulation rather than
something else. The nodal velocity is projected from the Gauss points, so at the two end
nodes the projection is one-sided and costs ~13 % where `v` varies as `1/r`, while the
interior agrees to better than 0.7 % — an artefact of the point-to-grid transfer, not of
the radius weighting. And with a small contrast the relative L2 would be insensitive: a
4 mol/m³ drop on a background of 1000 scores under 0.4 % even for a completely wrong
profile. `darcy_radial` in `analytical.jl` gives the matching velocity profile for anyone
who wants to inspect it.

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
