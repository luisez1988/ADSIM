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

## Analytical references

- **diffusion_1d** — 1D pure diffusion, sine-series (consolidation) solution.
- **advection_1d** — 1D advection–diffusion, Cleary & Adrian (1973) finite-column
  solution with a first-type inlet boundary.

> The analytical `D` (and `v`, `R`) in each case config must match the solver's
> **effective** transport coefficients for that input set, not just the raw
> material properties. If a case starts failing after a physics change, check
> whether the effective coefficient changed before assuming a regression.
