# Changelog

All notable changes to ADSIM will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed
- **Critical time step is now measured from the assembled operator rather than assumed
  from a closed form.** `calculate_critical_time_step` evaluated `λ_max` of `M̂⁻¹K`
  analytically as `4·coef/h_min²`, which is exact on a uniform *plane* quad mesh and only
  there. Two errors followed from it, and together they made the 2 mm axisymmetric
  carbonation run unstable at the shipped `courant_number = 0.9`:
  - `get_maximum_total_concentration` counted the initial field and `concentration_bc`
    but not the pressure-specified boundaries, which prescribe a concentration just as
    surely. The benchmark starts at 41 mol/m³ of air while its injection face holds CO₂ at
    113 kPa — 46.4 mol/m³ — so the advective bound was **13.1%** too permissive.
  - The closed form is blind to the axisymmetric measure. Both `M̂` and `K` carry
    `r_gp |det J|`, and away from the axis the radius cancels between them; in the element
    touching `r = 0` it does not. Measured against a dense eigensolve, `λ_max` on that mesh
    is **3.3%** above the plane value, and its eigenvector carries **83% of its energy on
    the axis column** — which is why the resulting instability appeared as a sawtooth
    against the axis rather than across the domain.

  Combined, the step was 16.8% too long: an effective Courant number of **1.051**, just
  over the stability boundary. The `4 mm` mesh came out at 0.871 and was stable, so the
  defect only surfaced on refinement. The clamp at `C_MIN` arrested the growing mode every
  step, which is why it presented as a bounded, plausible-looking field instead of a NaN.

  `λ_max` is now obtained by power iteration on the operator the solver actually
  assembles, reusing `assemble_lumped_mass_vector!` and
  `assemble_element_stiffness_matrices` — both of which already carry the axisymmetric
  measure. The closed forms are retained and cross-checked, and the step is taken as
  `min(measured, analytic)` so the measurement can only ever shorten it. Cost is ~50 ms
  once at startup.

  **This changes every axisymmetric result.** `2mm_size` moves from `dt = 9.17e-4` to
  `7.06e-4` s; plane runs and all twelve verification cases are unchanged. Any figure
  already generated from an axisymmetric run needs regenerating. The measurement also
  catches non-uniform *plane* meshes: `conduction_disk_2d` was running 10% over its true
  stability limit and is now corrected.
- Verification gained **relative L∞ and total-variation** metrics alongside relative L2,
  with TV gated by default at `1.15x` the exact profile (`tv_tolerance` per case). L2
  averages a localized sawtooth away across the profile; TV does not. All twelve cases
  score between 1.000 and 1.009.
- The non-negativity clamp now **counts** its firings and reports them per load step,
  warning when they exceed 0.1% of nodal updates. It previously warned once per gas for a
  whole run, so an unstable mode being arrested thousands of times a step read as a single
  rounding artefact at step 1.

### Added
- `--report-timestep` runs setup through the time step calculation and stops before the
  solver, printing the per-mechanism stability spectrum — measured `λ_max`, the closed-form
  estimate, their ratio, and the resulting limit. Answers "what will my step be and which
  term sets it" in seconds on any mesh, without committing to the run.
- `src/test/test_operator_spectrum.jl`, guarding the above: the power iteration against a
  dense eigensolve, the axisymmetric penalty being real, and the invariant that the step
  the solver uses never exceeds `2/λ_max`.
- Axisymmetric formulation, implementing the manuscript appendix. Selecting
  `solver_type = "2D-Axisymmetric"` treats the mesh as the meridian section of a body of
  revolution: `x` is the radius `r`, `y` is the axis `z`, and every domain and boundary
  integral gains a factor `r`. The dropdown already existed in GiD and the value already
  reached the solver — it was read into `solver_settings["dimension"]` and only printed.
  It now selects the formulation.
  - The whole change is carried by one cached array. `ShapeFunctionData` gained
    `weight[e,p]`, holding `|det J|` under plane strain and `r_gp |det J|` when
    axisymmetric, evaluated once at initialization from `r_gp = N_gp · R_nodes`. Every
    quadrature sum reads `get_weight` instead of `get_detJ`, so both formulations run the
    same code with no branch in the time loop and no extra arithmetic per step. That
    includes the cached element stiffness, and with it the diffusive and thermal
    conductive fluxes.
  - Boundary integrals use the ring measure `l_e * r̄_e` in place of the edge length. This
    also makes the external flux vanish on an edge lying on the axis with no special case,
    since that surface has zero area, so the axis needs no boundary condition.
  - Gravity must act along the axis, so `gravity_x_component` must be `0.0`. A radial
    component would vary with the azimuthal angle and contradict axisymmetry, so this is
    an error rather than a warning. A negative radius and an unknown `solver_type` are
    likewise rejected, and a checkpoint records its formulation and will not resume under
    the other one.
  - The stability limit and the reaction kinetics are unchanged, as the appendix derives.
  - Four verification cases added: `axisym_diffusion_1d` and `axisym_advection_1d` reuse
    the plane meshes and must reproduce them exactly (measured 5e-15), while
    `axisym_conduction_radial` and `axisym_darcy_radial` carry genuine radial transport and
    fail by 46× and 557× if run in plane mode — which is what shows the flag is acted on.

### Changed
- The lumped mass vector is now assembled as Eq. (lumped_mass_matrix) states it,
  `M_a = Σ_gp θ_g N_a(gp) w_gp W_gp`, instead of computing the element area and splitting
  it equally as `θ_g A_e / 4`. The two agree exactly for parallelograms — every 1-D strip
  in the verification suite — and differ for general quads, where the equal split gives
  every node the same mass regardless of element shape. Verified against independent 40×40
  Gauss integration: exact to 1e-16 on a parallelogram, a trapezoid, and axisymmetric
  elements both away from and touching the axis. The full plane suite is unchanged (8/8,
  matching the previous outputs to round-off).

### Fixed
- The gravity term had an inverted sign, in both places it appears. The solver was
  integrating `v = -(K/μ_g)(∇p + ρ_g g)` rather than Darcy's law, Eq. (Darcy_law),
  `v = -(K/μ_g)(∇p - ρ_g g)`. Since the gravity vector is supplied by the user already
  oriented — `gravity_y_component = -1.0` in every input file, i.e. pointing down — the
  two combined to drive a dense gas *upward*, and left a hydrostatic gas column out of
  equilibrium. The gravitational nodal flux now carries the leading minus of
  Eq. (gravitational_flux), and the gravitational Darcy velocity is `+(K/μ_g) ρ_g g`.
  The two must move together, and did.
  - The magnitude of the term was already correct and is unchanged: the Figure-3
    verification measured it against `|ρ_g g| / |∇p|` and agreed to 0.016 percentage
    points, with only the sign inverted. Expect results from any run with `gravity = 1`
    to change by roughly twice the size of the gravitational contribution.
  - `Gas_Seepage_Velocity` in the VTK output and the sensible heat carried by the gas
    both read the same Gauss-point velocity and are corrected with it.
  - Every shipped verification case sets `gravity = 0`, so the whole suite is bit-for-bit
    unchanged — which is also why nothing caught this. The sign convention of the input
    vector is now stated in `Problemtype/Calculation_Parameters_TOML_Instructions.md`,
    where previously only its orientation constraint under axisymmetry was documented.

### Added
- Node probing is implemented. Nodes selected in GiD under Calculation Data →
  Simulation probing are now followed through the run and written as one CSV time
  series per node, `output/<project>_node<ID>_stage<N>.csv`, holding every nodal
  state variable under the same names the VTK output uses. The node list was
  previously written to the calculation file and parsed by the solver, but never
  used.
  - `data_saving_interval` is the probe sampling interval. It was read into an
    unused local before; it is the field that already sits beside the node list
    in the GiD interface, and it is independent of the VTK cadence, so nodes can
    be sampled much more finely than the mesh is dumped.
  - Sampling does not affect the solution: a sample is taken at the first step at
    or past each scheduled probe time and stamped with the true current time,
    rather than shortening `dt` to land on it. The schedule advances in whole
    multiples of the interval, so it does not drift.
  - Each stage writes its own file; a restarted run begins at the restart time.
  - Invalid or repeated node ids are reported and skipped instead of aborting.
- Reaction enthalpy ΔH_r is now an input, set in GiD under Materials → Reactant
  properties and written to a `[reactants]` table of the material file. It was
  previously hardcoded at -113800 J/mol CO₂ in the solver.
- The solver logs the reaction enthalpy and the Arrhenius pair in use, and warns
  when the enthalpy is entered as a positive (endothermic) number, which would
  cool the soil instead of heating it.

### Fixed
- The standalone executable can be built again. `create_app` needs a package, so
  the repository now carries a `Project.toml` declaring the name, uuid and
  dependencies, and `src/ADSIM.jl` defines the `julia_main()` entry point that
  PackageCompiler calls. Neither existed, so `julia build_app.jl` could not have
  succeeded.
- The precompile execution file no longer targets `Reaction_test`, a project with
  no mesh in `src/data`, which made `main()` exit non-zero and abort the build. It
  now traces a case that exists, and does so inside a temporary workspace, so a
  build no longer writes VTK files, probe series or a checkpoint into
  `src/output/` (the checkpoint would have pushed the next real run into a new
  stage).

### Changed
- The carbonation reaction inputs are separated from the soil properties, and
  split by what they describe:
  - `[reactants]` holds the reaction itself — enthalpy, Arrhenius factor and
    activation energy — and applies to every material. In GiD this is the new
    Materials → Reactant properties page.
  - `[reaction."<name>"]` holds the mix design of one soil material, its lime
    content and residual lime, so different layers can still carry different
    amounts of lime. In GiD this is a "Reaction properties" page on each soil.

  All four keys previously sat in `[soil."<name>"]` alongside the physical
  properties.
- Material files written against the old layout still run: the reaction keys are
  read from the `[soil]` tables with a deprecation warning, and a missing
  `[reactants]` table falls back to the previous built-in enthalpy, so results
  are unchanged. A legacy file whose soils declare different kinetics cannot be
  represented by the single Arrhenius pair; the first material's values are used
  and the conflict is reported.

## [0.1.0] - 2025-11-30

### Added

#### Core Solver Features
- Fully explicit finite element solver for advection-diffusion-reaction problems in porous media
- Support for 2D and 3D unstructured meshes (triangular and tetrahedral elements)
- Multi-gas transport simulation with individual gas species tracking
- Advection, diffusion, and reaction kinetics modeling
- Dynamic gas field implementation for spatially varying properties
- Temperature-dependent transport properties
- Adaptive time-stepping with CFL condition enforcement

#### Reaction Modeling
- First-order reaction kinetics
- Carbonation reaction modeling (Ca(OH)₂ + CO₂ → CaCO₃)
- Degree of carbonation tracking
- Lime and calcium carbonate concentration tracking
- Heat generation from exothermic reactions

#### Input/Output
- TOML-based input format for materials and calculation parameters
- Custom mesh file format (.mesh) with boundary condition support
- VTK output for visualization in ParaView
- Comprehensive logging system with execution time tracking
- Initial and boundary condition specifications

#### GiD Problem Type
- Full GiD integration for pre-processing
- Interactive mesh generation interface
- Material property editor
- Calculation parameter configuration
- Automated TOML and mesh file generation from GiD
- Splash screen and about dialog

#### Build System
- PackageCompiler.jl integration for standalone executable generation
- Automated build script (`build_app.jl`)
- Cross-platform Julia dependencies management

#### Version Management
- Centralized VERSION file at project root
- Automated version update script (`scripts/update_version.jl`)
- Version display in log files, splash screen, and VTK output
- `--version` command-line flag support

#### Verification & Testing
- Advection test case
- Diffusion test case
- Reaction test case
- Multi-gas verification problem
- Pani et al. benchmark problem

#### Documentation
- Comprehensive README with project overview
- Materials TOML writing instructions
- Mesh file writing instructions
- Calculation parameters TOML instructions
- Dynamic gas fields implementation guide
- Fully explicit solver documentation
- Time step module documentation
- External assets guide
- Python script for ParaView data extraction

### Technical Details
- **Language**: Julia 1.8+
- **License**: MIT
- **Author**: Luis E. Zambrano-Cruzatty
- **Solver Type**: Fully explicit finite element method
- **Element Types**: Linear triangular (2D), linear tetrahedral (3D)
- **Shape Functions**: Implemented with local-to-global coordinate transformation

### Known Limitations (Alpha Release)
- Explicit time-stepping requires small time steps for stability
- Limited to linear elements only
- No adaptive mesh refinement
- Windows executable only (Linux/macOS builds not yet automated)
- GiD problem type requires GiD 15.0 or higher
- Documentation still in progress for some advanced features

### Dependencies
- Julia >= 1.8
- TOML.jl
- Dates (Julia stdlib)
- Printf (Julia stdlib)
- PackageCompiler.jl (for building)

---

## Release Notes for v0.1.0

This is the first alpha release of ADSIM (Advection-Diffusion for Soil Improvement and Modification). The software is functional for simulating gas transport and reaction in porous media but is still under active development.

**What works:**
- Core solver for multi-gas transport with reactions
- GiD pre-processor integration
- VTK output for ParaView visualization
- Several validation test cases

**What's in progress:**
- User documentation and tutorials
- Additional validation cases
- Performance optimization
- Multi-platform build automation

**Feedback Welcome:**
This is an alpha release intended for testing and feedback. Please report issues, bugs, or feature requests via GitHub Issues.

---

[Unreleased]: https://github.com/luisez1988/ADSIM/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/luisez1988/ADSIM/releases/tag/v0.1.0
