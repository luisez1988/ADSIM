# Changelog

All notable changes to ADSIM will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
