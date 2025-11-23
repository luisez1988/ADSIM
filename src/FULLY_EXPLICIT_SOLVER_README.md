# Fully Explicit Solver for Gas Diffusion

## Overview

This document describes the implementation of the fully explicit finite element solver for transient gas diffusion in porous media within ADSIM.

## Governing Equation

The transient diffusion of gases in porous media is governed by:

```
∂(θ_g C_g)/∂t = ∇ · (D_eff ∇C_g)
```

where:
- `θ_g = n(1 - S_r)` = gas volume fraction [-]
- `C_g` = gas concentration [mol/m³]
- `D_eff = D_g × τ` = effective diffusion coefficient [m²/s]
- `D_g` = gas diffusion coefficient [m²/s]
- `τ` = granular tortuosity [-]
- `n` = porosity [-]
- `S_r` = degree of saturation [-]

## Finite Element Formulation

### Semi-discrete Form

After spatial discretization using the Galerkin finite element method:

```
M dC/dt = F
```

where:
- `M` = lumped mass matrix (diagonal)
- `C` = nodal concentration vector
- `F` = diffusion flow vector = `-K × C`
- `K` = stiffness matrix

### Mass Matrix (Lumped)

The lumped mass matrix is diagonal:

```
M_i = ∫_Ω θ_g N_i dΩ ≈ θ_g × A_e / 4
```

For each element, the contribution is distributed equally to its 4 nodes.

### Stiffness Matrix

The stiffness matrix represents the diffusion operator:

```
K_ij = ∫_Ω D_eff ∇N_i · ∇N_j dΩ
```

Evaluated using 2×2 Gaussian quadrature (4 Gauss points).

## Time Integration

### Forward Euler (Explicit) Scheme

```
C^(n+1) = C^n + Δt × M^(-1) × F^n
C^(n+1) = C^n + Δt × M^(-1) × (-K × C^n)
```

Since M is diagonal (lumped), the inversion is trivial:

```
dC_i/dt = F_i / M_i
C_i^(n+1) = C_i^n + Δt × (F_i / M_i)
```

### Stability Criterion

The explicit scheme is conditionally stable. The time step must satisfy:

```
Δt ≤ Δt_crit
```

where `Δt_crit` is calculated from:

```
Δt_diffusion = (h_min² × τ) / (θ_g × D_max)
```

The actual time step used is:

```
Δt = C_N × Δt_crit
```

where `C_N` is the Courant number (typically 0.5-1.0).

## Algorithm Implementation

### Preprocessing (Once per simulation)

1. **Assemble lumped mass vector M**
   - Loop over all elements
   - Calculate element area using Gauss quadrature
   - Distribute mass contribution to nodes

2. **Assemble stiffness matrices K** (one per gas)
   - Loop over all elements
   - Calculate element stiffness using Gauss quadrature
   - Assemble into global matrix

### Time Stepping Loop

For each time step `n`:

For each gas species `g` (parallelized with `@threads`):

1. **Compute diffusion flow**
   ```julia
   F = -K × C^n
   ```

2. **Compute rate of change**
   ```julia
   dC/dt = F / M  (element-wise)
   ```

3. **Apply boundary conditions**
   - Zero rate at fixed concentration nodes
   - Maintain BC values

4. **Update concentrations**
   ```julia
   C^(n+1) = C^n + Δt × dC/dt
   ```

5. **Enforce non-negativity**
   ```julia
   C = max(C, 0)
   ```

6. **Write output** (if at output interval)

## Boundary Conditions

### Concentration Boundary Conditions

At nodes with fixed concentrations:
- Rate is set to zero: `dC/dt = 0`
- Concentration is maintained at BC value

### Implementation

```julia
for (node_id, concentrations) in mesh.concentration_bc
    dC_dt[node_id] = 0.0
    C[node_id] = concentrations[gas_idx]
end
```

## Parallelization

The solver leverages Julia's multi-threading capabilities:

1. **Gas species loop**: Each gas can be solved independently in parallel
2. **Element assembly**: Element loops use `@threads` where safe
3. **Matrix operations**: Uses optimized BLAS routines via `LinearAlgebra`

To enable threading, run Julia with:
```bash
julia -t auto kernel.jl TEST
```

## File Structure

- `fully_explicit_solver.jl` - Main solver implementation
  - `assemble_lumped_mass_vector!()` - Mass matrix assembly
  - `assemble_diffusion_matrix!()` - Stiffness matrix assembly
  - `compute_diffusion_flow!()` - Flow computation
  - `apply_boundary_conditions!()` - BC application
  - `fully_explicit_diffusion_solver()` - Main solver routine
  - `write_output_vtk()` - Output writing wrapper

## Usage Example

```julia
# In kernel.jl, after initialization:
fully_explicit_diffusion_solver(mesh, materials, calc_params, time_data, log_print)
```

## Testing

Run the test problem:

```bash
cd src
julia kernel.jl TEST
```

This will:
1. Read `data/TEST.mesh`, `data/TEST_mat.toml`, `data/TEST_calc.toml`
2. Solve the diffusion problem for CO2 and Air
3. Write VTK output files to `output/` directory
4. Generate a log file at `output/TEST.log`

## Verification

The solver can be verified against analytical solutions for:

1. **1D diffusion in semi-infinite domain** (complementary error function solution)
2. **1D diffusion with constant concentration boundaries** (steady-state linear profile)
3. **2D diffusion from point source** (fundamental solution)

## Future Extensions

Current implementation: **Diffusion only**

Planned extensions:
- [ ] Advection terms (Darcy flow)
- [ ] Chemical reactions (lime carbonation)
- [ ] Heat transfer
- [ ] Adaptive time stepping
- [ ] Implicit/semi-implicit schemes for stiff problems
- [ ] Nonlinear diffusion coefficients

## Performance Considerations

### Memory Usage

- Mass vector `M`: `O(N)` where N = number of nodes
- Stiffness matrix `K`: `O(N²)` - stored as dense matrix
- For large problems, consider sparse matrix storage

### Computational Cost

Per time step:
- Matrix-vector product: `O(N²)` operations
- For `NGases` gases and `Nsteps` time steps: `O(NGases × Nsteps × N²)`

### Optimization Tips

1. Use sparse matrices for large problems (future improvement)
2. Enable multi-threading for gas loop
3. Reduce output frequency for long simulations
4. Use appropriate Courant number (0.5-1.0)

## References

1. **Finite Element Method**
   - Zienkiewicz, O.C., Taylor, R.L., "The Finite Element Method", Vol. 1-3

2. **Diffusion in Porous Media**
   - Bear, J., "Dynamics of Fluids in Porous Media"
   - Millington, R.J., Quirk, J.P., "Permeability of porous solids"

3. **Time Integration**
   - Reddy, J.N., "An Introduction to the Finite Element Method"

## Troubleshooting

### Common Issues

1. **Negative concentrations**
   - Check Courant number (reduce if needed)
   - Verify boundary conditions
   - Check material properties (positive diffusion coefficients)

2. **Unstable solution (oscillations)**
   - Reduce Courant number
   - Check time step calculation
   - Verify mesh quality (no distorted elements)

3. **Slow convergence**
   - Enable multi-threading
   - Check output frequency
   - Consider reducing mesh refinement

## Contact

For questions or issues, contact: Luis Zambrano-Cruzatty
