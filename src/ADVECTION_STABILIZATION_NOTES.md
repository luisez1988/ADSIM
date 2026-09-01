# Streamline-upwind (SU) stabilization of the advective flux

## Why

ADSIM's advective (pressure-diffusion) flux and its thermal-pressure companion term are
discretized as plain consistent-Galerkin FEM (`fully_explicit_solver.jl`, the advective
Gauss-point loop inside `if calculate_advection`), with no upwinding or flux-limiting. On a
mesh where the local element Péclet number is large — as it is at the sharp CO2/Air front of
the carbonation problem — this loses the discrete maximum principle and produces spurious,
non-monotonic oscillation: computed concentrations undershoot below zero, and the solver's
non-negativity clamp silently absorbs it. Investigation (comparing degree-of-carbonation
profiles across mesh sizes, then checking the clamp-firing rate in the run logs) found the
clamp firing on ~1-2% of all nodal updates, every output interval, for the entire simulation
duration on both the `2mm_size` and `4mm_size` carbonation cases — not a startup transient, a
sustained artifact. Critically, CO2 and Air (which never reacts) clamp at nearly identical
rates, ruling out the reaction rate law as the primary driver and implicating the shared
advective discretization.

This is fixed by adding streamline-upwind (SU) stabilization — the variant of SUPG that
stabilizes only the spatial operator, not the mass matrix, which is what keeps it compatible
with this solver's lumped-mass, fully-explicit time integration.

## The mechanism (why plain Galerkin loses monotonicity)

This is a statement about the sign structure of the assembled Galerkin system, not an analogy
to finite differences — it holds on any mesh, in any dimension, for any element order.

For an interior basis function `Nᵢ` (compact support, vanishing on the domain boundary) and the
consistent Galerkin advective bilinear form `Bᵢⱼ = ∫Ω (∇Nᵢ·v) Nⱼ dΩ` (locally-frozen `v`, the
same linearization any local stability diagnostic makes), the divergence theorem gives

```
Bᵢⱼ + Bⱼᵢ = ∫Ω v·∇(NᵢNⱼ) dΩ = ∮∂Ω (NᵢNⱼ)(v·n) dΓ = 0     (Nᵢ = 0 on ∂Ω, for any j)
```

so `Bᵢⱼ = −Bⱼᵢ`: the consistent Galerkin advective operator is exactly skew-symmetric at every
interior node, and in particular `Bᵢᵢ = 0` — it contributes nothing to its own node's diagonal.
(The same follows from the non-conservative form `∫Nᵢ v·∇C dΩ`, so this isn't an artifact of
which weak form is used.) This is the discrete echo of a continuum fact: pure advection
transports a quantity without adding or removing it, so it can never by itself contribute the
kind of diagonal weight that keeps a discrete scheme bounded.

The diffusive stiffness matrix `K`, by contrast, is an M-matrix for these elements: positive
diagonal, non-positive off-diagonal, rows summing to zero. That alone gives the discrete
diffusion equation a maximum principle. Assembling diffusion + advection together, row `i`'s
entry to neighbor `j` is `Kᵢⱼ + Bᵢⱼ`. As long as `|Bᵢⱼ| ≤ |Kᵢⱼ|`, the combined entry stays ≤ 0
and the M-matrix property survives. Once the advective contribution on some edge exceeds the
diffusive one, that entry flips positive, the M-matrix property is lost, and the discrete
solution can overshoot or undershoot its true bounds — oscillate.

**Worked corollary (1D, exact).** For linear elements, uniform spacing `h`, constant `v, D`:
`Kᵢⱼ = −D/h`, `Bᵢⱼ = ±v/2`, so `|Bᵢⱼ| ≤ |Kᵢⱼ|` becomes the mesh Péclet number
`Pe ≡ vh/D ≤ 2`. The discrete nodal equation is a linear recursion
`(1−Pe/2)Cᵢ₊₁ − 2Cᵢ + (1+Pe/2)Cᵢ₋₁ = 0`, with characteristic roots `r=1` and
`r=(2+Pe)/(2−Pe)`. For `Pe<2` the ratio `r>0` (smooth decay); for `Pe>2` it flips negative,
so `rⁱ` alternates sign with `i` — a sawtooth riding the true exponential profile. (Some texts
define `Pe ≡ vh/2D`, giving the equivalent `Pe ≤ 1` form of the same criterion.)

## The fix

Perturb the Galerkin test function in the flow direction, `wᵢ → wᵢ + τ*(v·∇wᵢ)`, applied to the
spatial part of the residual only (full SUPG also perturbs the transient term, which would make
the mass matrix non-diagonal — incompatible with this solver's lumped-mass explicit update).
Keeping only the advective residual against itself:

```
∫_Ω τ*(v·∇wᵢ)(v·∇C) dΩ = ∫_Ω ∇wᵢ · (τ* v⊗v) · ∇C dΩ
```

an anisotropic diffusion term with tensor `τ*·v⊗v` — rank 1, aligned entirely with the flow
direction, exactly zero in the crosswind direction (unlike naive isotropic artificial
viscosity, which would also blur the solution transverse to the flow). Element form:

```
K^SU_ij = τ* ∫_Ωe (v·∇Nᵢ)(v·∇Nⱼ) dΩ
```

### Is this consistent?

Full SUPG (using the complete strong-form residual `R(C) = v·∇C − ∇·(D∇C) − source`) is exactly
consistent for any `h`: substituting the true solution makes the entire added term vanish,
since `R(C_exact) ≡ 0`. Keeping only the advective self-term, as here, does not vanish at the
true solution in general — a genuine, if small, perturbation at finite `h`. Two things make
this the right default anyway:

1. ADSIM uses bilinear C0 quads, and a bilinear field has zero curvature inside each
   undistorted (parallelogram) element (`∇·(D∇C_h) ≡ 0` there exactly, piecewise-constant
   gradient), so the dropped diffusive-residual term is identically zero on undistorted
   elements — the "full residual" and "advective-only" choices coincide exactly there, and
   dropping it is the standard simplification in the SUPG literature for this element order
   regardless (Donea & Huerta 2003, §2.4).
2. `τ* → 0` as `h → 0`, so the scheme is still convergent (asymptotically consistent) even
   where it doesn't coincide exactly.

### The stabilization parameter τ*

Classical "optimal" value (Christie, Griffiths, Mitchell & Zienkiewicz 1976; Heinrich,
Huyakorn, Zienkiewicz & Mitchell 1977), which makes the 1D discrete scheme nodally exact:

```
τ*_opt = (h / 2|v|) · [coth(Pe/2) − 2/Pe]        (→ h/2|v| as Pe→∞, → 0 as Pe→0)
```

Evaluating `coth` every Gauss point is wasteful and numerically delicate near `Pe = 0`. The
standard robust practical approximation (Brooks & Hughes 1982), implemented here:

```
τ* ≈ (h / 2|v|) · min(Pe/3, 1)
```

### Applied to ADSIM's species equation

The existing advective flux for species `g` is already the discrete weak form of `∇·(C_g v_g)`
with `v_g = −(k/μ_g)·R T·∇C_tot` (never materialized as a vector today; made explicit by this
work). Using this species' own effective diffusivity (matching the diffusion term exactly,
`θ_g D_g τ_soil` — `τ_soil` here is granular tortuosity, a material property, unrelated to the
stabilization parameter `τ*`, an unfortunate but unavoidable symbol collision):

```
D_g,eff = θ_g D_g τ_soil,   Pe_g = |v_g| h_e / D_g,eff,   τ*_g = (h_e/2|v_g|)·min(Pe_g/3, 1)
```

`D_g,eff` is set to zero when `diffusion = 0` in `[solver]`, matching what the assembled system
actually contains rather than what the material file lists — `axisym_darcy_radial` in the
verification suite is exactly this case, by construction.

The companion "thermal pressure" term needs its **own, independent** treatment — it shares the
element operator shape but is driven by the opposite gradient (roles of `T` and `C_tot`
swapped, the product-rule split of `∇P = R(T∇C_tot + C_tot∇T)`):

```
v_g^T = −(k/μ_g)·R·C_tot·∇T,   τ*_g,T = (h_e/2|v_g^T|)·min(Pe_g^T/3, 1)
```

`v_g` and `v_g^T` are not parallel in general, so `τ*_g ≠ τ*_g,T`. Axisymmetric formulation
needs no special treatment: `r_gp` enters purely through the shared quadrature measure
(`ShapeFunctions.get_weight`) that every term in the loop already consumes.

## Configuration

Opt-in via `[solver] advection_stabilization = 1` in the calculation TOML (default `0`, so
every existing calc.toml is unaffected). Implemented in
`fully_explicit_solver.jl`'s existing advective Gauss-point loop, gated by
`calculate_stabilization`; see `stab_tau` for the guarded τ* formula.

## References

1. Hughes, T.J.R., Brooks, A.N. (1979). *A multi-dimensional upwind scheme with no crosswind
   diffusion.* In *Finite Element Methods for Convection Dominated Flows*, ASME AMD Vol. 34.
2. Brooks, A.N., Hughes, T.J.R. (1982). *Streamline upwind/Petrov-Galerkin formulations for
   convection dominated flows...* Computer Methods in Applied Mechanics and Engineering,
   32(1-3), 199-259.
3. Donea, J., Huerta, A. (2003). *Finite Element Methods for Flow Problems.* Wiley.
4. Christie, I., Griffiths, D.F., Mitchell, A.R., Zienkiewicz, O.C. (1976). *Finite element
   methods for second order differential equations with significant first derivatives.*
   International Journal for Numerical Methods in Engineering, 10(6), 1389-1396.

Optional, for later refinement: Tezduyar, T.E., Osawa, Y. (2000). *Finite element
stabilization parameters computed from element matrices and vectors.* CMAME, 190(3-4),
411-430 — a more modern, robust multi-scale τ formula.
