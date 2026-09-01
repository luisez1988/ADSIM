# IMPES: implicit treatment of the advective term

Companion to `ADVECTION_STABILIZATION_NOTES.md` and `Instructions_coding.md`. Derives the
IMPES (IMplicit Pressure, Explicit Species) split implemented in `impes_solver.jl` and
selected by `[solver] time_integration = "IMPES"`.

## 1. Why

The advective term sets the time step on every carbonation run. From
`output/2mm_size_timestep_report.log` on the 2 mm mesh:

| Mechanism | λ measured | λ analytic | Δt [s] |
|---|---|---|---|
| Diffusive | 1.6435e+00 | 1.5916e+00 | 1.2169e+00 |
| **Advective** | **2.5463e+03** | **2.4659e+03** | **7.8545e-04** |
| Reactive | — | — | 6.2251e-01 |

Critical Δt = 7.85e-4 s, actual 7.07e-4 s, **1 698 640 steps** for 1200 s of physical time.
The advective operator is 1550× tighter than diffusion and 790× tighter than the reaction.
Enabling SU stabilization roughly halves the step again, because `coef_su ≤ coef_advection`
(Cauchy–Schwarz, see `build_element_coefficients` in `time_step.jl`) and the two enter the
harmonic sum separately.

The reason is structural. In ADSIM's formulation the advective flux is not written as
$\mathbf{v}\cdot\nabla C_\alpha$; it is written as a **diffusion operator acting on the
total concentration**, because the driving pressure is itself a function of that total:

$$
q^{a}_{\alpha,i} \;=\; \int_\Omega \frac{R\,k\,C_\alpha\,T}{\mu_\alpha}\,
\nabla N_i \cdot \nabla C_t \; d\Omega ,
\qquad C_t = \sum_\beta C_\beta .
$$

The bracketed coefficient is a **pressure diffusivity**

$$
D_P \;=\; \frac{R\,k\,C_\alpha\,T}{\mu_\alpha} \;=\; \frac{k\,P}{\mu_\alpha}\,\frac{C_\alpha}{C_t} .
$$

Backing it out of the measured step, $D_P = h^2/(2\Delta t) = (2.01\times10^{-3})^2 /
(2\cdot 7.85\times10^{-4}) \approx 2.6\times10^{-3}\ \mathrm{m^2/s}$, against a molecular
$\theta_g \tau D_g \approx 1.6\times10^{-6}\ \mathrm{m^2/s}$ — three orders of magnitude.
Forward Euler on a diffusion operator costs $\Delta t \propto h^2/D_P$, hence the step.

**Key point.** This is a stiff *feedback loop*, not a fast transport process: a
perturbation in $C_\alpha$ changes $C_t$, hence $P$, hence $\mathbf v$, hence $C_\alpha$
again. Physically the loop is pressure equilibration, which is fast and which we do not
need to resolve in time. Breaking it is exactly what IMPES does.

## 2. Notation

| Symbol | Meaning | Code |
|---|---|---|
| $\alpha,\beta$ | gas species index, $1\ldots N_g$ | `gas_idx` |
| $C_\alpha$ | species molar concentration [mol/m³ **gas**] | `C_g[:, α]` |
| $C_t=\sum_\beta C_\beta$ | total concentration | `total_concentration` |
| $P = R\,T\,C_t$ | absolute pressure (ideal gas; **diagnostic**) | `P` |
| $\theta_g = n(1-S_r)$, $\theta_w = n S_r$ | gas / water volume fractions (constant per element) | `props.θ_g` |
| $k$, $\mu_\alpha$, $\tau$, $D_\alpha$ | permeability, viscosity, tortuosity, diffusivity | `props.permeability`, … |
| $M_i = \sum_e\sum_p \theta_g N_i(\xi_p)\, w_p W_p$ | lumped mass | `M`, `assemble_lumped_mass_vector!` |
| $w_p$ | quadrature measure, $\lvert\det J\rvert$ or $r_{gp}\lvert\det J\rvert$ | `ShapeFunctions.get_weight(e,p)` |
| $n$ (superscript) | time level | |

## 3. Governing equations as currently assembled

Strong form for species $\alpha$:

$$
\frac{\partial (\theta_g C_\alpha)}{\partial t}
= \nabla\!\cdot\!\big(\theta_g \tau D_\alpha \nabla C_\alpha\big)
- \nabla\!\cdot\!\big(C_\alpha \mathbf v_\alpha\big)
+ \theta_w r_\alpha ,
$$

$$
\mathbf v_\alpha = -\frac{k}{\mu_\alpha}\big(\nabla P - \rho_g \mathbf g\big),
\qquad
\nabla P = R\big(T\,\nabla C_t + C_t\,\nabla T\big).
$$

Lumped-mass Galerkin semi-discretization (Eq. `FEM_compact`; `fully_explicit_solver.jl`
rate assembly):

$$
M_i \dot C_{\alpha,i}
= q^{b}_{\alpha,i} - q^{d}_{\alpha,i} - q^{a}_{\alpha,i}
- q^{T}_{\alpha,i} - q^{g}_{\alpha,i} + q^{r}_{\alpha,i}
$$

with, summed over elements and the four Gauss points,

$$
\begin{aligned}
q^{d}_{\alpha,i} &= \theta_g \tau D_\alpha \big(\mathbf K\, C_\alpha\big)_i
&&\text{diffusion}\\[2pt]
q^{a}_{\alpha,i} &= \int \frac{R k C_\alpha T}{\mu_\alpha}\,\nabla N_i\cdot\nabla C_t \,d\Omega
&&\text{advection}\\[2pt]
q^{T}_{\alpha,i} &= \int \frac{R k C_\alpha C_t}{\mu_\alpha}\,\nabla N_i\cdot\nabla T \,d\Omega
&&\text{thermal pressure}\\[2pt]
q^{g}_{\alpha,i} &= -\int \frac{k C_\alpha \rho_g}{\mu_\alpha}\,\nabla N_i\cdot \mathbf g \,d\Omega
&&\text{gravity}\\[2pt]
q^{su}_{\alpha,i} &= \int \tau^{*}_{\alpha}\,
(\mathbf v_\alpha\!\cdot\!\nabla N_i)(\mathbf v_\alpha\!\cdot\!\nabla C^h_\alpha)\,d\Omega
&&\text{SU stabilization}
\end{aligned}
$$

`q^su` is accumulated into `q_aux` alongside `q^a` in the code; it is kept separate here
because IMPES treats the two differently. The reaction enters only for CO₂, as
$q^{r}_{i} = -M_i(\theta_w/\theta_g)\,r_i$.

## 4. The IMPES split

### 4.1 Step 1 — implicit equation for $C_t$

Evaluate $\nabla C_t$ at level $n{+}1$ in $q^{a}$ and lag every coefficient at level $n$.
Each species equation now contains the shared unknown field $C_t^{n+1}$ but advances only
its own $C_\alpha$, so the single closed equation available for $C_t^{n+1}$ is obtained by
**summing the species equations already being solved**. Using

$$
\sum_\alpha \frac{R k C_\alpha^{n} T^{n}}{\mu_\alpha}
\;=\; R\,k\,T^{n} \sum_\alpha \frac{C_\alpha^{n}}{\mu_\alpha}
\;\equiv\; \Lambda^{n}
\qquad \textbf{(total mobility)}
$$

gives one scalar equation for one scalar field:

$$
\boxed{\;
\Big(\tfrac{1}{\Delta t}\mathbf M + \mathbf A^{n}\Big) C_t^{n+1}
= \tfrac{1}{\Delta t}\mathbf M\, C_t^{n} + \mathbf b^{n}
\;}
\tag{IMPES-P}
$$

$$
\big(\mathbf A^{n} u\big)_i = \int \Lambda^{n}\,\nabla N_i\cdot\nabla u \,d\Omega ,
\qquad
b^{n}_i = \sum_\alpha\Big[\,q^{b}_{\alpha,i} - q^{d,n}_{\alpha,i} - q^{T,n}_{\alpha,i}
- q^{g,n}_{\alpha,i} - q^{su,n}_{\alpha,i}\Big] + q^{r,n}_i .
$$

**This is not a new governing equation.** It is the sum of the mass balances already in
the code, and $P = R T C_t$ remains a diagnostic recovered afterwards, exactly as at
present. No new primary unknown, no new boundary-condition type, no new checkpoint state:
the implicit unknown is literally `sum(C_g, dims=2)`.

$\mathbf M$ is diagonal and positive; $\mathbf A^{n}$ is the cached geometric operator
scaled by $\Lambda^{n}\ge 0$, hence symmetric positive semi-definite. The system matrix is
therefore **SPD**, and backward Euler on it is **unconditionally stable** — the $h^2/D_P$
limit is gone by construction.

### 4.2 Step 2 — explicit species transport

With $C_t^{n+1}$ known, each species advances exactly as today, with $\nabla C_t^{n+1}$
substituted for $\nabla C_t^{n}$ **in the advective flux only**:

$$
\frac{M_i\big(C^{n+1}_{\alpha,i} - C^{n}_{\alpha,i}\big)}{\Delta t}
= q^{b}_{\alpha,i} - q^{d,n}_{\alpha,i} - q^{a}_{\alpha,i}\!\big(C_t^{n+1}\big)
- q^{T,n}_{\alpha,i} - q^{g,n}_{\alpha,i} - q^{su,n}_{\alpha,i}
+ q^{r,n}_{\alpha,i} - \lambda_i w_{\alpha,i}
\tag{IMPES-S}
$$

$$
q^{a}_{\alpha,i}\!\big(C_t^{n+1}\big)
= \int \frac{R k C^{n}_\alpha T^{n}}{\mu_\alpha}\,\nabla N_i\cdot\nabla C_t^{n+1}\, d\Omega .
$$

The assembly is the *same loop* as the present code; only the vector being contracted
changes from `C_t` to the just-solved `C_t^{n+1}`.

### 4.3 Exact volume consistency

Sum (IMPES-S) over $\alpha$. Every explicit term reassembles $b^n$ by definition, and the
advective terms collapse through the identity that defined $\Lambda^n$:

$$
\sum_\alpha q^{a}_{\alpha,i}\!\big(C_t^{n+1}\big)
= \int \Big(\sum_\alpha \tfrac{R k C^n_\alpha T^n}{\mu_\alpha}\Big)
\nabla N_i\cdot\nabla C_t^{n+1}\,d\Omega
= \big(\mathbf A^{n} C_t^{n+1}\big)_i .
$$

So the sum of (IMPES-S) **is** (IMPES-P), and therefore at every node where the total row
was not replaced,

$$
\sum_\alpha C^{n+1}_{\alpha,i} \;=\; C^{n+1}_{t,i}
\qquad\text{exactly, to round-off.}
$$

There is **no IMPES volume-balance drift** in this model. The property holds because the
phase is single and the closure $P = RTC_t$ is linear in $C_t$; black-oil IMPES, where the
closure is a nonlinear volume constraint, does not enjoy it. Two conditions must be
honoured in code or it is lost:

1. Steps 1 and 2 must use the **same Gauss-point values** of $C^n_\alpha$, $T^n$, $k$,
   $\mu_\alpha$ and the same quadrature measure. This is why $\mathbf A^n$ is assembled
   per Gauss point (element-averaging $\Lambda$ to reuse the cached `K_e` would break it).
2. $q^{su}$ contracts $\nabla C_\alpha$, not $\nabla C_t$, so it is per-species and
   belongs in $b^n$ as $\sum_\alpha q^{su,n}_{\alpha,i}$ — **not** as an implicit operator
   on $C_t$.

The reaction throttle and the `C_MIN` clamp break the identity when they fire. That is
intended and unchanged; the residual is monitored as a diagnostic.

### 4.4 Boundary conditions

The existing dictionaries map directly onto (IMPES-P). No new BC type is introduced.

| Existing BC | Row treatment in (IMPES-P) |
|---|---|
| `absolute_pressure_bc[j]` | Dirichlet $C^{n+1}_{t,j} = P_{bc}(t^{n+1})/(R T_j)$ — the same target already imposed on `total_concentration[j]` |
| every species prescribed at $j$ | Dirichlet $C^{n+1}_{t,j} = \sum_\alpha C^{bc}_\alpha$ |
| some species prescribed at $j$ | free — natural condition |
| `uniform_flow_bc` | enters $b^n$ through $\sum_\alpha q^{b}_{\alpha,i}$; no row change |

Dirichlet rows are applied by **symmetric elimination** (zero the row and the column, move
the known contribution to the right-hand side) so the matrix stays SPD for CG.

### 4.5 The Lagrange pressure BC is retained unchanged

At an `absolute_pressure_bc` node the total row of (IMPES-P) was *replaced*, so §4.3 does
not apply there: the species must be redistributed to match the pinned total. The required
correction is

$$
\lambda_i = \frac{M_i}{\Delta t}\Big(\textstyle\sum_\alpha \tilde C^{n+1}_{\alpha,i}
- C^{n+1}_{t,i}\Big),
\qquad
w_{\alpha,i} = \frac{P^{b}_{\alpha,i}\,C^{n}_{\alpha,i}}
{\sum_\beta P^{b}_{\beta,i}\,C^{n}_{\beta,i}} ,
$$

where $\tilde C^{n+1}$ is the uncorrected (IMPES-S) update and $P^{b}$ is the
`P_boundary` free-gas gate. This is **algebraically identical** to the corrected explicit
solver: $M_j\big(C^n_t + \Delta t\,\dot C_{t,j} - C^{target}_t\big)/\Delta t$ is the same
residual, and `free_conc_sum` / `n_free_nodes` already implement $w_\alpha$ including the
`C_TOL` degenerate fallback. It is ported verbatim, not replaced.

## 5. Stability of what remains explicit

Linearize (IMPES-S) about $C_\alpha$ with $C_t^{n+1}$ frozen:

$$
\frac{\delta q^{a}_{\alpha,i}}{\delta C_\alpha}
= \int \nabla N_i \cdot \Big(\tfrac{R k T}{\mu_\alpha}\nabla C_t^{n+1}\Big)\,\delta C_\alpha \,d\Omega
= -\int \nabla N_i \cdot \mathbf v_\alpha \, \delta C_\alpha \,d\Omega ,
\qquad
\mathbf v_\alpha = -\frac{kRT}{\mu_\alpha}\nabla C_t^{n+1}.
$$

A **first-order advection operator**, not a diffusion one. The remaining explicit limits
are

$$
\Delta t_{\text{diff}} = \frac{h^2}{2 D_{\max}\tau},
\qquad
\Delta t_{\text{adv}} = \frac{\theta_g h}{|\mathbf v|_{\max}},
\qquad
\Delta t_{su} = \frac{h^2}{2\max_e \tau^{*}|\mathbf v|^2} = \frac{h}{|\mathbf v|_{\max}},
\qquad
\Delta t_{\text{rxn}} = \frac{1}{2\kappa_{\max}},
$$

the $\Delta t_{su}$ form following from $\tau^{*}|\mathbf v|^2 \le h|\mathbf v|/2$, and
combined as before:

$$
\frac{1}{\Delta t_{\text{comb}}} = \frac{1}{\Delta t_{\text{diff}}}
+ \frac{1}{\Delta t_{\text{adv}}} + \frac{1}{\Delta t_{su}},
\qquad
\Delta t = C_N \min\big(\Delta t_{\text{comb}},\ \Delta t_{\text{rxn}}\big).
$$

The scaling changed from $O(h^2/D_P)$ to $O(h/|\mathbf v|)$. Two consequences:

- **$|\mathbf v|$ is a runtime quantity.** The closed forms in `time_step.jl` bound the
  advective limit from `C_max`/`T_ref`, which no longer describes the operator. Δt is
  therefore recomputed each step from the Gauss-point velocities in `v_gp_cache`.
- **SU stabilization matters more, not less.** Once the pressure feedback is implicit the
  species step is a genuinely advection-dominated explicit scheme, and unstabilized
  Galerkin has no discrete maximum principle there (see `ADVECTION_STABILIZATION_NOTES.md`
  §"The mechanism"). Running IMPES with `advection_stabilization = 0` is permitted but
  warned about at startup.

### 5.1 Two controls the linear limits do not supply

Linear stability is necessary and not sufficient here. Both of the following were found by
running `advection_1d` under IMPES; without either, the run diverges.

**The clamp becomes a positive feedback.** In the explicit solver the `C_MIN` clamp is a
cosmetic floor: a species driven slightly negative is reset to zero, creating a little
mass, and nothing downstream notices within the step. Under IMPES the pressure equation
responds to that created mass *immediately* — higher total, higher pressure, higher
velocity, more overshoot, more created mass. On `advection_1d` this took four steps to go
from a 0.6 % artefact to $|\mathbf v| = 2\times10^4$ m/s. The step is therefore bounded by
a **positivity condition**

$$
\Delta t \le \tfrac12 \min_{i,\alpha \,:\, \dot C_{\alpha,i} < 0}
\frac{C_{\alpha,i}}{-\dot C_{\alpha,i}} ,
$$

evaluated from the previous step's rates and skipping nodes already essentially empty
(where a small negative rate is round-off and dividing by it would collapse the step for
nothing). This is not redundant with the Courant condition: it is what keeps the clamp
silent, and the clamp is what the feedback needs.

**The first step cannot be the asymptotic step.** A run begins on an initial condition the
implicit solve resolves in one step but the explicit species step cannot — `advection_1d`
starts with a 40 mol/m³ jump across a single element. Taking the full diffusion-limited
step (70 s, Fourier 0.9 — *linearly stable*) across that jump moves the minority species
by 60 mol/m³ from a value of 41, i.e. straight through zero and into the feedback above.
The solver therefore **starts at the explicit-equivalent step** — the harmonic limit
*including* the advective and SU terms, the step the explicit integrator would have used —
and lets the step grow by at most **1.1× per step**. The ramp is geometric, so it costs a
few dozen steps, and it resolves the startup transient at the resolution that transient
actually needs. 1.1 rather than something larger because the species step is first order
in time: on `advection_1d`, 1.2 left 0.25 % of avoidable error at the first output that
1.1 removes, for about twenty extra steps out of a hundred.

### 5.2 Measured outcome

`advection_1d`, 1 m column, h = 0.05 m, 2000 s:

| | Explicit | IMPES |
|---|---|---|
| Critical Δt (start) | 7.63e-2 s (advective) | 78.1 s (diffusive) |
| Steps taken | 29 136 | **136** |
| rel-L2 at t = 250 s | 1.707 % | 1.936 % |
| rel-L2 at t = 2000 s | 0.516 % | **0.307 %** |
| Total variation | 1.000× | 1.000× |

**214× fewer steps.** IMPES is slightly less accurate at the first output, where its step
is still ramping, and *more* accurate at every later output, where the explicit run's
error is dominated by its own spatial discretization rather than by time resolution.

## 6. Linear system and matrix-free PCG

$\mathbf S = \tfrac{1}{\Delta t}\mathbf M + \mathbf A^{n}$ is SPD, so a
Jacobi-preconditioned conjugate gradient needs no global matrix and no new dependency —
consistent with the rest of ADSIM, which forms no global system anywhere.

**Element operator**, rebuilt each step because $\Lambda^n$ depends on composition and
temperature:

$$
\mathbf A_e = \sum_{p=1}^{4} \Lambda^n(\xi_p)\,
\big(\nabla\mathbf N\big)_p \big(\nabla\mathbf N\big)_p^{\!\top} w_p W_p ,
\qquad
\Lambda^n(\xi_p) = R\,k\,T_{gp}\sum_\alpha \frac{C_{\alpha,gp}}{\mu_\alpha}
$$

stored as one $4\times4$ per element (128 B/element). $w_p$ carries $r_{gp}$ on an
axisymmetric run, so both formulations are handled without a special case — the same
guarantee `apply_element_stiffness!` relies on.

**Matvec** $\;\mathbf y = (\mathbf M/\Delta t)\odot\mathbf x + \sum_e \mathbf A_e \mathbf x_e$,
using the same fixed-chunk scatter as the flux loop.

**Preconditioner** $\;d_i = M_i/\Delta t + \sum_e \mathbf A_e[i_{loc},i_{loc}]$, accumulated
during assembly.

**Warm start** $\;x_0 = C_t^{n}$. The field moves little per step, so the iteration count
should be small and roughly constant; it will climb as Δt grows, which is expected.
Convergence is on relative residual with a hard iteration cap and a logged warning on
failure — never a silent bad answer.

## 7. Algorithm, one step

```
 1. refresh transient BC values at t_{n+1}
 2. assemble level-n explicit fluxes: q^d, q^T, q^g, q^su, q^b, q^r
 3. b^n = Sum_alpha(explicit terms) + M C_t^n / dt
    A^n  = per-Gauss-point element operators
    apply Dirichlet rows by symmetric elimination
 4. PCG  ->  C_t^{n+1}
 5. q^a_alpha(C_t^{n+1})  per species
 6. velocities from grad P^{n+1} = R(T^n grad C_t^{n+1} + C_t^{n+1} grad T^n) -> v_gp_cache
 7. Lagrange multipliers at absolute_pressure_bc nodes            (Sec. 4.5)
 8. species update + C_MIN clamp
 9. reaction, energy equation, BC re-imposition, P = R T C_t, probes, VTK
10. next dt from step 6's velocities and step 8's rates           (Sec. 5)
```

Note that step 6 precedes step 8, unlike the explicit solver, where the velocity is built
from $P^{k}$ *before* the concentration update so that the whole right-hand side sits at
level $k$. Under IMPES $P^{n+1}$ is available first and is the correct level to drive both
the species flux and the thermal advection.

## 8. Caveats and open points

**Lagged SU velocity.** $\tau^{*}_\alpha$ and $\mathbf a_\alpha$ depend on
$\mathbf v_\alpha \propto \nabla C_t$, but $q^{su}$ must enter $b^n$, which is formed
*before* $C_t^{n+1}$ exists. The SU velocity is therefore evaluated at $\nabla C_t^{n}$ in
both $b^n$ and (IMPES-S), which keeps §4.3 exact at the cost of a one-step lag in the
stabilization parameter only.

**Frozen mobility, and no Picard iteration.** $\Lambda^{n}$ is lagged, which linearizes
the one genuine nonlinearity of the pressure equation. The dependence is mild —
composition and temperature change slowly relative to pressure equilibration — and a
single pass is standard IMPES practice.

Relinearizing it is *not* simply a matter of repeating steps 3–5. $\Lambda$ is built from
$C_\alpha^{n}$, and the species field is not updated until after the pressure solve, so a
second pass over the current loop rebuilds an identical $\Lambda$ and changes nothing. A
real Picard iteration would have to carry a tentative species update — rates, throttle,
Lagrange correction and clamp — into a scratch field, rebuild $\Lambda$ from that, and
re-solve, committing only on the final pass. That is a restructuring of the step, not a
loop around it, and there is no evidence yet that it is needed. An
`impes_picard_iterations` key was implemented, found to be a no-op for exactly this
reason, and removed rather than shipped.

**Viscosity convention (pre-existing).** The transported flux uses the *species*
viscosity $\mu_\alpha$, so

$$
\Lambda = R k T \sum_\alpha \frac{C_\alpha}{\mu_\alpha} = \frac{R k T\, C_t}{\bar\mu_h},
\qquad \frac{1}{\bar\mu_h} = \sum_\alpha \frac{x_\alpha}{\mu_\alpha}
\quad\text{(harmonic mole-fraction mean)},
$$

whereas the reported/cached Darcy velocity uses the arithmetic mean
$\bar\mu_a = \sum_\alpha x_\alpha \mu_\alpha$. These differ; by AM–HM,
$\bar\mu_h \le \bar\mu_a$ always. For an equimolar CO₂/air mixture
($\mu = 1.47$ and $1.81\times10^{-5}$ Pa·s) the gap is ~1.1 %, so it is not currently
material — but it is a real inconsistency between the transported flux and the diagnostic
velocity, and it predates this work. The IMPES solver reproduces the existing convention
on both sides rather than silently changing the physics. Worth resolving separately.

**Reported time step.** With adaptive Δt, `--report-timestep` and the startup stability
report describe neither the initial nor the asymptotic step: they report the *limit* the
solution-independent mechanisms impose. The run starts below it (§5.1) and, once the
transient has passed, sits at whichever of the Courant, diffusion, reaction and positivity
limits binds. The advective and SU rows are retained in that report, marked implicit,
because they are the numbers that justify the change. Each output logs the step range
actually taken, the limit that set it, and the PCG iteration count.

**Volume-consistency floor.** The diagnostic of §4.3 cannot read below the accuracy of the
pressure solve: a relative residual of $10^{-10}$ leaves a solution error of order
$10^{-10}\kappa$, so a reading around $10^{-8}$ on a poorly conditioned mesh is the PCG
tolerance, not a formulation defect. On `advection_1d` the converged value is $10^{-11}$.
A reading that *grows* with time, or that jumps when the front moves, is the signal worth
chasing.

**Debug hook.** `ADSIM_IMPES_DEBUG=1` prints per-step Δt, max |v|, the Courant limit, the
largest species increment and the largest $|\Sigma C_\alpha - C_t|$ for the first dozen
steps. That is the startup transient, which is where every failure found so far began.

## References

- Aziz & Settari, *Petroleum Reservoir Simulation*, 1979 — the IMPES split and its
  volume-balance error in the general (multiphase) case.
- Brooks & Hughes, *Comput. Methods Appl. Mech. Engrg.* **32** (1982) 199–259 — SU/SUPG
  and the $\tau^{*}$ used here.
- `ADVECTION_STABILIZATION_NOTES.md` — why the species step needs stabilization at all.
- `Instructions_coding.md` — lumped mass, flow vectors, quadrature, pressure BC.
