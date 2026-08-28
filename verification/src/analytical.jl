#______________________________________________________
# analytical.jl
#
# Closed-form reference solutions used to verify ADSIM.
# These mirror the formulas in the verification notebooks
# (Cleary & Adrian, 1973). erfc is implemented here so the
# suite has no external package dependency.
#______________________________________________________

module Analytical

export erfc, diffusion_series, advection_diffusion, conduction_cosine, linear_profile, carbonation_isothermal, evaluate
export bessel_j0, bessel_j1, besselj0_zeros, conduction_disk

# ---- erf / erfc -----------------------------------------------------------
# Abramowitz & Stegun 7.1.26 rational approximation (|error| < 1.5e-7),
# extended to negative x by symmetry. Plenty accurate for a 2% L2 gate.
function _erf(x::Float64)
    s = sign(x)
    x = abs(x)
    t = 1.0 / (1.0 + 0.3275911 * x)
    y = 1.0 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t -
               0.284496736) * t + 0.254829592) * t * exp(-x * x)
    return s * y
end

"""erfc(x) — complementary error function."""
erfc(x::Real) = 1.0 - _erf(Float64(x))

# ---- Pure diffusion (Fourier / consolidation series) ----------------------
"""
    diffusion_series(x, t; C0, Ci, D, R=1.0, L=1.0, terms=500)

1D diffusion into a semi-infinite/finite column with a fixed-value boundary,
written as the classic sine series (matches `A_1D` in the diffusion notebook).
`x` is distance from the driving boundary.
"""
function diffusion_series(x::Float64, t::Float64; C0, Ci, D, R=1.0, L=1.0, terms::Int=500)
    t <= 0 && return Ci
    A = 0.0
    @inbounds for m in 0:terms-1
        M = (2m + 1) * pi / 2
        A += 2 * sin(M * x / L) * exp(-(M^2) * D * R * t / (L^2)) / M
    end
    return Ci + (C0 - Ci) * (1 - A)
end

# ---- Diffusion + advection (Cleary & Adrian, 1973) ------------------------
"""
    advection_diffusion(x, t; C0, Ci, D, v, R=1.0, L=1.0)

1D advection-diffusion with a first-type inlet boundary and a finite column
of length `L`. Mirrors `A_1D(x,t,v,D,R,L)` in the advection notebook.
`x` is distance from the inlet boundary.
"""
function advection_diffusion(x::Float64, t::Float64; C0, Ci, D, v, R=1.0, L=1.0)
    t <= 0 && return Ci
    den = 2 * sqrt(D * R * t)
    term1 = 0.5 * erfc((R * x - v * t) / den)
    term2 = 0.5 * exp(v * x / D) * erfc((R * x + v * t) / den)
    term3 = 0.5 * (2 + v * (2L - x) / D + v^2 * t / (D * R)) *
            exp(v * L / D) * erfc((R * (2L - x) + v * t) / den)
    term4 = -sqrt(v^2 * t / (pi * D * R)) *
            exp(v * L / D - R / (4 * D * t) * (2L - x + v * t / R)^2)
    A = term1 + term2 + term3 + term4
    return Ci + (C0 - Ci) * A
end

# ---- Darcy velocity (steady, uniform) -------------------------------------
"""
    darcy_velocity(; K, mu, dC, L, T=298.0, R=8.314)

Steady Darcy flux driven by a fixed concentration difference across a column.
The gas is ideal (P = C R T), so a difference `dC` over length `L` sets a
pressure gradient `R T dC / L`, and Darcy's law gives `v = (K/mu) R T dC / L`.
Independent of position and time (returns the uniform magnitude).
"""
darcy_velocity(; K, mu, dC, L, T=298.0, R=8.314) = (K / mu) * R * T * abs(dC) / L

# ---- Steady radial Darcy flow (axisymmetric) ------------------------------
"""
    darcy_radial(r; K, mu, C_w, C_e, r_w, r_e, T=298.0, R=8.314)

Magnitude of the steady radial Darcy velocity in an annulus `r_w <= r <= r_e` whose two
cylindrical surfaces are held at total concentrations `C_w` and `C_e`.

Unlike the plane `darcy_velocity` this varies with position, and that is the point: it is
the only reference in the suite that fails outright if the axisymmetric radius factor is
missing from the advective kernels. In a column aligned with the axis the radius cancels
between mass and stiffness; here it does not, so a missing `r_gp` changes the answer instead
of dividing out.

Steady state with an ideal gas (`p = C R T`) requires the same flux through every cylinder,
so with the cylindrical divergence of Eq. (ax_divergence),

    d/dr ( r C dC/dr ) = 0    =>    r C dC/dr = A    (constant)

Integrating gives `C²` linear in `ln r`,

    C(r) = sqrt( C_w² + 2 A ln(r/r_w) ),   A = (C_e² - C_w²) / (2 ln(r_e/r_w))

and Darcy's law then gives a velocity falling off essentially as `1/r`:

    |v_r| = (K/mu) R T |A| / ( C(r) r )

The departure from an exact `1/r` law is the factor `C(r)`, which varies by
`|C_e - C_w| / C` across the annulus - a fraction of a percent for the concentrations used
here, but carried exactly rather than linearised.

# Arguments (keyword)
- `K`: intrinsic permeability [m²]
- `mu`: gas dynamic viscosity [Pa.s]
- `C_w`, `C_e`: prescribed total concentration at the inner and outer radius [mol/m³]
- `r_w`, `r_e`: inner and outer radius [m]
- `T`: temperature [K]
- `R`: universal gas constant [J/(mol.K)]

# Returns
- Magnitude of the radial Darcy velocity at radius `r` [m/s]
"""
function darcy_radial(r::Float64; K, mu, C_w, C_e, r_w, r_e, T=298.0, R=8.314)
    A = (C_e^2 - C_w^2) / (2 * log(r_e / r_w))
    # Clamp to the meshed annulus: nodes sit exactly on the two radii, so round-off in a
    # coordinate must not push the logarithm outside the range the profile is defined on.
    rc = clamp(r, r_w, r_e)
    C = sqrt(C_w^2 + 2 * A * log(rc / r_w))
    return (K / mu) * R * T * abs(A) / (C * rc)
end

"""
    radial_log_profile(r; C_w, C_e, r_w, r_e)

Steady total concentration between two cylinders held at `C_w` and `C_e`, for a single ideal
gas driven by its own pressure gradient with no molecular diffusion:

    C(r) = sqrt( C_w² + 2 A ln(r/r_w) ),   A = (C_e² - C_w²) / (2 ln(r_e/r_w))

This is the integrated form of `r C dC/dr = constant`, the statement that the same molar flux
crosses every cylinder. It is scored in preference to the velocity of `darcy_radial` because
`C` is a primary nodal unknown, whereas the nodal velocity is projected from the Gauss points
and is therefore one-sided at the two end nodes - where `v` varies steeply, that projection
carries an error of its own that has nothing to do with the axisymmetric weighting.

The corresponding plane-strain profile has `C²` linear in `r` rather than in `ln r`, which is
a large difference in shape over a decade of radius, so this discriminates sharply between the
two formulations rather than merely tolerating both.
"""
function radial_log_profile(r::Float64; C_w, C_e, r_w, r_e)
    A = (C_e^2 - C_w^2) / (2 * log(r_e / r_w))
    rc = clamp(r, r_w, r_e)
    return sqrt(C_w^2 + 2 * A * log(rc / r_w))
end

# ---- Heat conduction in an insulated bar ---------------------------------
"""
    conduction_cosine(x, t; T0, A, alpha, L=1.0)

Temperature in a 1-D bar insulated at both ends, starting from a single cosine
mode. cos(pi x / L) is the fundamental eigenfunction of the Laplacian under
zero-flux boundaries, so it decays without changing shape:

    T(x,t) = T0 + A cos(pi x / L) exp(-alpha pi^2 t / L^2)

with `alpha = lambda_e / C_mix` the thermal diffusivity of the mixture. Zero flux
is the natural boundary condition of the weak form, so an unset thermal boundary
is exactly an insulated end — no boundary condition machinery is required, which
is what lets this case run before thermal boundary conditions exist.

# Arguments (keyword)
- `T0`: mean temperature the bar relaxes to [K]
- `A`: initial amplitude of the mode [K]
- `alpha`: thermal diffusivity lambda_e / C_mix [m²/s]
- `L`: bar length [m]

# Returns
- Temperature at position `x` and time `t` [K]
"""
function conduction_cosine(x::Float64, t::Float64; T0, A, alpha, L=1.0)
    return T0 + A * cos(pi * x / L) * exp(-alpha * pi^2 * t / L^2)
end

# ---- Steady linear profile ------------------------------------------------
"""
    linear_profile(x; v0, vL, L=1.0)

Straight line from `v0` at x = 0 to `vL` at x = L, independent of time.

This is the steady state of pure conduction between two prescribed temperatures with
no source: the second derivative vanishes, so the profile is exactly linear. Being
exact rather than a truncated series, it is a strict check — the discrete operator
has to return zero residual on a linear field, so any error in the conductivity, the
gradient operator or the Dirichlet gate is visible immediately.
"""
linear_profile(x::Float64; v0, vL, L=1.0) = v0 + (vL - v0) * x / L

# ---- Isothermal carbonation (0-D, constant CO2) ---------------------------
"""
    carbonation_isothermal(t; k_o, E, beta, T, theta_w, C_co2, p_r, R=8.3145)

Degree of carbonation of a specimen held at constant temperature and constant
CO2 concentration. This is Eq. (5.4) of
`calibration_reactionrate/integrator_mathematics.md`:

    DoC(t) = (1 - p_r) (1 - exp(-λ t))

Whenever λ is constant the rate law integrates in closed form, and the result is
exact — the reference document records it as "exact to round-off at any step
count". λ is constant here because the CO2 field is pinned by a Dirichlet
boundary on every node and the temperature does not move (the case sets the
reaction enthalpy to zero, so the reaction releases no heat).

λ is assembled exactly as the solver assembles it, so the two are the same
number rather than two expressions that happen to agree:

    λ = k_T · a · θ_w · K_H(T) · R · T · C_co2

The solver computes `r = k_T a C_aq (A_s - A_r)` with `C_aq = K_H R T C_g` per
m³ of water, then `dA_s/dt = -θ_w r` per m³ of total volume. The notebook writes
`λ = k_T θ_w K_H p_CO2`; since `p_CO2 = C_co2 R T` for an ideal gas the two are
identical.

# Arguments (keyword, all mirroring the material file)
- `k_o`: Arrhenius factor [m³ mol⁻¹ s⁻¹]
- `E`: activation energy [J/mol]
- `beta`: interfacial-area coefficient [m³/mol]; 0 disables the area factor
- `T`: absolute temperature [K]
- `theta_w`: volumetric water content, porosity × saturation [-]
- `C_co2`: prescribed CO2 gas concentration [mol per m³ of gas]
- `p_r`: residual lime as a fraction of the initial lime [-]

# Returns
- Degree of carbonation [-], approaching `1 - p_r` rather than 1
"""
function carbonation_isothermal(t::Float64; k_o, E, beta, T, theta_w, C_co2, p_r, R=8.3145)
    t <= 0 && return 0.0

    # Henry solubility, van 't Hoff form (solver: henry_solubility)
    K_H_ref = 3.3e-4
    T_ref = 298.15
    K_H = K_H_ref * exp(2400.0 * (1.0 / T - 1.0 / T_ref))

    # Interfacial-area factor (solver: interfacial_area_factor). Evaluated at the
    # reference temperature by construction, and constant here because C_co2 is.
    a = 1.0
    if beta != 0.0
        P_atm = 101325.0
        x_co2_atm = 420e-6
        K_H_at_ref = K_H_ref
        C_aq_ref = K_H_at_ref * R * T_ref * C_co2
        C_aq_atm = K_H_at_ref * x_co2_atm * P_atm
        a = exp(beta * (C_aq_ref - C_aq_atm))
    end

    k_T = k_o * exp(-E / (R * T))
    λ = k_T * a * theta_w * K_H * R * T * C_co2

    return (1.0 - p_r) * (1.0 - exp(-λ * t))
end

# ---- Bessel functions J0, J1 and the zeros of J0 --------------------------
# Abramowitz & Stegun 9.4.1-9.4.6 polynomial/rational approximations, kept here
# for the same reason as erfc above: the suite must not need an external package.
# |error| < 5e-8 for J0 and < 1.3e-8 for J1, orders of magnitude below the
# discretisation error any case gates on.

"""bessel_j0(x) — Bessel function of the first kind, order 0."""
function bessel_j0(x::Float64)
    ax = abs(x)
    if ax < 3.0
        y = (x / 3.0)^2
        return 1.0 + y * (-2.2499997 + y * (1.2656208 + y * (-0.3163866 +
               y * (0.0444479 + y * (-0.0039444 + y * 0.0002100)))))
    end
    z = 3.0 / ax
    f = 0.79788456 + z * (-0.00000077 + z * (-0.00552740 + z * (-0.00009512 +
        z * (0.00137237 + z * (-0.00072805 + z * 0.00014476)))))
    θ = ax - 0.78539816 + z * (-0.04166397 + z * (-0.00003954 + z * (0.00262573 +
        z * (-0.00054125 + z * (-0.00029333 + z * 0.00013558)))))
    return f * cos(θ) / sqrt(ax)
end

"""bessel_j1(x) — Bessel function of the first kind, order 1."""
function bessel_j1(x::Float64)
    ax = abs(x)
    if ax < 3.0
        y = (x / 3.0)^2
        r = x * (0.5 + y * (-0.56249985 + y * (0.21093573 + y * (-0.03954289 +
            y * (0.00443319 + y * (-0.00031761 + y * 0.00001109))))))
        return r
    end
    z = 3.0 / ax
    f = 0.79788456 + z * (0.00000156 + z * (0.01659667 + z * (0.00017105 +
        z * (-0.00249511 + z * (0.00113653 + z * (-0.00020033))))))
    θ = ax - 2.35619449 + z * (0.12499612 + z * (0.00005650 + z * (-0.00637879 +
        z * (0.00074348 + z * (0.00079824 + z * (-0.00029166))))))
    r = f * cos(θ) / sqrt(ax)
    return x < 0 ? -r : r
end

"""
    besselj0_zeros(n) -> Vector{Float64}

First `n` positive zeros of J0. Started from McMahon's asymptotic expansion
(A&S 9.5.12) and polished with Newton's method on J0' = -J1, which converges in
two or three iterations and lands each root to about 1e-8 — the accuracy floor
of the J0 approximation itself.
"""
function besselj0_zeros(n::Int)
    j = zeros(Float64, n)
    for s in 1:n
        β = (s - 0.25) * pi
        x = β + 1.0 / (8β) - 31.0 / (384 * β^3)
        for _ in 1:8
            dx = bessel_j0(x) / bessel_j1(x)   # -J0/J0' with J0' = -J1
            x += dx
            abs(dx) < 1e-13 && break
        end
        j[s] = x
    end
    return j
end

# ---- Transient conduction in a solid disk ---------------------------------
"""
    conduction_disk(r, t; T0, Ts, alpha, R=1.0, terms=200)

Temperature in a solid disk (or infinite cylinder) of radius `R`, initially
uniform at `T0`, with its whole rim held at `Ts` from t = 0 on. Separation of
variables in polar coordinates gives the Fourier-Bessel series

    (T - Ts)/(T0 - Ts) = sum_n  2/(j_n J_1(j_n))  J_0(j_n r/R)  exp(-alpha j_n^2 t/R^2)

with `j_n` the n-th positive zero of J0 (Carslaw & Jaeger, §7.6). The fundamental
mode decays with time constant `R^2/(alpha j_1^2)`.

This is the suite's only genuinely two-dimensional reference. Every other case
runs on an axis-aligned strip, where the isoparametric Jacobian is diagonal and
a whole class of mapping errors cancels; here the elements are skewed general
quads, and the solution has to stay axisymmetric on a mesh that is not.

# Arguments (keyword)
- `T0`: uniform initial temperature [K]
- `Ts`: prescribed rim temperature [K]
- `alpha`: thermal diffusivity lambda_e / C_mix [m²/s]
- `R`: disk radius [m]
- `terms`: number of series terms; the n-th decays as exp(-j_n^2 alpha t/R^2),
  so the default is far more than the earliest compared snapshot needs

# Returns
- Temperature at radius `r` and time `t` [K]
"""
function conduction_disk(r::Float64, t::Float64; T0, Ts, alpha, R=1.0, terms::Int=200)
    t <= 0 && return T0
    j = besselj0_zeros(terms)
    s = 0.0
    @inbounds for n in 1:terms
        jn = j[n]
        # Coefficients fall off only as jn^-1/2, so it is the exponential that
        # converges the series: past the point where jn^2 alpha t/R^2 is large the
        # remaining terms cannot matter, and the loop can stop.
        e = exp(-alpha * jn^2 * t / R^2)
        e < 1e-16 && break
        s += 2.0 * bessel_j0(jn * r / R) / (jn * bessel_j1(jn)) * e
    end
    return Ts + (T0 - Ts) * s
end

"""
    evaluate(spec::Dict, x, t) -> Float64

Dispatch on `spec["type"]` (a value read from a case `.toml`) and evaluate the
matching analytical solution at position `x` and time `t`. Add a branch here to
register a new analytical model for a new verification case.
"""
function evaluate(spec::AbstractDict, x::Float64, t::Float64)
    g(k, d=nothing) = get(spec, k, d)
    typ = spec["type"]
    if typ == "diffusion_series"
        return diffusion_series(x, t; C0=g("C0"), Ci=g("Ci"), D=g("D"),
                                R=g("R", 1.0), L=g("L", 1.0), terms=Int(g("terms", 500)))
    elseif typ == "advection_diffusion"
        return advection_diffusion(x, t; C0=g("C0"), Ci=g("Ci"), D=g("D"),
                                   v=g("v"), R=g("R", 1.0), L=g("L", 1.0))
    elseif typ == "radial_log_profile"
        # x is the nodal radius, supplied by line_profile with axis = :x and origin = 0.
        return radial_log_profile(x; C_w=g("C_w"), C_e=g("C_e"), r_w=g("r_w"), r_e=g("r_e"))
    elseif typ == "darcy_radial"
        # x is the nodal radius, supplied by line_profile with axis = :x and origin = 0.
        return darcy_radial(x; K=g("K"), mu=g("mu"), C_w=g("C_w"), C_e=g("C_e"),
                            r_w=g("r_w"), r_e=g("r_e"), T=g("T", 298.0), R=g("R", 8.314))
    elseif typ == "darcy_velocity"
        # Uniform in space and time; x and t are ignored.
        return darcy_velocity(; K=g("K"), mu=g("mu"), dC=g("dC"),
                              L=g("L", 1.0), T=g("T", 298.0), R=g("R", 8.314))
    elseif typ == "conduction_cosine"
        return conduction_cosine(x, t; T0=g("T0"), A=g("A"), alpha=g("alpha"), L=g("L", 1.0))
    elseif typ == "conduction_disk"
        # x is the nodal radius, supplied by line_profile with axis = :radius.
        return conduction_disk(x, t; T0=g("T0"), Ts=g("Ts"), alpha=g("alpha"),
                               R=g("R", 1.0), terms=Int(g("terms", 200)))
    elseif typ == "linear_profile"
        # Steady state: time is ignored.
        return linear_profile(x; v0=g("v0"), vL=g("vL"), L=g("L", 1.0))
    elseif typ == "carbonation_isothermal"
        # Uniform in space; x is ignored, t is the probe sample time.
        return carbonation_isothermal(t; k_o=g("k_o"), E=g("E", 0.0), beta=g("beta", 0.0),
                                      T=g("T", 298.15), theta_w=g("theta_w"),
                                      C_co2=g("C_co2"), p_r=g("p_r", 0.0),
                                      R=g("R", 8.3145))
    else
        error("Unknown analytical type '$typ'. Add it to Analytical.evaluate.")
    end
end

end # module
