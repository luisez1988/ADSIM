#______________________________________________________
# analytical.jl
#
# Closed-form reference solutions for ADSIM verification.
# Same formulas as the verification notebooks (Cleary & Adrian, 1973).
# erfc implemented here: no external package dependency.
#______________________________________________________

module Analytical

export erfc, diffusion_series, advection_diffusion, conduction_cosine, linear_profile, carbonation_isothermal, carbonation_fdm, evaluate
export bessel_j0, bessel_j1, besselj0_zeros, conduction_disk

# ---- erf / erfc -----------------------------------------------------------
# Abramowitz & Stegun 7.1.26 rational approximation (|error| < 1.5e-7),
# extended to negative x by symmetry. Sufficient for a 2% L2 gate.
function _erf(x::Float64)
    s = sign(x)
    x = abs(x)
    t = 1.0 / (1.0 + 0.3275911 * x)
    y = 1.0 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t -
               0.284496736) * t + 0.254829592) * t * exp(-x * x)
    return s * y
end

"""erfc(x): complementary error function."""
erfc(x::Real) = 1.0 - _erf(Float64(x))

# ---- Pure diffusion (Fourier / consolidation series) ----------------------
"""
    diffusion_series(x, t; C0, Ci, D, R=1.0, L=1.0, terms=500)

Evaluate 1D diffusion into a column with a fixed-value boundary, sine series form
(same as `A_1D` in the diffusion notebook). `x`: distance from the driving boundary.
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

Evaluate 1D advection-diffusion in a finite column of length `L` with a first-type inlet
boundary (same as `A_1D(x,t,v,D,R,L)` in the advection notebook). `x`: distance from the inlet.
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

Return the uniform steady Darcy flux magnitude for a concentration difference `dC` over
length `L`. Ideal gas (P = C R T): pressure gradient `R T dC / L`, `v = (K/mu) R T dC / L`.
"""
darcy_velocity(; K, mu, dC, L, T=298.0, R=8.314) = (K / mu) * R * T * abs(dC) / L

# ---- Steady radial Darcy flow (axisymmetric) ------------------------------
"""
    darcy_radial(r; K, mu, C_w, C_e, r_w, r_e, T=298.0, R=8.314)

Return the steady radial Darcy velocity magnitude in an annulus `r_w <= r <= r_e` with
total concentrations `C_w` and `C_e` on the two cylindrical surfaces.

Position-dependent, unlike `darcy_velocity`. Only reference in the suite sensitive to a
missing axisymmetric radius factor in the advective kernels: in an axis-aligned column the
radius cancels between mass and stiffness, here it does not.

Steady state, ideal gas (`p = C R T`): same flux through every cylinder. With the
cylindrical divergence of Eq. (ax_divergence),

    d/dr ( r C dC/dr ) = 0    =>    r C dC/dr = A    (constant)

Integrating gives `C²` linear in `ln r`,

    C(r) = sqrt( C_w² + 2 A ln(r/r_w) ),   A = (C_e² - C_w²) / (2 ln(r_e/r_w))

Darcy's law then gives a velocity decreasing approximately as `1/r`:

    |v_r| = (K/mu) R T |A| / ( C(r) r )

Deviation from `1/r`: factor `C(r)`, varying by `|C_e - C_w| / C` across the annulus (a
fraction of a percent here). Kept exact, not linearised.

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
    # Clamp to the meshed annulus: nodes on the two radii, round-off must not move the
    # logarithm argument outside the profile range
    rc = clamp(r, r_w, r_e)
    C = sqrt(C_w^2 + 2 * A * log(rc / r_w))
    return (K / mu) * R * T * abs(A) / (C * rc)
end

"""
    radial_log_profile(r; C_w, C_e, r_w, r_e)

Return the steady total concentration between two cylinders at `C_w` and `C_e`, single
ideal gas driven by its own pressure gradient, no molecular diffusion:

    C(r) = sqrt( C_w² + 2 A ln(r/r_w) ),   A = (C_e² - C_w²) / (2 ln(r_e/r_w))

Integrated form of `r C dC/dr = constant` (same molar flux through every cylinder).
Preferred over the `darcy_radial` velocity: `C` is a primary nodal unknown, while the
nodal velocity is projected from the Gauss points and one-sided at the end nodes, with an
error unrelated to the axisymmetric weighting where `v` varies steeply.

Plane-strain profile: `C²` linear in `r` instead of `ln r`. Large shape difference over a
decade of radius, therefore a sharp discrimination between the two formulations.
"""
function radial_log_profile(r::Float64; C_w, C_e, r_w, r_e)
    A = (C_e^2 - C_w^2) / (2 * log(r_e / r_w))
    rc = clamp(r, r_w, r_e)
    return sqrt(C_w^2 + 2 * A * log(rc / r_w))
end

# ---- Heat conduction in an insulated bar ---------------------------------
"""
    conduction_cosine(x, t; T0, A, alpha, L=1.0)

Evaluate the temperature [K] of a 1-D bar insulated at both ends, initial single cosine
mode. cos(pi x / L) is the fundamental Laplacian eigenfunction under zero-flux
boundaries, therefore it decays without changing shape:

    T(x,t) = T0 + A cos(pi x / L) exp(-alpha pi^2 t / L^2)

`alpha = lambda_e / C_mix`: mixture thermal diffusivity. Zero flux is the natural
condition of the weak form, therefore unset thermal boundaries are insulated ends (no
thermal BCs required).

# Arguments (keyword)
- `T0`: mean temperature the bar relaxes to [K]
- `A`: initial amplitude of the mode [K]
- `alpha`: thermal diffusivity lambda_e / C_mix [m²/s]
- `L`: bar length [m]
"""
function conduction_cosine(x::Float64, t::Float64; T0, A, alpha, L=1.0)
    return T0 + A * cos(pi * x / L) * exp(-alpha * pi^2 * t / L^2)
end

# ---- Steady linear profile ------------------------------------------------
"""
    linear_profile(x; v0, vL, L=1.0)

Return the straight line from `v0` at x = 0 to `vL` at x = L, time independent.

Steady pure conduction between two prescribed temperatures without source (zero second
derivative). Exact instead of a truncated series: the discrete operator must give zero
residual on a linear field, exposing errors in conductivity, gradient operator or
Dirichlet gate.
"""
linear_profile(x::Float64; v0, vL, L=1.0) = v0 + (vL - v0) * x / L

# ---- Isothermal carbonation (0-D, constant CO2) ---------------------------
#
# Two references for the same rate law, sharing the coefficient helpers below:
#   carbonation_isothermal: closed form, exact for a constant coefficient
#   carbonation_fdm:        explicit finite-difference integration of the ODE
# FDM is the general one (no constant-coefficient requirement). With a constant
# coefficient both agree to round-off, a consistency check of the two codings.

"""
    _henry_solubility(T)

Henry solubility, van 't Hoff form. Same as `henry_solubility` in the solver.
"""
function _henry_solubility(T)
    K_H_ref = 3.3e-4
    T_ref = 298.15
    return K_H_ref * exp(2400.0 * (1.0 / T - 1.0 / T_ref))
end

"""
    _area_factor(C_total, T, beta; R=8.3145)

Interfacial-area factor, pressure-ratio form. Same as `interfacial_area_factor` in the
solver, floor included:

    a = exp[b (max{p/P_atm - 1, 0})],   p = C_total R T

Driver: TOTAL gas pressure over all species (pressure on the menisci; the calibration
chamber pressure is a total). Equal to the CO2 partial pressure for single-gas CO2,
hence the `C_total = C_co2` default of `carbonation_fdm`.

`beta` DIMENSIONLESS. The earlier dissolved-concentration form had units of m3/mol and is
not interchangeable (see solver docstring).
"""
function _area_factor(C_total, T, beta; R=8.3145)
    beta == 0.0 && return 1.0
    P_atm = 101325.0
    p = C_total * R * T
    return exp(beta * max(p / P_atm - 1.0, 0.0))
end

"""
    _carbonation_lambda(; k_o, E, beta, T, theta_w, C_co2, C_total, R=8.3145)

Decay constant of the lime ODE, assembled in the same way as the solver:

    lambda = k_T * a(C_total) * theta_w * K_H(T) * R * T * C_co2

Solver: `r = k_T a C_aq (A_s - A_r)`, `C_aq = K_H R T C_g,CO2` per m3 of water, then
`dA_s/dt = -theta_w r` per m3 of total volume.

`C_aq` (dissolved CO2) follows the CO2 partial pressure; `a` (gas-liquid interface)
follows the total gas pressure.
"""
function _carbonation_lambda(; k_o, E, beta, T, theta_w, C_co2, C_total, R=8.3145)
    k_T = k_o * exp(-E / (R * T))
    a = _area_factor(C_total, T, beta; R=R)
    return k_T * a * theta_w * _henry_solubility(T) * R * T * C_co2
end

"""
    carbonation_isothermal(t; k_o, E, beta, T, theta_w, C_co2, p_r, R=8.3145)

Return the degree of carbonation at constant temperature and CO2 concentration,
Eq. (5.4) of `calibration_reactionrate/integrator_mathematics.md`:

    DoC(t) = (1 - p_r) (1 - exp(-lambda t))

Exact closed form for constant lambda ("exact to round-off at any step count" in the
reference document). lambda constant here: CO2 prescribed by a Dirichlet boundary on
every node, constant temperature (zero reaction enthalpy, no heat released).

# Arguments (keyword, same as the material file)
- `k_o`: Arrhenius factor [m3 mol^-1 s^-1]
- `E`: activation energy [J/mol]
- `beta`: interfacial-area coefficient [-], dimensionless; 0 disables the factor
- `T`: absolute temperature [K]
- `theta_w`: volumetric water content, porosity x saturation [-]
- `C_co2`: prescribed CO2 gas concentration [mol per m3 of gas]; drives C_aq
- `C_total`: prescribed TOTAL gas concentration [mol per m3 of gas]; drives the
  interfacial-area factor. Defaults to `C_co2` (single-gas case)
- `p_r`: residual lime as a fraction of the initial lime [-]

# Returns
- Degree of carbonation [-], approaching `1 - p_r` instead of 1
"""
function carbonation_isothermal(t::Float64; k_o, E, beta, T, theta_w, C_co2, p_r,
                                C_total=nothing, R=8.3145)
    t <= 0 && return 0.0
    C_tot = C_total === nothing ? C_co2 : C_total
    lam = _carbonation_lambda(; k_o=k_o, E=E, beta=beta, T=T, theta_w=theta_w,
                                C_co2=C_co2, C_total=C_tot, R=R)
    return (1.0 - p_r) * (1.0 - exp(-lam * t))
end

"""
    carbonation_fdm(t; k_o, E, beta, T, theta_w, C_co2, p_r, R=8.3145, substeps=20000)

Integrate the degree of carbonation of the same rate law numerically instead of the
closed form.

Normalised lime ODE, u = A_s / A_s0, u(0) = 1:

    du/dt = -lambda(u, t) (u - p_r),      DoC = 1 - u

Classical RK4, `substeps` uniform steps from 0 to `t`. lambda rebuilt each stage from
`_carbonation_lambda` (single definition of the law, no constant-lambda assumption).
Agrees to round-off with `carbonation_isothermal` on constant-coefficient cases; remains
valid when the CO2 field or temperature changes.

RK4 instead of forward Euler: a first-order reference error would mask the ADSIM
forward-Euler lime update error. Integration error at the default `substeps` far below the
0.5 % gate.

# Arguments
Same as [`carbonation_isothermal`](@ref), including the `C_total` default, plus:
- `substeps`: RK4 steps to reach `t` [-]

# Returns
- Degree of carbonation [-], approaching `1 - p_r` instead of 1
"""
function carbonation_fdm(t::Float64; k_o, E, beta, T, theta_w, C_co2, p_r,
                         C_total=nothing, R=8.3145, substeps::Int=20000)
    t <= 0 && return 0.0

    # Single-gas default: total = CO2. Mixed-gas cases must pass C_total
    C_tot = C_total === nothing ? C_co2 : C_total

    # Closure: coefficient rebuilt at every stage, no constant-coefficient assumption
    du = (u, tau) -> -_carbonation_lambda(; k_o=k_o, E=E, beta=beta, T=T,
                                            theta_w=theta_w, C_co2=C_co2,
                                            C_total=C_tot, R=R) * (u - p_r)

    h = t / substeps
    u = 1.0
    tau = 0.0
    for _ in 1:substeps
        k1 = du(u, tau)
        k2 = du(u + 0.5 * h * k1, tau + 0.5 * h)
        k3 = du(u + 0.5 * h * k2, tau + 0.5 * h)
        k4 = du(u + h * k3, tau + h)
        u += (h / 6.0) * (k1 + 2k2 + 2k3 + k4)
        tau += h
    end

    return 1.0 - u
end

# ---- Bessel functions J0, J1 and the zeros of J0 --------------------------
# Abramowitz & Stegun 9.4.1-9.4.6 polynomial/rational approximations (no external
# package, as for erfc). |error| < 5e-8 for J0 and < 1.3e-8 for J1, orders of
# magnitude below the gated discretisation errors.

"""bessel_j0(x): Bessel function of the first kind, order 0."""
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

"""bessel_j1(x): Bessel function of the first kind, order 1."""
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

Return the first `n` positive zeros of J0. McMahon asymptotic expansion (A&S 9.5.12)
as the initial guess, refined by Newton's method with J0' = -J1 (two or three
iterations, ~1e-8 accuracy, the floor of the J0 approximation).
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

Evaluate the temperature [K] of a solid disk (or infinite cylinder) of radius `R`,
initially uniform at `T0`, rim at `Ts` from t = 0. Separation of variables in polar
coordinates gives the Fourier-Bessel series

    (T - Ts)/(T0 - Ts) = sum_n  2/(j_n J_1(j_n))  J_0(j_n r/R)  exp(-alpha j_n^2 t/R^2)

`j_n`: n-th positive zero of J0 (Carslaw & Jaeger, §7.6). Fundamental mode time
constant: `R^2/(alpha j_1^2)`.

Only two-dimensional reference of the suite. Other cases use axis-aligned strips
(diagonal isoparametric Jacobian, mapping errors cancel); here the elements are skewed
general quads and the solution must remain axisymmetric on a non-axisymmetric mesh.

# Arguments (keyword)
- `T0`: uniform initial temperature [K]
- `Ts`: prescribed rim temperature [K]
- `alpha`: thermal diffusivity lambda_e / C_mix [m²/s]
- `R`: disk radius [m]
- `terms`: number of series terms; the n-th decays as exp(-j_n^2 alpha t/R^2),
  default well above the requirement of the earliest compared snapshot
"""
function conduction_disk(r::Float64, t::Float64; T0, Ts, alpha, R=1.0, terms::Int=200)
    t <= 0 && return T0
    j = besselj0_zeros(terms)
    s = 0.0
    @inbounds for n in 1:terms
        jn = j[n]
        # Coefficients decay only as jn^-1/2; convergence from the exponential.
        # Stop once jn^2 alpha t/R^2 is large.
        e = exp(-alpha * jn^2 * t / R^2)
        e < 1e-16 && break
        s += 2.0 * bessel_j0(jn * r / R) / (jn * bessel_j1(jn)) * e
    end
    return Ts + (T0 - Ts) * s
end

"""
    evaluate(spec::Dict, x, t) -> Float64

Evaluate the analytical solution selected by `spec["type"]` (case `.toml`) at position
`x` and time `t`. New analytical models: add a branch here.
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
                                      C_co2=g("C_co2"), C_total=g("C_total", g("C_co2")),
                                      p_r=g("p_r", 0.0), R=g("R", 8.3145))
    elseif typ == "carbonation_fdm"
        # Same law as carbonation_isothermal, numerically integrated
        return carbonation_fdm(t; k_o=g("k_o"), E=g("E", 0.0), beta=g("beta", 0.0),
                               T=g("T", 298.15), theta_w=g("theta_w"),
                               C_co2=g("C_co2"), C_total=g("C_total", g("C_co2")),
                               p_r=g("p_r", 0.0), R=g("R", 8.3145),
                               substeps=Int(g("substeps", 20000)))
    else
        error("Unknown analytical type '$typ'. Add it to Analytical.evaluate.")
    end
end

end # module
