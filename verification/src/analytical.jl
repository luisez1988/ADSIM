#______________________________________________________
# analytical.jl
#
# Closed-form reference solutions used to verify ADSIM.
# These mirror the formulas in the verification notebooks
# (Cleary & Adrian, 1973). erfc is implemented here so the
# suite has no external package dependency.
#______________________________________________________

module Analytical

export erfc, diffusion_series, advection_diffusion, evaluate

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
    elseif typ == "darcy_velocity"
        # Uniform in space and time; x and t are ignored.
        return darcy_velocity(; K=g("K"), mu=g("mu"), dC=g("dC"),
                              L=g("L", 1.0), T=g("T", 298.0), R=g("R", 8.314))
    else
        error("Unknown analytical type '$typ'. Add it to Analytical.evaluate.")
    end
end

end # module
