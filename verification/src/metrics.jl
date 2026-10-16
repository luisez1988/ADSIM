#______________________________________________________
# metrics.jl
#
# Error measures of a numerical profile against an analytical
# one. Pass/fail gated on the relative L2 norm; the other
# measures are informative.
#______________________________________________________

module Metrics

export relative_l2, rmse, max_abs_error, total_variation, tv_ratio, relative_linf

"""
    relative_l2(num, exact) -> Float64

Compute the relative discrete L2 error ||num - exact||_2 / ||exact||_2.
Absolute L2 norm when `exact` is ~0.
"""
function relative_l2(num::AbstractVector, exact::AbstractVector)
    length(num) == length(exact) || error("length mismatch: $(length(num)) vs $(length(exact))")
    num_norm = sqrt(sum(abs2, num .- exact))
    den = sqrt(sum(abs2, exact))
    return den > 1e-30 ? num_norm / den : num_norm
end

"""rmse(num, exact): root-mean-square error."""
rmse(num::AbstractVector, exact::AbstractVector) = sqrt(sum(abs2, num .- exact) / length(num))

"""max_abs_error(num, exact): maximum pointwise |num - exact|."""
max_abs_error(num::AbstractVector, exact::AbstractVector) = maximum(abs.(num .- exact))

"""
    relative_linf(num, exact) -> Float64

Compute the maximum pointwise error relative to the scale of the exact solution.

L2 averages the error over the whole profile; an error on a few nodes barely changes it.
L∞ complements it for localized defects.
"""
function relative_linf(num::AbstractVector, exact::AbstractVector)
    den = maximum(abs, exact)
    e = max_abs_error(num, exact)
    return den > 1e-30 ? e / den : e
end

"""
    total_variation(v) -> Float64

Compute the discrete total variation Σ|v[i] - v[i-1]| along the profile.

Monotone profile: TV = |v[end] - v[1]|. Each node-to-node reversal adds twice its size,
therefore TV detects oscillations on a smooth solution.
"""
function total_variation(v::AbstractVector)
    length(v) < 2 && return 0.0
    s = 0.0
    prev = first(v)
    for x in Iterators.drop(v, 1)
        s += abs(x - prev)
        prev = x
    end
    return s
end

"""
    tv_ratio(num, exact) -> Float64

Compute the total variation of the numerical profile over that of the exact one.

1.0 for a solution as smooth as the exact one; increases without bound with node-to-node
oscillation, also when the oscillation does not change the L2 norm.

Detects cases such as the axisymmetric instability: a bounded sawtooth at the symmetry
axis (bounded by the non-negativity clamp each step), averaged out by the relative-L2 gate.

Measured along the extracted profile only. An axial profile cannot detect a radial
sawtooth; radial-profile cases (`axisym_conduction_radial`, `axisym_darcy_radial`) cover it.

# Flat reference profiles
`darcy_1d` compares a uniform Darcy velocity (exact TV = 0), giving `Inf` for any
discretisation error. Denominator floored at 1% of the solution scale: against a flat
reference, drift up to ~1% of scale over the profile is not oscillation. Measured
variation of that case: 0.26% of scale (smooth monotone drift), score 0.26 vs. gate 1.15.

No loss of sensitivity: a sawtooth of 0.1% of scale over 20 nodes gives TV = 4% of scale
(score 4.0); 1% gives 40. For profiles with real variation the floor is orders of
magnitude below `TV(exact)`.
"""
function tv_ratio(num::AbstractVector, exact::AbstractVector)
    tv_e = total_variation(exact)
    tv_n = total_variation(num)
    scale = isempty(exact) ? 0.0 : maximum(abs, exact)

    den = max(tv_e, 1e-2 * scale)
    den > 1e-30 && return tv_n / den
    return tv_n > 1e-30 ? Inf : 1.0
end

end # module
