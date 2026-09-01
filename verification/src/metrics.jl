#______________________________________________________
# metrics.jl
#
# Error measures for comparing a numerical profile against
# an analytical one. Pass/fail is gated on the relative L2
# norm; the others are reported for insight.
#______________________________________________________

module Metrics

export relative_l2, rmse, max_abs_error, total_variation, tv_ratio, relative_linf

"""
    relative_l2(num, exact) -> Float64

Relative discrete L2 error: ||num - exact||_2 / ||exact||_2.
Falls back to the absolute L2 norm when `exact` is essentially zero.
"""
function relative_l2(num::AbstractVector, exact::AbstractVector)
    length(num) == length(exact) || error("length mismatch: $(length(num)) vs $(length(exact))")
    num_norm = sqrt(sum(abs2, num .- exact))
    den = sqrt(sum(abs2, exact))
    return den > 1e-30 ? num_norm / den : num_norm
end

"""rmse(num, exact) — root-mean-square error."""
rmse(num::AbstractVector, exact::AbstractVector) = sqrt(sum(abs2, num .- exact) / length(num))

"""max_abs_error(num, exact) — worst-case pointwise |num - exact|."""
max_abs_error(num::AbstractVector, exact::AbstractVector) = maximum(abs.(num .- exact))

"""
    relative_linf(num, exact) -> Float64

Worst-case pointwise error relative to the scale of the exact solution.

L2 divides the error budget across the whole profile, so an error confined to a handful
of nodes barely moves it. L∞ does not, which is what makes it the right companion gate
for a defect that is localized rather than distributed.
"""
function relative_linf(num::AbstractVector, exact::AbstractVector)
    den = maximum(abs, exact)
    e = max_abs_error(num, exact)
    return den > 1e-30 ? e / den : e
end

"""
    total_variation(v) -> Float64

Discrete total variation Σ|v[i] - v[i-1]| along the profile.

A monotone profile has TV = |v[end] - v[1]| exactly; every node-to-node reversal adds
twice the size of the reversal. That makes TV the natural detector for an oscillation
riding on a smooth solution.
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

Total variation of the numerical profile divided by that of the exact one.

Reads as 1.0 for a solution that is as smooth as it should be and grows without bound as
node-to-node oscillation develops, independently of whether that oscillation is large
enough to move the L2 norm.

This is the metric that would have caught the axisymmetric instability. That run was
producing a bounded sawtooth against the axis of symmetry - bounded because the solver's
non-negativity clamp arrested it every step - which a relative-L2 gate averaged away
across the profile while the field was visibly wrong.

Note the direction it can see: the comparison is made along the profile the case
extracts, so it detects oscillation ALONG that line only. A case whose profile runs
axially cannot see a sawtooth in the radial direction; that is what the radial-profile
cases (`axisym_conduction_radial`, `axisym_darcy_radial`) are for.

# Flat reference profiles
`darcy_1d` compares a uniform Darcy velocity, whose exact total variation is identically
zero, so a bare ratio reports `Inf` for any discretisation error at all. The denominator
is therefore floored at 1% of the solution scale, which reads as: against a flat
reference, drift of up to about 1% of scale across the whole profile is not oscillation.
That case's measured variation is 0.26% of scale - a smooth monotone drift, not a
wiggle - and it scores 0.26 against a gate of 1.15.

The floor costs nothing in sensitivity. A sawtooth alternating by just 0.1% of scale
across 20 nodes has TV = 4% of scale and scores 4.0; one alternating by 1% scores 40.
On a profile with genuine variation the floor sits orders of magnitude below
`TV(exact)` and never binds.
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
