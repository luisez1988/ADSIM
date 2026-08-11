#______________________________________________________
# metrics.jl
#
# Error measures for comparing a numerical profile against
# an analytical one. Pass/fail is gated on the relative L2
# norm; the others are reported for insight.
#______________________________________________________

module Metrics

export relative_l2, rmse, max_abs_error

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

end # module
