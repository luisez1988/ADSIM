"""
    shape_functions.jl

Module for isoparametric shape functions and Jacobian calculations for 4-noded quad elements.
Uses 2x2 Gaussian quadrature (4 Gauss points) for integration.
"""

module ShapeFunctions

using Base.Threads

export ShapeFunctionData, initialize_shape_functions!, get_N, get_B, get_invJ, get_detJ, get_dN_dx,
       get_weight, is_axisymmetric

"""
    ShapeFunctionData

Structure to store precomputed shape function values and element-specific Jacobian data.

# Fields
- `N::Vector{Vector{Float64}}`: Shape functions at each Gauss point [4 Gauss points][4 nodes]
- `B::Vector{Matrix{Float64}}`: Derivatives of shape functions at each Gauss point [4 Gauss points][4 nodes, 2 coords]
- `invJ::Array{Float64, 4}`: Inverse Jacobian at each element and Gauss point [Nelements, 4 Gauss points, 2, 2]
- `detJ::Matrix{Float64}`: Determinant of Jacobian at each element and Gauss point [Nelements, 4 Gauss points]
- `weight::Matrix{Float64}`: Quadrature measure at each element and Gauss point [Nelements, 4 Gauss points]
- `dN_dx::Array{Float64, 4}`: Shape function derivatives in physical coordinates, B · J^-1,
  at each element and Gauss point [Nelements, 4 Gauss points, 4 nodes, 2 coords]
- `axisymmetric::Bool`: True when `weight` carries the radius factor of Appendix `app:axisymmetric`
- `gauss_points::Matrix{Float64}`: Gauss point coordinates in isoparametric space [4 points, 2 coords (ξ, η)]
- `gauss_weights::Vector{Float64}`: Gauss point weights [4 points]

`dN_dx` is the quantity every flux term actually consumes. It is stored rather than recomputed
from `B` and `invJ` at each use, because the product is element- and Gauss-point-specific but
time-invariant, so recomputing it inside the time loop is pure waste.

`weight` is what carries the plane/axisymmetric distinction, and it is the ONLY thing that does.
The axisymmetric appendix reduces to a single substitution: wherever the plane formulation
multiplies an integrand by the area measure `|det J| W_gp`, the axisymmetric one multiplies it by
`r_gp |det J| W_gp`, leaving every gradient and stiffness contraction untouched. Folding `r_gp`
into a cached measure at initialization therefore makes every quadrature sum in the solver correct
in both formulations with no branch, no second code path and no extra arithmetic per time step:

    weight[e, p] = axisymmetric ? r_gp * |det J| : |det J|

`detJ` is kept separately and still means the Jacobian determinant. The split is deliberate and is
what makes the change auditable: a site reading `get_weight` is a quadrature sum and must carry the
radius, a site reading `get_detJ` is geometry and must not.
"""
mutable struct ShapeFunctionData
    N::Vector{Vector{Float64}}
    B::Vector{Matrix{Float64}}
    invJ::Array{Float64, 4}
    detJ::Matrix{Float64}
    weight::Matrix{Float64}
    dN_dx::Array{Float64, 4}
    axisymmetric::Bool
    gauss_points::Matrix{Float64}
    gauss_weights::Vector{Float64}
end

# Global storage for shape function data
global shape_funcs::Union{ShapeFunctionData, Nothing} = nothing

"""
    compute_shape_functions(ξ::Float64, η::Float64)

Compute the shape functions for a 4-noded quad element at isoparametric coordinates (ξ, η).

# Arguments
- `ξ::Float64`: Isoparametric coordinate ξ
- `η::Float64`: Isoparametric coordinate η

# Returns
- `Vector{Float64}`: Shape function values at the 4 nodes [N₁, N₂, N₃, N₄]

# Shape functions
N₁ = 1/4 * (1 - ξ) * (1 - η)
N₂ = 1/4 * (1 + ξ) * (1 - η)
N₃ = 1/4 * (1 + ξ) * (1 + η)
N₄ = 1/4 * (1 - ξ) * (1 + η)
"""
function compute_shape_functions(ξ::Float64, η::Float64)
    N = zeros(4)
    N[4] = 0.25 * (1.0 - ξ) * (1.0 - η)
    N[1] = 0.25 * (1.0 + ξ) * (1.0 - η)
    N[2] = 0.25 * (1.0 + ξ) * (1.0 + η)
    N[3] = 0.25 * (1.0 - ξ) * (1.0 + η)
    return N
end

"""
    compute_shape_function_derivatives(ξ::Float64, η::Float64)

Compute the derivatives of shape functions for a 4-noded quad element 
in isoparametric coordinates (ξ, η).

# Arguments
- `ξ::Float64`: Isoparametric coordinate ξ
- `η::Float64`: Isoparametric coordinate η

# Returns
- `Matrix{Float64}`: Derivatives of shape functions [4 nodes, 2 coords]
  - Column 1: ∂N/∂ξ for each node
  - Column 2: ∂N/∂η for each node

# Derivatives
∂N₁/∂ξ = -1/4 * (1 - η),  ∂N₁/∂η = -1/4 * (1 - ξ)
∂N₂/∂ξ =  1/4 * (1 - η),  ∂N₂/∂η = -1/4 * (1 + ξ)
∂N₃/∂ξ =  1/4 * (1 + η),  ∂N₃/∂η =  1/4 * (1 + ξ)
∂N₄/∂ξ = -1/4 * (1 + η),  ∂N₄/∂η =  1/4 * (1 - ξ)
"""
function compute_shape_function_derivatives(ξ::Float64, η::Float64)
    B = zeros(4, 2)
    
    # ∂N/∂ξ (column 1)
    B[4, 1] = -0.25 * (1.0 - η)
    B[1, 1] =  0.25 * (1.0 - η)
    B[2, 1] =  0.25 * (1.0 + η)
    B[3, 1] = -0.25 * (1.0 + η)
    
    # ∂N/∂η (column 2)
    B[4, 2] = -0.25 * (1.0 - ξ)
    B[1, 2] = -0.25 * (1.0 + ξ)
    B[2, 2] =  0.25 * (1.0 + ξ)
    B[3, 2] =  0.25 * (1.0 - ξ)
    
    return B
end

"""
    compute_jacobian(B::Matrix{Float64}, X_nodes::Matrix{Float64})

Compute the Jacobian matrix for the transformation from isoparametric to physical coordinates.

# Arguments
- `B::Matrix{Float64}`: Derivatives of shape functions [4 nodes, 2 coords]
- `X_nodes::Matrix{Float64}`: Physical coordinates of element nodes [4 nodes, 2 coords (x, y)]

# Returns
- `Matrix{Float64}`: Jacobian matrix [2×2]

# Formula
J = Bᵀ · X_nodes
"""
function compute_jacobian(B::Matrix{Float64}, X_nodes::Matrix{Float64})
    # J = B^T * X_nodes
    # B is [4 nodes, 2 coords], so B^T is [2, 4]
    # X_nodes is [4 nodes, 2 coords]
    # Result is [2, 2]
    return B' * X_nodes
end

"""
    inverse_and_determinant(J::Matrix{Float64})

Compute the inverse and determinant of a 2×2 Jacobian matrix.

# Arguments
- `J::Matrix{Float64}`: Jacobian matrix [2×2]

# Returns
- `Tuple{Matrix{Float64}, Float64}`: (Inverse Jacobian [2×2], Determinant)
"""
function inverse_and_determinant(J::Matrix{Float64})
    # For a 2x2 matrix: det(J) = J[1,1]*J[2,2] - J[1,2]*J[2,1]
    detJ = J[1, 1] * J[2, 2] - J[1, 2] * J[2, 1]
    
    if abs(detJ) < 1e-12
        error("Singular Jacobian detected: det(J) = $detJ")
    end
    
    # Inverse of 2x2 matrix
    invJ = zeros(2, 2)
    invJ[1, 1] =  J[2, 2] / detJ
    invJ[1, 2] = -J[1, 2] / detJ
    invJ[2, 1] = -J[2, 1] / detJ
    invJ[2, 2] =  J[1, 1] / detJ
    
    return invJ, detJ
end

"""
    initialize_shape_functions!(mesh; axisymmetric = false)

Initialize shape function data structure with precomputed values at Gauss points
and element-specific Jacobian data. Stores data globally within the module.

# Arguments
- `mesh`: Mesh data structure containing:
  - `num_elements`: Number of elements
  - `elements`: Element connectivity matrix [num_elements, 4 nodes]
  - `coordinates`: Node coordinates [num_nodes, 2 coords]
- `axisymmetric::Bool`: When true, fold the Gauss point radius into the cached quadrature
  measure per Appendix `app:axisymmetric`. The first coordinate is the radius `r` and the
  second is the axial coordinate `z`; this pairing is fixed by Eq. `ax_gravity`, which puts
  gravity along the axis of revolution, and ADSIM's gravity vector already runs along `y`.

# Notes
Uses 2×2 Gaussian quadrature (4 Gauss points) with coordinates:
- ξ₁, η₁ = -0.577350269, -0.577350269
- ξ₂, η₂ = -0.577350269,  0.577350269
- ξ₃, η₃ =  0.577350269,  0.577350269
- ξ₄, η₄ =  0.577350269, -0.577350269

All weights equal 1.0 for 2×2 quadrature.
Data is stored in global variable `shape_funcs` accessible within the module.
"""
function initialize_shape_functions!(mesh; axisymmetric::Bool = false)
    # The radius is the first coordinate. A negative one is meaningless on a meridian
    # section and would flip the sign of every quadrature weight in its element, so it is
    # rejected here rather than left to produce a plausible-looking wrong answer.
    if axisymmetric
        r_min = minimum(@view mesh.coordinates[:, 1])
        if r_min < 0.0
            error("Axisymmetric analysis requires every node at radius r = x >= 0, but the " *
                  "mesh has a node at x = $r_min. The meridian section must lie on one side " *
                  "of the axis of revolution (Appendix: Axisymmetric formulation).")
        end
    end

    # Define Gauss points for 2x2 quadrature
    # ξ = ±1/√3 ≈ ±0.577350269
    gp = 1.0 / sqrt(3.0)
    gauss_points = [
        -gp  -gp;   # Point 1
        -gp   gp;   # Point 2
         gp   gp;   # Point 3
         gp  -gp    # Point 4
    ]
    
    # Weights for 2x2 Gaussian quadrature (all equal to 1.0)
    gauss_weights = ones(4)
    
    # Precompute shape functions and derivatives at each Gauss point
    N_gauss = Vector{Vector{Float64}}(undef, 4)
    B_gauss = Vector{Matrix{Float64}}(undef, 4)
    
    for p in 1:4
        ξ = gauss_points[p, 1]
        η = gauss_points[p, 2]
        N_gauss[p] = compute_shape_functions(ξ, η)
        B_gauss[p] = compute_shape_function_derivatives(ξ, η)
    end
    
    # Initialize arrays for element-specific Jacobian data
    Nelements = mesh.num_elements
    invJ_elements = zeros(Nelements, 4, 2, 2)  # [element, gauss point, 2x2 matrix]
    detJ_elements = zeros(Nelements, 4)         # [element, gauss point]
    weight_elements = zeros(Nelements, 4)       # [element, gauss point]
    dN_dx_elements = zeros(Nelements, 4, 4, 2)  # [element, gauss point, 4 nodes, 2 coords]

    # Compute Jacobian data for each element at each Gauss point (parallelized)
    @threads for e in 1:Nelements
        # Get node indices for this element
        node_ids = mesh.elements[e, :]
        
        # Get physical coordinates of element nodes [4 nodes, 2 coords]
        X_nodes = mesh.coordinates[node_ids, :]
        
        # Compute Jacobian at each Gauss point
        for p in 1:4
            J = compute_jacobian(B_gauss[p], X_nodes)
            invJ, detJ = inverse_and_determinant(J)
            
            # Check for degenerate elements
            if abs(detJ) < 1e-12
                error("Element $e, Gauss point $p: Jacobian determinant = $detJ < 1e-12. Element is degenerate or has bad geometry.")
            end
            
            invJ_elements[e, p, :, :] = invJ
            detJ_elements[e, p] = detJ

            # Quadrature measure. Plane strain integrates over the area element; the
            # axisymmetric formulation integrates over the solid of revolution, whose
            # measure carries the radius (Eq. ax_measures, with the common 2*pi cancelled).
            # r_gp is interpolated by the same shape functions that map the element,
            # because the radius IS one of the two physical coordinates (Eq. ax_radius).
            if axisymmetric
                r_gp = 0.0
                for i in 1:4
                    r_gp += N_gauss[p][i] * X_nodes[i, 1]
                end
                weight_elements[e, p] = r_gp * abs(detJ)
            else
                weight_elements[e, p] = abs(detJ)
            end

            # Physical-coordinate derivatives dN/dx. The chain rule gives
            # dN_i/dx_j = sum_k B[i,k] * dxi_k/dx_j, and compute_jacobian returns
            # J[k,j] = dx_j/dxi_k, so the matrix of dxi_k/dx_j is the TRANSPOSE of
            # inv(J), not inv(J) itself. The two coincide whenever J is diagonal -
            # axis-aligned rectangular elements - which is why every 1D verification
            # case passed while this was missing; on a skewed quad the element failed
            # the gradient patch test, keeping the right principal conductivities but
            # rotating their directions. Cached so the time loop never re-forms it.
            dN_dx_elements[e, p, :, :] = B_gauss[p] * transpose(invJ)
        end
    end

    # An element with every node on the axis has zero measure everywhere and would
    # contribute nothing while still owning nodes, leaving them with no equation. The
    # appendix notes that r_gp > 0 for any element with at least one node off the axis,
    # so a vanishing weight means a degenerate mesh rather than a legitimate limit.
    if axisymmetric
        for e in 1:Nelements
            if maximum(@view weight_elements[e, :]) <= 0.0
                error("Element $e has zero measure in the axisymmetric formulation: all four " *
                      "of its nodes lie on the axis r = 0. Remove it from the meridian section.")
            end
        end
    end

    # Store globally within module
    global shape_funcs = ShapeFunctionData(N_gauss, B_gauss, invJ_elements, detJ_elements,
                                            weight_elements, dN_dx_elements, axisymmetric,
                                            gauss_points, gauss_weights)

    return shape_funcs
end

"""
    get_N(p::Int)

Get precomputed shape functions at Gauss point p.

# Arguments
- `p::Int`: Gauss point index (1-4)

# Returns
- `Vector{Float64}`: Shape function values at the 4 nodes
"""
function get_N(p::Int)
    if shape_funcs === nothing
        error("Shape functions not initialized. Call initialize_shape_functions! first.")
    end
    return shape_funcs.N[p]
end

"""
    get_B(p::Int)

Get precomputed shape function derivatives at Gauss point p.

# Arguments
- `p::Int`: Gauss point index (1-4)

# Returns
- `Matrix{Float64}`: Shape function derivatives [4 nodes, 2 coords]
"""
function get_B(p::Int)
    if shape_funcs === nothing
        error("Shape functions not initialized. Call initialize_shape_functions! first.")
    end
    return shape_funcs.B[p]
end

"""
    get_invJ(e::Int, p::Int)

Get precomputed inverse Jacobian for element e at Gauss point p.

# Arguments
- `e::Int`: Element index
- `p::Int`: Gauss point index (1-4)

# Returns
- A 2×2 view of the inverse Jacobian matrix

Returns a view rather than a copy: this is called once per element per Gauss point per gas per
time step, so slicing here would allocate in the innermost loop.
"""
function get_invJ(e::Int, p::Int)
    if shape_funcs === nothing
        error("Shape functions not initialized. Call initialize_shape_functions! first.")
    end
    return @view shape_funcs.invJ[e, p, :, :]
end

"""
    get_dN_dx(e::Int, p::Int)

Get the precomputed physical-coordinate shape function derivatives dN/dx = B · J^-1 for
element e at Gauss point p.

# Arguments
- `e::Int`: Element index
- `p::Int`: Gauss point index (1-4)

# Returns
- A [4 nodes, 2 coords] view of the derivative matrix

This replaces the `B * invJ` product previously formed at each use site. Returns a view, so it
allocates nothing.
"""
function get_dN_dx(e::Int, p::Int)
    if shape_funcs === nothing
        error("Shape functions not initialized. Call initialize_shape_functions! first.")
    end
    return @view shape_funcs.dN_dx[e, p, :, :]
end

"""
    get_detJ(e::Int, p::Int)

Get precomputed Jacobian determinant for element e at Gauss point p.

# Arguments
- `e::Int`: Element index
- `p::Int`: Gauss point index (1-4)

# Returns
- `Float64`: Determinant of Jacobian
"""
function get_detJ(e::Int, p::Int)
    if shape_funcs === nothing
        error("Shape functions not initialized. Call initialize_shape_functions! first.")
    end
    return shape_funcs.detJ[e, p]
end

"""
    get_weight(e::Int, p::Int)

Get the precomputed quadrature measure for element e at Gauss point p: `|det J|` under plane
strain, `r_gp |det J|` under the axisymmetric formulation.

# Arguments
- `e::Int`: Element index
- `p::Int`: Gauss point index (1-4)

# Returns
- `Float64`: Measure to multiply the Gauss weight by

This is what every quadrature sum in the solver must use — lumped masses, all nodal fluxes, and
both sums of the point-to-grid velocity transfer. `get_detJ` remains the Jacobian determinant and
must not be substituted for it: on an axisymmetric run the two differ by the radius factor, and
using the wrong one produces a plane-strain answer on a mesh the user declared axisymmetric.
"""
function get_weight(e::Int, p::Int)
    if shape_funcs === nothing
        error("Shape functions not initialized. Call initialize_shape_functions! first.")
    end
    return shape_funcs.weight[e, p]
end

"""
    is_axisymmetric()

True when the cached quadrature measure carries the radius factor.
"""
function is_axisymmetric()
    if shape_funcs === nothing
        error("Shape functions not initialized. Call initialize_shape_functions! first.")
    end
    return shape_funcs.axisymmetric
end

end # module ShapeFunctions
