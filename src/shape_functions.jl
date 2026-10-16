"""
    shape_functions.jl

Isoparametric shape functions and Jacobians for 4-noded quad elements.
2x2 Gaussian quadrature (4 Gauss points).
"""

module ShapeFunctions

using Base.Threads

export ShapeFunctionData, initialize_shape_functions!, get_N, get_B, get_invJ, get_detJ, get_dN_dx,
       get_weight, is_axisymmetric

"""
    ShapeFunctionData

Precomputed shape function values and element Jacobian data.

# Fields
- `N::Vector{Vector{Float64}}`: Shape functions at each Gauss point [4 Gauss points][4 nodes]
- `B::Vector{Matrix{Float64}}`: Derivatives of shape functions at each Gauss point [4 Gauss points][4 nodes, 2 coords]
- `invJ::Array{Float64, 4}`: Inverse Jacobian at each element and Gauss point [Nelements, 4 Gauss points, 2, 2]
- `detJ::Matrix{Float64}`: Determinant of Jacobian at each element and Gauss point [Nelements, 4 Gauss points]
- `weight::Matrix{Float64}`: Quadrature measure at each element and Gauss point [Nelements, 4 Gauss points]
- `dN_dx::Array{Float64, 4}`: Shape function derivatives in physical coordinates, B · J^-1,
  at each element and Gauss point [Nelements, 4 Gauss points, 4 nodes, 2 coords]
- `axisymmetric::Bool`: True when `weight` includes the radius factor of Appendix `app:axisymmetric`
- `gauss_points::Matrix{Float64}`: Gauss point coordinates in isoparametric space [4 points, 2 coords (ξ, η)]
- `gauss_weights::Vector{Float64}`: Gauss point weights [4 points]

`dN_dx` is used by every flux term. Stored instead of recomputed from `B` and `invJ`:
element- and Gauss-point-specific but time-invariant.

`weight` is the ONLY field with the plane/axisymmetric distinction. The axisymmetric
formulation replaces the area measure `|det J| W_gp` by `r_gp |det J| W_gp`, with gradients
and stiffness contractions unchanged. `r_gp` included in the cached measure at
initialization: every quadrature sum valid for both formulations, no branch, no extra
arithmetic per time step.

    weight[e, p] = axisymmetric ? r_gp * |det J| : |det J|

`detJ` stored separately (Jacobian determinant only). `get_weight` sites are quadrature
sums (with radius); `get_detJ` sites are geometry (without radius).
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

Compute the 4-node quad shape functions [N₁, N₂, N₃, N₄] at isoparametric coordinates (ξ, η).

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

Compute the 4-node quad shape function derivatives [4 nodes, 2 coords] at (ξ, η).
Column 1: ∂N/∂ξ, column 2: ∂N/∂η.

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

Compute the 2×2 Jacobian of the isoparametric-to-physical mapping.

# Arguments
- `B`: Shape function derivatives [4 nodes, 2 coords]
- `X_nodes`: Physical node coordinates [4 nodes, 2 coords (x, y)]

# Formula
J = Bᵀ · X_nodes
"""
function compute_jacobian(B::Matrix{Float64}, X_nodes::Matrix{Float64})
    # [2, 4] * [4, 2] -> [2, 2]
    return B' * X_nodes
end

"""
    inverse_and_determinant(J::Matrix{Float64})

Return `(inv(J), det(J))` of a 2×2 Jacobian. Errors on a singular Jacobian.
"""
function inverse_and_determinant(J::Matrix{Float64})
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

Precompute Gauss-point shape functions and element Jacobian data. Mutates the module
global `shape_funcs`.

# Arguments
- `mesh`: Uses `num_elements`, `elements` [num_elements, 4] and `coordinates` [num_nodes, 2]
- `axisymmetric`: Include the Gauss-point radius in the cached quadrature measure
  (Appendix `app:axisymmetric`). First coordinate = radius `r`, second = axial `z`;
  fixed by Eq. `ax_gravity` (gravity along the axis of revolution, i.e. `y`).

# Notes
2×2 Gaussian quadrature (4 Gauss points):
- ξ₁, η₁ = -0.577350269, -0.577350269
- ξ₂, η₂ = -0.577350269,  0.577350269
- ξ₃, η₃ =  0.577350269,  0.577350269
- ξ₄, η₄ =  0.577350269, -0.577350269

All weights equal 1.0.
"""
function initialize_shape_functions!(mesh; axisymmetric::Bool = false)
    # Radius = first coordinate. Negative radius is not valid on a meridian section and
    # flips the sign of the element quadrature weights; rejected here.
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

            # Quadrature measure. Plane strain: area element. Axisymmetric: solid of
            # revolution, measure includes the radius (Eq. ax_measures, 2*pi cancelled).
            # r_gp interpolated with the element shape functions, since the radius is a
            # physical coordinate (Eq. ax_radius).
            if axisymmetric
                r_gp = 0.0
                for i in 1:4
                    r_gp += N_gauss[p][i] * X_nodes[i, 1]
                end
                weight_elements[e, p] = r_gp * abs(detJ)
            else
                weight_elements[e, p] = abs(detJ)
            end

            # Physical-coordinate derivatives dN/dx. Chain rule:
            # dN_i/dx_j = sum_k B[i,k] * dxi_k/dx_j, with J[k,j] = dx_j/dxi_k from
            # compute_jacobian, therefore dxi_k/dx_j = TRANSPOSE of inv(J). Identical for
            # diagonal J (axis-aligned rectangles, all 1D verification cases). Without the
            # transpose, skewed quads fail the gradient patch test (correct principal
            # conductivities, rotated directions). Cached for the time loop.
            dN_dx_elements[e, p, :, :] = B_gauss[p] * transpose(invJ)
        end
    end

    # Element with all nodes on the axis: zero measure, nodes without an equation.
    # r_gp > 0 for any element with at least one node off the axis (appendix), therefore
    # zero weight = degenerate mesh.
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

Return the precomputed nodal shape function values at Gauss point `p` (1-4).
"""
function get_N(p::Int)
    if shape_funcs === nothing
        error("Shape functions not initialized. Call initialize_shape_functions! first.")
    end
    return shape_funcs.N[p]
end

"""
    get_B(p::Int)

Return the precomputed isoparametric derivatives [4 nodes, 2 coords] at Gauss point `p` (1-4).
"""
function get_B(p::Int)
    if shape_funcs === nothing
        error("Shape functions not initialized. Call initialize_shape_functions! first.")
    end
    return shape_funcs.B[p]
end

"""
    get_invJ(e::Int, p::Int)

Return a 2×2 view of the inverse Jacobian of element `e` at Gauss point `p` (1-4).

View instead of copy: called per element, Gauss point, gas and time step; a slice would
allocate in the innermost loop.
"""
function get_invJ(e::Int, p::Int)
    if shape_funcs === nothing
        error("Shape functions not initialized. Call initialize_shape_functions! first.")
    end
    return @view shape_funcs.invJ[e, p, :, :]
end

"""
    get_dN_dx(e::Int, p::Int)

Return a [4 nodes, 2 coords] view of the physical derivatives dN/dx = B · J^-1 of element
`e` at Gauss point `p` (1-4). No allocation.
"""
function get_dN_dx(e::Int, p::Int)
    if shape_funcs === nothing
        error("Shape functions not initialized. Call initialize_shape_functions! first.")
    end
    return @view shape_funcs.dN_dx[e, p, :, :]
end

"""
    get_detJ(e::Int, p::Int)

Return the Jacobian determinant of element `e` at Gauss point `p` (1-4).
"""
function get_detJ(e::Int, p::Int)
    if shape_funcs === nothing
        error("Shape functions not initialized. Call initialize_shape_functions! first.")
    end
    return shape_funcs.detJ[e, p]
end

"""
    get_weight(e::Int, p::Int)

Return the quadrature measure of element `e` at Gauss point `p` (1-4), to be multiplied by
the Gauss weight: `|det J|` for plane strain, `r_gp |det J|` for axisymmetry.

Required by every quadrature sum in the solver: lumped masses, nodal fluxes, both sums of
the point-to-grid velocity transfer. Not interchangeable with `get_detJ`: they differ by the
radius factor under axisymmetry, and `get_detJ` gives a plane-strain result.
"""
function get_weight(e::Int, p::Int)
    if shape_funcs === nothing
        error("Shape functions not initialized. Call initialize_shape_functions! first.")
    end
    return shape_funcs.weight[e, p]
end

"""
    is_axisymmetric()

Return `true` when the cached quadrature measure includes the radius factor.
"""
function is_axisymmetric()
    if shape_funcs === nothing
        error("Shape functions not initialized. Call initialize_shape_functions! first.")
    end
    return shape_funcs.axisymmetric
end

end # module ShapeFunctions
