# Precomputed 1D Gauss–Legendre constants 

const GAUSS_LEGENDRE_RULES_CONST = (
    (SVector(0.0), SVector(2.0)),
    (SVector(-1/√3, 1/√3), SVector(1.0, 1.0)),
    (SVector(-√(3/5), 0.0, √(3/5)), SVector(5/9, 8/9, 5/9)),
    (SVector(-√(3/7 + 2/7*√(6/5)), -√(3/7 - 2/7*√(6/5)), √(3/7 - 2/7*√(6/5)), √(3/7 + 2/7*√(6/5))),
     SVector((18-√30)/36, (18+√30)/36, (18+√30)/36, (18-√30)/36)),
    (SVector(-1/3*√(5+2*√(10/7)), -1/3*√(5-2*√(10/7)), 0.0, 1/3*√(5-2*√(10/7)), 1/3*√(5+2*√(10/7))),
     SVector((322-13*√70)/900, (322+13*√70)/900, 128/225, (322+13*√70)/900, (322-13*√70)/900)),
    (SVector(-0.932469514203152, -0.661209386466265, -0.238619186083197, 0.238619186083197, 0.661209386466265, 0.932469514203152),
     SVector(0.171324492379170, 0.360761573048139, 0.467913934572691, 0.467913934572691, 0.360761573048139, 0.171324492379170))
)

"""
    get_gauss_legendre_rule_vec(order::Int) -> (Vector{Float64}, Vector{Float64})

Type-stable 1D Gauss–Legendre rule. Always returns Vectors.
"""
function get_gauss_legendre_rule_vec(order::Int)
    if 1 <= order <= length(GAUSS_LEGENDRE_RULES_CONST)
        pS, wS = GAUSS_LEGENDRE_RULES_CONST[order]

        return (Vector{Float64}(pS), Vector{Float64}(wS))
    else
        p, w = gausslegendre(order)           
        P = Vector{Float64}(undef, length(p))
        W = Vector{Float64}(undef, length(w))
        @inbounds @simd for i in eachindex(P)
            P[i] = Float64(p[i])
            W[i] = Float64(w[i])
        end
        return (P, W)
    end
end


function get_gauss_legendre_rule(::Val{O}) where {O}
    if 1 <= O <= length(GAUSS_LEGENDRE_RULES_CONST)
        return GAUSS_LEGENDRE_RULES_CONST[O]    # (SVector{O}, SVector{O})
    else
        p, w = gausslegendre(O)
        return (SVector{O,Float64}(p), SVector{O,Float64}(w))
    end
end

# Tensor-product rule builders 

function create_tensor_product_rules(order::Int, ::Val{2})
    p1d, w1d = get_gauss_legendre_rule_vec(order)  # Vector{Float64}, Vector{Float64}
    n = order * order
    points  = Vector{SVector{2,Float64}}(undef, n)
    weights = Vector{Float64}(undef, n)
    idx = 1
    @inbounds for j in 1:order, i in 1:order
        xi = p1d[i]; yi = p1d[j]
        points[idx]  = SVector{2,Float64}(xi, yi)
        weights[idx] = w1d[i] * w1d[j]
        idx += 1
    end
    return (points, weights) :: Tuple{Vector{SVector{2,Float64}}, Vector{Float64}}
end

function create_tensor_product_rules(order::Int, ::Val{3})
    p1d, w1d = get_gauss_legendre_rule_vec(order)
    n = order * order * order
    points  = Vector{SVector{3,Float64}}(undef, n)
    weights = Vector{Float64}(undef, n)
    idx = 1
    @inbounds for k in 1:order, j in 1:order, i in 1:order
        xi = p1d[i]; yi = p1d[j]; zi = p1d[k]
        points[idx]  = SVector{3,Float64}(xi, yi, zi)
        weights[idx] = w1d[i] * w1d[j] * w1d[k]
        idx += 1
    end
    return (points, weights) :: Tuple{Vector{SVector{3,Float64}}, Vector{Float64}}
end

const MAX_PRECOMPUTED_ORDER = 6

const QUADRATURE_RULES_QUAD_CONST =
    Dict{Int, Tuple{Vector{SVector{2,Float64}}, Vector{Float64}}}(
        (o => create_tensor_product_rules(o, Val(2)) for o in 1:MAX_PRECOMPUTED_ORDER)...
    )

const QUADRATURE_RULES_HEX_CONST =
    Dict{Int, Tuple{Vector{SVector{3,Float64}}, Vector{Float64}}}(
        (o => create_tensor_product_rules(o, Val(3)) for o in 1:MAX_PRECOMPUTED_ORDER)...
    )

const QUADRATURE_RULES_QUAD_CACHE =
    Dict{Int, Tuple{Vector{SVector{2,Float64}}, Vector{Float64}}}()

const QUADRATURE_RULES_HEX_CACHE  =
    Dict{Int, Tuple{Vector{SVector{3,Float64}}, Vector{Float64}}}()

# Final integration rule dispatch 

function int_rule(::LinearElement, order::Int)
    p1d, w1d = get_gauss_legendre_rule_vec(order)     # Vector{Float64}, Vector{Float64}
    pts = Vector{SVector{1,Float64}}(undef, length(p1d))
    @inbounds @simd for i in eachindex(p1d)
        pts[i] = SVector{1,Float64}(p1d[i])
    end
    
    return (pts, w1d) :: Tuple{Vector{SVector{1,Float64}}, Vector{Float64}}
end

function int_rule(::QuadrilateralElement, order::Int)
    if haskey(QUADRATURE_RULES_QUAD_CONST, order)
        return QUADRATURE_RULES_QUAD_CONST[order]
    elseif haskey(QUADRATURE_RULES_QUAD_CACHE, order)
        return QUADRATURE_RULES_QUAD_CACHE[order]
    else
        rule = create_tensor_product_rules(order, Val(2))
        QUADRATURE_RULES_QUAD_CACHE[order] = rule
        return rule
    end
end

function int_rule(::HexahedralElement, order::Int)
    if haskey(QUADRATURE_RULES_HEX_CONST, order)
        return QUADRATURE_RULES_HEX_CONST[order]
    elseif haskey(QUADRATURE_RULES_HEX_CACHE, order)
        return QUADRATURE_RULES_HEX_CACHE[order]
    else
        rule = create_tensor_product_rules(order, Val(3))
        QUADRATURE_RULES_HEX_CACHE[order] = rule
        return rule
    end
end
