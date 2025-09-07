
# =============================================
# Tetrahedral 
# =============================================

function _tet_nodes(::Val{P})::SVector{((P+1)*(P+2)*(P+3)÷6), NTuple{4,Int}} where {P}
    N = (P + 1) * (P + 2) * (P + 3) ÷ 6
    nodes = MVector{N, NTuple{4,Int}}(undef)
    idx = 1
    @inbounds for i in 0:P, j in 0:(P - i), k in 0:(P - i - j)
        l = P - i - j - k
        nodes[idx] = (l, k, j, i)   # (a,b,c,d) for (λ1,λ2,λ3,λ4)
        idx += 1
    end
    return SVector{N, NTuple{4,Int}}(nodes)
end

function product_term(a::Int, p::Int, λ::T) where {T<:Real}
    a == 0 && return one(T)
    pT = T(p)
    res = one(T)
    @inbounds @simd for k in 0:(a-1)
        res *= (pT*λ - T(k)) / T(a - k)
    end
    return res
end

function d_product_term(a::Int, p::Int, λ::T) where {T<:Real}
    a == 0 && return zero(T)
    pT = T(p)
    s  = zero(T)
    @inbounds for j in 0:(a - 1)
        term = one(T)
        # no `continue` in a @simd loop
        @inbounds @simd for k in 0:(a - 1)
            num   = pT*λ - T(k)
            denom = T(a - k)
            # skip k==j by multiplying by 1 (branch allowed; `continue` is not)
            factor = (k == j) ? one(T) : (num / denom)
            term *= factor
        end
        s += pT * term / T(a - j)
    end
    return s
end

function shape_functions_TetN(ξ_in::T, η_in::T, ζ_in::T, ::Val{P}) where {T<:Real, P}
    Nn = (P + 1) * (P + 2) * (P + 3) ÷ 6
    nodes = _tet_nodes(Val(P))                     # SVector{Nn, NTuple{4,Int}}

    # ensure arithmetic stays in T
    oneT = one(T)
    ξ, η, ζ = ξ_in, η_in, ζ_in
    λ1 = oneT - ξ - η - ζ
    λ2 = ξ
    λ3 = η
    λ4 = ζ

    N      = MVector{Nn, T}(undef)
    dN_dξ  = MVector{Nn, T}(undef)
    dN_dη  = MVector{Nn, T}(undef)
    dN_dζ  = MVector{Nn, T}(undef)

    @inbounds for (idx, (a,b,c,d)) in pairs(nodes)
        # values
        P1 = product_term(a, P, λ1)
        P2 = product_term(b, P, λ2)
        P3 = product_term(c, P, λ3)
        P4 = product_term(d, P, λ4)
        val = P1 * P2 * P3 * P4
        N[idx] = val

        # derivatives (dλ1/dξ = dλ1/dη = dλ1/dζ = -1; others per coordinate)
        dP1 = d_product_term(a, P, λ1)
        dP2 = d_product_term(b, P, λ2)
        dP3 = d_product_term(c, P, λ3)
        dP4 = d_product_term(d, P, λ4)

        # ∂/∂ξ: λ1'=-1, λ2'=+1, λ3'=0, λ4'=0
        dN_dξ[idx] = (-dP1)*P2*P3*P4 + P1*( dP2)*P3*P4

        # ∂/∂η: λ1'=-1, λ2'=0, λ3'=+1, λ4'=0
        dN_dη[idx] = (-dP1)*P2*P3*P4 + P1*( 0   )*P3*P4 + P1*P2*( dP3)*P4

        # ∂/∂ζ: λ1'=-1, λ2'=0, λ3'=0, λ4'=+1
        dN_dζ[idx] = (-dP1)*P2*P3*P4 + P1*( 0   )*P3*P4 + P1*P2*( 0   )*P4 + P1*P2*P3*( dP4)
    end

    return SVector{Nn,T}(N), SVector{Nn,T}(dN_dξ), SVector{Nn,T}(dN_dη), SVector{Nn,T}(dN_dζ)
end

shape_functions(::Tet4,  ξ::T, η::T, ζ::T) where {T<:Real} = shape_functions_TetN(ξ, η, ζ, Val(1))
shape_functions(::Tet10, ξ::T, η::T, ζ::T) where {T<:Real} = shape_functions_TetN(ξ, η, ζ, Val(2))
shape_functions(::Tet20, ξ::T, η::T, ζ::T) where {T<:Real} = shape_functions_TetN(ξ, η, ζ, Val(3))
shape_functions(::Tet35, ξ::T, η::T, ζ::T) where {T<:Real} = shape_functions_TetN(ξ, η, ζ, Val(4))
shape_functions(::Tet56, ξ::T, η::T, ζ::T) where {T<:Real} = shape_functions_TetN(ξ, η, ζ, Val(5))
