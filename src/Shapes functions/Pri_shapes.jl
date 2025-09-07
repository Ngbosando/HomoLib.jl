# 1D equally-spaced nodes in [-1,1] 
linspace_svector(::Type{T}, ::Val{P}) where {T, P} =
    SVector{P+1,T}(ntuple(j -> T(-1 + 2*(j-1)/P), P+1))

# Triangle Bernstein basis and derivatives 
function _triangle_bernstein(::Val{P}, ξ::T, η::T) where {P, T<:Real}
    L1 = one(T) - ξ - η
    L2 = ξ
    L3 = η
    Ntri = (P+1)*(P+2) ÷ 2
    N     = MVector{Ntri,T}(undef)
    dN_ξ  = MVector{Ntri,T}(undef)
    dN_η  = MVector{Ntri,T}(undef)

    idx = 1
    @inbounds for i in 0:P
        for j in 0:(P - i)
            k = P - i - j
            # multinomial coefficient
            C = T(factorial(P) ÷ (factorial(i)*factorial(j)*factorial(k)))
            L1k = k == 0 ? one(T) : L1^k
            L2i = i == 0 ? one(T) : L2^i
            L3j = j == 0 ? one(T) : L3^j
            Nij = C * L1k * L2i * L3j

            # d/dξ: dL1/dξ = -1, dL2/dξ = 1, dL3/dξ = 0
            dξ = C * ( (k==0 ? zero(T) : (-k)*L1^(k-1)) * L2i * L3j +
                       (i==0 ? zero(T) :  i *L1k * L2^(i-1) * L3j) )

            # d/dη: dL1/dη = -1, dL2/dη = 0, dL3/dη = 1
            dη = C * ( (k==0 ? zero(T) : (-k)*L1^(k-1)) * L2i * L3j +
                       (j==0 ? zero(T) :  j *L1k * L2i * L3^(j-1)) )

            N[idx]    = Nij
            dN_ξ[idx] = dξ
            dN_η[idx] = dη
            idx += 1
        end
    end
    return SVector{Ntri,T}(N), SVector{Ntri,T}(dN_ξ), SVector{Ntri,T}(dN_η)
end

# Identity reordering placeholder 
function gmsh_ordering_prism(::Val{P}) where {P}
    Ntri = (P+1)*(P+2) ÷ 2
    NN   = (P+1) * Ntri
    return SVector{NN,Int}(ntuple(i->i, NN))
end

# ===============================
# Generic Prism of order P
# ===============================
function shapeFunctions_PrismN(ξ::T, η::T, ζ::T, ::Val{P}) where {P, T<:Real}
# 1D Lagrange on [-1,1] with P+1 nodes (static)
nodes_ζ = linspace_svector(T, Val(P))
Lζ, dLζ = lagrange_shape_functions(nodes_ζ, ζ)        # each SVector{P+1,T}

# Triangle (order P) static basis
Ntri, dNtri_ξ, dNtri_η = _triangle_bernstein(Val(P), ξ, η)  # each SVector{Ntri,T}

Ntri_len = length(Ntri)
NN = (P+1) * Ntri_len

N     = MVector{NN,T}(undef)
dN_ξ  = MVector{NN,T}(undef)
dN_η  = MVector{NN,T}(undef)
dN_ζ  = MVector{NN,T}(undef)

idx = 1
@inbounds for k in 1:(P+1)
    lk  = Lζ[k]
    dlk = dLζ[k]
    for i in 1:Ntri_len
        nval = Ntri[i]
        N[idx]    = nval       * lk
        dN_ξ[idx] = dNtri_ξ[i] * lk
        dN_η[idx] = dNtri_η[i] * lk
        dN_ζ[idx] = nval       * dlk
        idx += 1
    end
end

# Apply Gmsh ordering (identity by default; replace if you have the perm)
perm = gmsh_ordering_prism(Val(P))
return SVector{NN,T}(N)[perm],
        SVector{NN,T}(dN_ξ)[perm],
        SVector{NN,T}(dN_η)[perm],
        SVector{NN,T}(dN_ζ)[perm]
end

# ===============================
# Pri6 
# ===============================
function shape_functions(::Pri6, ξ::T, η::T, ζ::T) where {T<:Real}
# classic linear prism shape functions (6 nodes)
half = T(0.5)


N = MVector{6,T}(undef)
N[1] = (1 - ξ - η)*(1 - ζ)*half
N[2] = ξ*(1 - ζ)*half
N[3] = η*(1 - ζ)*half
N[4] = (1 - ξ - η)*(1 + ζ)*half
N[5] = ξ*(1 + ζ)*half
N[6] = η*(1 + ζ)*half

dξ = MVector{6,T}(undef)
dη = MVector{6,T}(undef)
dζ = MVector{6,T}(undef)

dξ[1] = -(1 - ζ)*half; dη[1] = -(1 - ζ)*half; dζ[1] = -(1 - ξ - η)*half
dξ[2] =  (1 - ζ)*half; dη[2] =  zero(T);        dζ[2] = -ξ*half
dξ[3] =  zero(T);        dη[3] =  (1 - ζ)*half; dζ[3] = -η*half
dξ[4] = -(1 + ζ)*half; dη[4] = -(1 + ζ)*half; dζ[4] =  (1 - ξ - η)*half
dξ[5] =  (1 + ζ)*half; dη[5] =  zero(T);        dζ[5] =  ξ*half
dξ[6] =  zero(T);        dη[6] =  (1 + ζ)*half; dζ[6] =  η*half

return SVector{6,T}(N), SVector{6,T}(dξ), SVector{6,T}(dη), SVector{6,T}(dζ)
end

# ===============================
# Pri18 (quadratic, P=2) 
# ===============================
shape_functions(::Pri18, ξ::T, η::T, ζ::T) where {T<:Real} =
shapeFunctions_PrismN(ξ, η, ζ, Val(2))

# ===============================
# Higher orders 
# ===============================
shape_functions(::Pri36, ξ::T, η::T, ζ::T) where {T<:Real} =
shapeFunctions_PrismN(ξ, η, ζ, Val(3))
shape_functions(::Pri56, ξ::T, η::T, ζ::T) where {T<:Real} =
shapeFunctions_PrismN(ξ, η, ζ, Val(4))
shape_functions(::Pri78, ξ::T, η::T, ζ::T) where {T<:Real} =
shapeFunctions_PrismN(ξ, η, ζ, Val(5))

# ===============================
# reference coords generator 
# ===============================
function reference_coords_prism(::Val{P}) where {P}
    pts = Vector{SVector{3,Float64}}()
    for k in 0:P
        ζ = k / P
        for i in 0:P, j in 0:(P - i)
            ξ = i / P
            η = j / P
            push!(pts, SVector{3,Float64}(ξ, η, ζ))
        end
    end
    return pts
end
