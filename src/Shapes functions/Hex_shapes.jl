# ==============================================================================
# Hexahedral Elements 
# ==============================================================================

const HEX_ORDERS = 2:6
const HEX_NODES = Dict(p => SVector{p, Float64}(range(-1.0, 1.0, length=p)) for p in HEX_ORDERS)

function lagrange_shape_functions(nodes::SVector{P, T}, ξ::T) where {P, T<:Real}
    L = MMatrix{P, 1, T}(undef)
    dL = MMatrix{P, 1, T}(undef)

    @inbounds for i in 1:P
        Li = one(T)
        dLi = zero(T)
        for j in 1:P
            if j != i
                Li *= (ξ - nodes[j]) / (nodes[i] - nodes[j])
                term = one(T) / (nodes[i] - nodes[j])
                prod_val = one(T)
                for k in 1:P
                    if k != i && k != j
                        prod_val *= (ξ - nodes[k]) / (nodes[i] - nodes[k])
                    end
                end
                dLi += term * prod_val
            end
        end
        L[i] = Li
        dL[i] = dLi
    end
    return SVector(L), SVector(dL)
end

function shape_functions_HexN(ξ::T, η::T, ζ::T, ::Val{P}) where {T<:Real, P}
    nodes_1d = HEX_NODES[P]
    
    Lξ, dLξ = lagrange_shape_functions(nodes_1d, ξ)
    Lη, dLη = lagrange_shape_functions(nodes_1d, η)
    Lζ, dLζ = lagrange_shape_functions(nodes_1d, ζ)

    n_nodes = P * P * P
    
    # Use a 3D mutable MArray on the stack to build values without allocation.
    # This makes the 3D tensor product structure clear.
    N_mat     = MArray{Tuple{P, P, P}, T}(undef)
    dN_dξ_mat = MArray{Tuple{P, P, P}, T}(undef)
    dN_dη_mat = MArray{Tuple{P, P, P}, T}(undef)
    dN_dζ_mat = MArray{Tuple{P, P, P}, T}(undef)

    # Your original loop structure, filling the 3D MArray
    @inbounds for k in 1:P, j in 1:P, i in 1:P
        lk = Lζ[k]; lj = Lη[j]; li = Lξ[i]
        
        N_mat[i, j, k]     = li * lj * lk
        dN_dξ_mat[i, j, k] = dLξ[i] * lj * lk
        dN_dη_mat[i, j, k] = li * dLη[j] * lk
        dN_dζ_mat[i, j, k] = li * lj * dLζ[k]
    end
    
    # Convert the 3D MArrays to flat SVectors (a free operation).
    # The flattening order matches your original `idx` counter.
    N     = SVector{n_nodes, T}(N_mat)
    dN_dξ = SVector{n_nodes, T}(dN_dξ_mat)
    dN_dη = SVector{n_nodes, T}(dN_dη_mat)
    dN_dζ = SVector{n_nodes, T}(dN_dζ_mat)
    
    # Return the results. NO reordering is applied, as you requested.
    return N, dN_dξ, dN_dη, dN_dζ
end

shape_functions(::Hex8,   ξ::T, η::T, ζ::T) where T = shape_functions_HexN(ξ, η, ζ, Val(2))
shape_functions(::Hex27,  ξ::T, η::T, ζ::T) where T = shape_functions_HexN(ξ, η, ζ, Val(3))
shape_functions(::Hex64,  ξ::T, η::T, ζ::T) where T = shape_functions_HexN(ξ, η, ζ, Val(4))
shape_functions(::Hex125, ξ::T, η::T, ζ::T) where T = shape_functions_HexN(ξ, η, ζ, Val(5))
shape_functions(::Hex216, ξ::T, η::T, ζ::T) where T = shape_functions_HexN(ξ, η, ζ, Val(6))