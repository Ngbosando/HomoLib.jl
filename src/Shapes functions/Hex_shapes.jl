


function shapeFunctions_HexN(ξ::T, η::T, ζ::T, count::Int) where T<:Real
    nodes_1d = collect(range(-1.0, 1.0, length=count))
    Lξ, dLξ = lagrange_shape_functions(nodes_1d, ξ)
    Lη, dLη = lagrange_shape_functions(nodes_1d, η)
    Lζ, dLζ = lagrange_shape_functions(nodes_1d, ζ)

    n = count^3
    N      = Vector{T}(undef, n)
    dN_dξ  = similar(N)
    dN_dη  = similar(N)
    dN_dζ  = similar(N)

    idx = 1
    @inbounds for k in 1:count, j in 1:count, i in 1:count
        lk = Lζ[k]; lj = Lη[j]; li = Lξ[i]
        N[idx]     = li * lj * lk
        dN_dξ[idx] = dLξ[i] * lj * lk
        dN_dη[idx] = li * dLη[j] * lk
        dN_dζ[idx] = li * lj * dLζ[k]
        idx += 1
    end
    

    # return N[reorder], dN_dξ[reorder], dN_dη[reorder], dN_dζ[reorder]
    return N, dN_dξ, dN_dη, dN_dζ

end

# Dispatch functions for specific Hex types
function shape_functions(::Hex8, ξ::T, η::T, ζ::T) where T<:Real
    return shapeFunctions_HexN(ξ, η, ζ, 2)
end

function shape_functions(::Hex27, ξ::T, η::T, ζ::T) where T<:Real
    return shapeFunctions_HexN(ξ, η, ζ, 3)
end

function shape_functions(::Hex64, ξ::T, η::T, ζ::T) where T<:Real
    return shapeFunctions_HexN(ξ, η, ζ, 4)
end

function shape_functions(::Hex125, ξ::T, η::T, ζ::T) where T<:Real
    return shapeFunctions_HexN(ξ, η, ζ, 5)
end
function shape_functions(::Hex216, ξ::T, η::T, ζ::T) where T<:Real
    N, dN_dξ, dN_dη, dN_dζ = shapeFunctions_HexN(ξ, η, ζ, 6)
    return N, dN_dξ, dN_dη, dN_dζ
end