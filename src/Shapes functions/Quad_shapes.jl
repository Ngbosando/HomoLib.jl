# =============================================
# Quadrilateral Element Structs
# =============================================



# =============================================
# Quad4 Implementation (Bilinear)
# =============================================
function quad_nodes_indices_gmsh_order(p)
    indices = []

    # 1. Corner nodes
    push!(indices, (0, 0))     # 1
    push!(indices, (p, 0))     # 2
    push!(indices, (p, p))     # 3
    push!(indices, (0, p))     # 4

    # 2. Edge 1-2 (bottom)
    for i in 1:(p - 1)
        push!(indices, (i, 0))
    end

    # 3. Edge 2-3 (right)
    for j in 1:(p - 1)
        push!(indices, (p, j))
    end

    # 4. Edge 3-4 (top)
    for i in (p - 1):-1:1
        push!(indices, (i, p))
    end

    # 5. Edge 4-1 (left)
    for j in (p - 1):-1:1
        push!(indices, (0, j))
    end

    # 6. Interior
    for j in 1:(p - 1)
        for i in 1:(p - 1)
            push!(indices, (i, j))
        end
    end

    return indices
end


function lagrange_shape_functions(nodes::Vector{T}, ξ::T) where T<:Real
    n = length(nodes)
    L = zeros(T, n)
    dL = zeros(T, n)

  

    for i in 1:n
        Li = one(T)
        dLi = zero(T)
        for j in 1:n
            if j != i
                Li *= (ξ - nodes[j]) / (nodes[i] - nodes[j])
                term = one(T) / (nodes[i] - nodes[j])
                prod = one(T)
                for k in 1:n
                    if k != i && k != j
                        prod *= (ξ - nodes[k]) / (nodes[i] - nodes[k])
                    end
                end
                dLi += term * prod
            end
        end
        L[i] = Li
        dL[i] = dLi
    end
    return L, dL
end

function gmsh_ordering(n::Int)
    idx = Vector{Int}(undef, n^2)
    stack = [(0, 0, n)]
    pos = 1
    while !isempty(stack)
        (i0, j0, nlayer) = popfirst!(stack)
        if nlayer == 1
            # Handle single node: assign once
            idx[pos] = j0 * n + i0 + 1
            pos += 1
        else
            # Corners (counter-clockwise)
            idx[pos] = j0 * n + i0 + 1;          pos += 1  # Bottom-left
            idx[pos] = j0 * n + (i0+nlayer-1) + 1; pos += 1  # Bottom-right
            idx[pos] = (j0+nlayer-1)*n + (i0+nlayer-1) + 1; pos += 1  # Top-right
            idx[pos] = (j0+nlayer-1)*n + i0 + 1;          pos += 1  # Top-left

            if nlayer > 2
                # Edges (excluding corners)
                # Bottom edge (left to right)
                for i in 1:(nlayer-2)
                    idx[pos] = j0 * n + (i0+i) + 1; pos += 1
                end
                # Right edge (bottom to top)
                for j in 1:(nlayer-2)
                    idx[pos] = (j0+j)*n + (i0+nlayer-1) + 1; pos += 1
                end
                # Top edge (right to left)
                for i in (nlayer-2):-1:1
                    idx[pos] = (j0+nlayer-1)*n + (i0+i) + 1; pos += 1
                end
                # Left edge (top to bottom)
                for j in (nlayer-2):-1:1
                    idx[pos] = (j0+j)*n + i0 + 1; pos += 1
                end

                # Recursively process interior
                pushfirst!(stack, (i0+1, j0+1, nlayer-2))
            end
        end
    end
    return idx
end



function shapeFunctions_QuadN(ξ::T, η::T,order::Int) where T<:Real
    # @assert 2 <= order <= 5 "Only supports orders 2 to 5"

    nodes_1d = range(-1.0, 1.0, length=order) |> collect

    Lξ, dLξ = lagrange_shape_functions(nodes_1d, ξ)
    Lη, dLη = lagrange_shape_functions(nodes_1d, η)
   
    n = order^2
    N      = Vector{T}(undef, n)
    dN_dξ  = similar(N)
    dN_dη  = similar(N)

    idx = 1
    @inbounds for j in 1:order
        Lj = Lη[j]; dLj = dLη[j]
        for i in 1:order
           li = Lξ[i]; dli = dLξ[i]
            N[idx]     = li * Lj
            dN_dξ[idx] = dli * Lj
            dN_dη[idx] = li * dLj
            idx += 1
        end
    end
   
    reorder = gmsh_ordering(order)
    return N[reorder], dN_dξ[reorder], dN_dη[reorder]
end




function shape_functions(::Quad4, ξ::T, η::T) where T<:Real
    N, dN_dξ, dN_dη = shapeFunctions_QuadN(ξ, η, 2)
    return N, dN_dξ, dN_dη
end

function shape_functions(::Quad9, ξ::T, η::T) where T<:Real
    N, dN_dξ, dN_dη = shapeFunctions_QuadN(ξ, η, 3)
    return N, dN_dξ, dN_dη
end
function shape_functions(::Quad16, ξ::T, η::T) where T<:Real
    N, dN_dξ, dN_dη = shapeFunctions_QuadN(ξ, η, 4)
    return N, dN_dξ, dN_dη
end
function shape_functions(::Quad25, ξ::T, η::T) where T<:Real
    N, dN_dξ, dN_dη = shapeFunctions_QuadN(ξ, η, 5)
    return N, dN_dξ, dN_dη
end

function shape_functions(::Quad36, ξ::T, η::T) where T<:Real
    N, dN_dξ, dN_dη = shapeFunctions_QuadN(ξ, η, 6)
    return N, dN_dξ, dN_dη
end
