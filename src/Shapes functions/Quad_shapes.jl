# =============================================
# Quadrilateral Elements (Your Logic, Optimized)
# =============================================

function gmsh_ordering(p::Int)
    n = p + 1 # Convert polynomial order to nodes per edge
    idx = Vector{Int}(undef, n^2)
    stack = [(0, 0, n)]
    pos = 1
    while !isempty(stack)
        (i0, j0, nlayer) = popfirst!(stack)
        if nlayer == 1
            idx[pos] = j0 * n + i0 + 1
            pos += 1
        else
            # Corners (counter-clockwise)
            idx[pos] = j0 * n + i0 + 1;                   pos += 1  # Bottom-left
            idx[pos] = j0 * n + (i0+nlayer-1) + 1;         pos += 1  # Bottom-right
            idx[pos] = (j0+nlayer-1)*n + (i0+nlayer-1) + 1; pos += 1  # Top-right
            idx[pos] = (j0+nlayer-1)*n + i0 + 1;           pos += 1  # Top-left

            if nlayer > 2
                # Edges
                for i in 1:(nlayer-2); idx[pos] = j0*n + (i0+i) + 1; pos += 1; end # Bottom
                for j in 1:(nlayer-2); idx[pos] = (j0+j)*n + (i0+nlayer-1) + 1; pos += 1; end # Right
                for i in (nlayer-2):-1:1; idx[pos] = (j0+nlayer-1)*n + (i0+i) + 1; pos += 1; end # Top
                for j in (nlayer-2):-1:1; idx[pos] = (j0+j)*n + i0 + 1; pos += 1; end # Left
                # Recursively process interior
                pushfirst!(stack, (i0+1, j0+1, nlayer-2))
            end
        end
    end
    return idx
end

const QUAD_ORDERS = 2:6
const QUAD_NODES = Dict(p => SVector{p, Float64}(range(-1.0, 1.0, length=p)) for p in QUAD_ORDERS)
const QUAD_GMSH_REORDER_MAP = Dict(p => SVector{(p)^2, Int}(gmsh_ordering(p-1)) for p in QUAD_ORDERS)

function shape_functions_QuadN(ξ::T, η::T, ::Val{P}) where {T<:Real, P}
    nodes_1d = QUAD_NODES[P]
    reorder_map = QUAD_GMSH_REORDER_MAP[P]

    # Calls the refactored, non-allocating version of your Lagrange function
    Lξ, dLξ = lagrange_shape_functions(nodes_1d, ξ)
    Lη, dLη = lagrange_shape_functions(nodes_1d, η)
    
    n_nodes = P * P
    
    # Use MMatrix to build values without heap allocations
    N_lex_mat = MMatrix{P, P, T}(undef)
    dN_dξ_lex_mat = MMatrix{P, P, T}(undef)
    dN_dη_lex_mat = MMatrix{P, P, T}(undef)

    # Your original loop structure for filling
    for j in 1:P
        for i in 1:P
            lex_idx = (j-1)*P + i
            N_lex_mat[lex_idx]     = Lξ[i] * Lη[j]
            dN_dξ_lex_mat[lex_idx] = dLξ[i] * Lη[j]
            dN_dη_lex_mat[lex_idx] = Lξ[i] * dLη[j]
        end
    end

    N_lex     = SVector{n_nodes, T}(N_lex_mat)
    dN_dξ_lex = SVector{n_nodes, T}(dN_dξ_lex_mat)
    dN_dη_lex = SVector{n_nodes, T}(dN_dη_lex_mat)
    dN_dζ_lex = SVector{n_nodes, T}( zero(SVector{n_nodes, T}))
    
    # Reorder using the map generated from YOUR gmsh_ordering function
    return (N_lex[reorder_map], dN_dξ_lex[reorder_map], dN_dη_lex[reorder_map], dN_dζ_lex)
end

shape_functions(::Quad4, ξ::T, η::T) where T = shape_functions_QuadN(ξ, η, Val(2))
shape_functions(::Quad9, ξ::T, η::T) where T = shape_functions_QuadN(ξ, η, Val(3))
shape_functions(::Quad16, ξ::T, η::T) where T = shape_functions_QuadN(ξ, η, Val(4))
shape_functions(::Quad25, ξ::T, η::T) where T = shape_functions_QuadN(ξ, η, Val(5))
shape_functions(::Quad36, ξ::T, η::T) where T = shape_functions_QuadN(ξ, η, Val(6))