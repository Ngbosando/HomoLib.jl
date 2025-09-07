# =============================================
# Triangular Elements 
# =============================================

# Helper for barycentric coordinates
ζ(ξ, η) = 1 - ξ - η

# Tri3 & Tri6 

function shape_functions(::Tri3, ξ::T, η::T) where {T<:Real}
    N = SVector{3, T}(ζ(ξ, η), ξ, η)
    dN_dξ = SVector{3, T}(-1, 1, 0)
    dN_dη = SVector{3, T}(-1, 0, 1)
    dN_dζ = SVector{3, T}( zero(SVector{3, T}))
    return (N, dN_dξ, dN_dη, dN_dζ)
end

function shape_functions(::Tri6, ξ::T, η::T) where {T<:Real}
    z = ζ(ξ, η)
    N = SVector{6, T}(
        z * (2z - 1), ξ * (2ξ - 1), η * (2η - 1),
        4z * ξ, 4ξ * η, 4η * z
    )
    dN_dξ = SVector{6, T}(-3 + 4ξ + 4η, 4ξ - 1, 0, 4z - 4ξ, 4η, -4η)
    dN_dη = SVector{6, T}(-3 + 4ξ + 4η, 0, 4η - 1, -4ξ, 4ξ, 4z - 4η)
    dN_dζ = SVector{6, T}( zero(SVector{6, T}))
    return (N, dN_dξ, dN_dη, dN_dζ)
end


"""
    generate_gmsh_tri_map(p::Int)

Generates the permutation vector to reorder nodes from a simple lexicographic
order to the GMSH standard order for a triangle of polynomial order `p`.
"""
function generate_gmsh_tri_map(p::Int)
    n_nodes = (p + 1) * (p + 2) ÷ 2
    # Map from barycentric powers (i, j, k) to lexicographic index
    lex_map = Dict{Tuple{Int, Int, Int}, Int}()
    idx = 1
    for k in 0:p, j in 0:(p-k)
        i = p - j - k
        lex_map[(i, j, k)] = idx
        idx += 1
    end

    # Build the GMSH ordering by looking up indices in the lex_map
    gmsh_order = Int[]
    # Vertices
    push!(gmsh_order, lex_map[(p, 0, 0)])
    push!(gmsh_order, lex_map[(0, p, 0)])
    push!(gmsh_order, lex_map[(0, 0, p)])
    # Edges
    for i in 1:p-1; push!(gmsh_order, lex_map[(p-i, i, 0)]); end # Edge 1-2
    for i in 1:p-1; push!(gmsh_order, lex_map[(0, p-i, i)]); end # Edge 2-3
    for i in 1:p-1; push!(gmsh_order, lex_map[(i, 0, p-i)]); end # Edge 3-1
    # Interior
    if p > 2
        for k in 1:p-2, j in 1:(p-k-1)
            i = p - j - k
            push!(gmsh_order, lex_map[(i, j, k)])
        end
    end
    return SVector{n_nodes, Int}(gmsh_order)
end

const TRI_GMSH_REORDER_MAP = Dict(
    1 => generate_gmsh_tri_map(1),
    2 => generate_gmsh_tri_map(2),
    3 => generate_gmsh_tri_map(3), # Tri10
    4 => generate_gmsh_tri_map(4), # Tri15
    5 => generate_gmsh_tri_map(5)  # Tri21
)

"""
    tri_poly_shape_functions(ξ, η, ::Val{P})

Computes shape functions for a triangle of order P in a systematic
(lexicographic) order using ForwardDiff for exact derivatives.
"""
function tri_poly_shape_functions(ξ::T, η::T, ::Val{P}) where {T, P}
    n_nodes = (P + 1) * (P + 2) ÷ 2
    
    # Use ForwardDiff to get values and derivatives simultaneously
    ξ_dual = ForwardDiff.Dual(ξ, (one(T), zero(T)))
    η_dual = ForwardDiff.Dual(η, (zero(T), one(T)))
    L1, L2, L3 = 1 - ξ_dual - η_dual, ξ_dual, η_dual

    # Compute values in lexicographic order
    N_lex_duals = MVector{n_nodes, typeof(ξ_dual)}(undef)
    idx = 1
    for k in 0:P, j in 0:(P-k)
        i = P - j - k
        # Using precomputed multinomial coefficients is faster than calling factorial
        # For simplicity here, we compute it on the fly.
        coeff = factorial(P) ÷ (factorial(i) * factorial(j) * factorial(k))
        N_lex_duals[idx] = coeff * L1^i * L2^j * L3^k
        idx += 1
    end
    
    # Extract values and partials
    N_lex = SVector{n_nodes, T}(ForwardDiff.value.(N_lex_duals))
    dN_dξ_lex = SVector{n_nodes, T}(ForwardDiff.partials.(N_lex_duals, 1))
    dN_dη_lex = SVector{n_nodes, T}(ForwardDiff.partials.(N_lex_duals, 2))
    dN_dζ_lex = SVector{n_nodes, T}( zero(SVector{n_nodes, T}))
    return N_lex, dN_dξ_lex, dN_dη_lex, dN_dζ_lex
end

function shape_functions(::Tri10, ξ::T, η::T) where {T<:Real}
    P = 3 # Cubic element
    reorder_map = TRI_GMSH_REORDER_MAP[P]
    N_lex, dN_dξ_lex, dN_dη_lex, dN_dζ_lex = tri_poly_shape_functions(ξ, η, Val(P))
     
    # Reorder the results to match GMSH numbering
    return (N_lex[reorder_map], dN_dξ_lex[reorder_map], dN_dη_lex[reorder_map], dN_dζ_lex)
end

function shape_functions(::Tri15, ξ::T, η::T) where {T<:Real}
    P = 4 # Quartic element
    reorder_map = TRI_GMSH_REORDER_MAP[P]
    N_lex, dN_dξ_lex, dN_dη_lex, dN_dζ_lex = tri_poly_shape_functions(ξ, η, Val(P))
  
    # Reorder the results to match GMSH numbering
    return (N_lex[reorder_map], dN_dξ_lex[reorder_map], dN_dη_lex[reorder_map], dN_dζ_lex)
end

function shape_functions(::Tri21, ξ::T, η::T) where {T<:Real}
    P = 5 # Quintic element
    reorder_map = TRI_GMSH_REORDER_MAP[P]
    N_lex, dN_dξ_lex, dN_dη_lex, dN_dζ_lex = tri_poly_shape_functions(ξ, η, Val(P))
    
    # Reorder the results to match GMSH numbering
    return (N_lex[reorder_map], dN_dξ_lex[reorder_map], dN_dη_lex[reorder_map], dN_dζ_lex)
end