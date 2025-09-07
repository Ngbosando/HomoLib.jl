# ==============================================================================
# Line Elements 
# ==============================================================================

function gmsh_ordering_line(count::Int)::Vector{Int}
    internal = count > 2 ? collect(2:count-1) : Int[]
    return [1, count, internal...]
end

const LINE_ORDERS = 2:6
const LINE_NODES = Dict(p => SVector{p, Float64}(range(-1.0, 1.0, length=p)) for p in LINE_ORDERS)
const LINE_GMSH_REORDER_MAP = Dict(p => SVector{p, Int}(gmsh_ordering_line(p)) for p in LINE_ORDERS)

function shape_functions_LineN(ξ::T, ::Val{P}) where {T<:Real, P}
    nodes_1d = LINE_NODES[P]
    reorder_map = LINE_GMSH_REORDER_MAP[P]
    
    # Calls the fast, non-allocating helper
    L_lex, dL_lex = lagrange_shape_functions(nodes_1d, ξ)
 
    # Reorder using the map generated from YOUR gmsh_ordering_line function
    N     = L_lex[reorder_map]
    dN_dξ = dL_lex[reorder_map]
    
    # Return the standard 4-tuple for compatibility. For 1D elements,
    # the η and ζ derivatives are nothing.
    return (N, dN_dξ, nothing, nothing)
end

shape_functions(::Lin2, ξ::T) where T = shape_functions_LineN(ξ, Val(2))
shape_functions(::Lin3, ξ::T) where T = shape_functions_LineN(ξ, Val(3))
shape_functions(::Lin4, ξ::T) where T = shape_functions_LineN(ξ, Val(4))
shape_functions(::Lin5, ξ::T) where T = shape_functions_LineN(ξ, Val(5))
shape_functions(::Lin6, ξ::T) where T = shape_functions_LineN(ξ, Val(6))