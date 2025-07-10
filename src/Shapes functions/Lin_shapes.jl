
# =============================================
# Lin2 Implementation 
# =============================================


function gmsh_ordering_line(count::Int)::Vector{Int}
    ordering = Int[]
    push!(ordering, 1)                       # start node
    push!(ordering, count)                   # end node
    if count > 2
        for i in 2:count-1
            push!(ordering, i)               # internal edge nodes
        end
    end
    return ordering
end

function shapeFunctions_LineN(ξ::T, count::Int) where T<:Real
    nodes_1d = collect(range(-1.0, 1.0, length=count))
    L, dL = lagrange_shape_functions(nodes_1d, ξ)
    reorder = gmsh_ordering_line(count)
  
    return L[reorder], dL[reorder]
end

# Dispatch functions for specific Line types
function shape_functions(::Lin2, ξ::T) where T<:Real
    return shapeFunctions_LineN(ξ, 2)
end

function shape_functions(::Lin3, ξ::T) where T<:Real
    return shapeFunctions_LineN(ξ, 3)
end

function shape_functions(::Lin4, ξ::T) where T<:Real
    return shapeFunctions_LineN(ξ, 4)
end

function shape_functions(::Lin5, ξ::T) where T<:Real
    return shapeFunctions_LineN(ξ, 5)
end

function shape_functions(::Lin6, ξ::T) where T<:Real
    return shapeFunctions_LineN(ξ, 6)
end