
# =============================================
# Lin2 Implementation 
# =============================================


function shape_functions(::Lin1, ξ::T) where T<:Real
    return [
        0.5 * (1 - ξ),
        0.5 * (1 + ξ)
    ]
end
# =============================================
# Lin2 Implementation 
# =============================================
# function shape_functions(::Lin2, ξ::T) where T<:Real
#     return [
#         0.5 * ξ^2 - 0.5 * ξ,
#         1.0 - ξ^2,
#         0.5 * ξ^2 + 0.5 * ξ
#     ]
# end

# # =============================================
# # Lin3 Implementation 
# # =============================================

# function shape_functions(::Lin3, ξ::T) where T<:Real
#     return [
#         -0.5625 * ξ^3 + 0.5625 * ξ^2 + 0.0625 * ξ - 0.0625,
#         0.5625 * ξ^3 + 0.5625 * ξ^2 - 0.0625 * ξ - 0.0625,
#         1.6875 * ξ^3 - 0.5625 * ξ^2 - 1.6875 * ξ + 0.5625,
#         -1.6875 * ξ^3 - 0.5625 * ξ^2 + 1.6875 * ξ + 0.5625
#     ]
# end

# # =============================================
# # Lin4 Implementation 
# # =============================================

# function shape_functions(::Lin4, ξ::T) where T<:Real
#     return [
#         0.666666666666667 * ξ^4 - 0.666666666666667 * ξ^3 - 0.166666666666667 * ξ^2 + 0.166666666666667 * ξ,
#         0.666666666666667 * ξ^4 + 0.666666666666667 * ξ^3 - 0.166666666666667 * ξ^2 - 0.166666666666667 * ξ,
#         -2.66666666666667 * ξ^4 + 1.33333333333333 * ξ^3 + 2.66666666666667 * ξ^2 - 1.33333333333333 * ξ,
#         -2.66666666666667 * ξ^4 - 1.33333333333333 * ξ^3 + 2.66666666666667 * ξ^2 + 1.33333333333333 * ξ,
#         4.0 * ξ^4 - 5.0 * ξ^2 + 1.0
#     ]
# end

# # =============================================
# # Lin5 Implementation 
# # =============================================


# function shape_functions(::Lin5, ξ::T) where T<:Real
#     return [
#         -0.813802083333333 * ξ^5 + 0.813802083333333 * ξ^4 + 0.325520833333333 * ξ^3 - 0.325520833333333 * ξ^2 - 0.01171875 * ξ + 0.01171875,
#         0.813802083333333 * ξ^5 + 0.813802083333333 * ξ^4 - 0.325520833333334 * ξ^3 - 0.325520833333333 * ξ^2 + 0.01171875 * ξ + 0.01171875,
#         4.06901041666667 * ξ^5 - 2.44140625 * ξ^4 - 4.23177083333333 * ξ^3 + 2.5390625 * ξ^2 + 0.162760416666667 * ξ - 0.09765625,
#         -4.06901041666667 * ξ^5 - 2.44140625 * ξ^4 + 4.23177083333333 * ξ^3 + 2.5390625 * ξ^2 - 0.162760416666667 * ξ - 0.09765625,
#         -8.13802083333334 * ξ^5 + 1.62760416666667 * ξ^4 + 11.0677083333333 * ξ^3 - 2.21354166666667 * ξ^2 - 2.9296875 * ξ + 0.5859375,
#         8.13802083333333 * ξ^5 + 1.62760416666667 * ξ^4 - 11.0677083333333 * ξ^3 - 2.21354166666667 * ξ^2 + 2.9296875 * ξ + 0.5859375
#     ]
# end




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