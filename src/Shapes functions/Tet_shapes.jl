# =============================================
# Tetrahedral Element Structs
# =============================================

# function generate_tet_nodes(p::Int)
#     nodes = NTuple{4,Int}[]
#     for a in 0:p
#         for b in 0:(p - a)
#             for c in 0:(p - a - b)
#                 d = p - a - b - c
#                 push!(nodes, (a, b, c, d))
#             end
#         end
#     end
#     return nodes
# end

# function product_term(exponent, p, λ)
#     exponent == 0 && return one(λ)
#     term = one(λ)
#     for k in 0:exponent-1
#         term *= (p*λ - k) / (exponent - k)
#     end
#     return term
# end
# function shapeFunctions_TetN(ξ::T, η::T, ζ::T, p::Int) where T<:Real
#     nodes = generate_tet_nodes(p)

#     function N_vec(x)
#         ξ, η, ζ = x
#         λ1 = 1 - ξ - η - ζ
#         λ2 = ξ
#         λ3 = η
#         λ4 = ζ
#         [prod([
#             product_term(a, p, λ1),
#             product_term(b, p, λ2),
#             product_term(c, p, λ3),
#             product_term(d, p, λ4)
#         ]) for (a, b, c, d) in nodes]
#     end

#     x = [ξ, η, ζ]
#     N = N_vec(x)
#     J = ForwardDiff.jacobian(N_vec, x)

#     dN_dξ = J[:, 1]
#     dN_dη = J[:, 2]
#     dN_dζ = J[:, 3]

#     return N, dN_dξ, dN_dη, dN_dζ
# end


# # Dispatch functions for specific TET types
# function shape_functions(::Tet4, ξ::T, η::T, ζ::T) where T<:Real

#     return shapeFunctions_TetN(ξ, η, ζ, 1)
# end

# function shape_functions(::Tet10, ξ::T, η::T, ζ::T) where T<:Real
#     return shapeFunctions_TetN(ξ, η, ζ, 2)
# end

# function shape_functions(::Tet20, ξ::T, η::T, ζ::T) where T<:Real
#     return shapeFunctions_TetN(ξ, η, ζ, 3)
# end

# function shape_functions(::Tet35, ξ::T, η::T, ζ::T) where T<:Real
#     return shapeFunctions_TetN(ξ, η, ζ, 4)
# end

# function shape_functions(::Tet56, ξ::T, η::T, ζ::T) where T<:Real 
#     return shapeFunctions_TetN(ξ, η, ζ, 5)
# end


function generate_tet_nodes(p::Int)
    nodes = NTuple{4, Int}[]
    for a in 0:p
        for b in 0:(p - a)
            for c in 0:(p - a - b)
                d = p - a - b - c
                push!(nodes, (a, b, c, d))
            end
        end
    end
    return nodes
end

function product_term(a::Int, p::Int, λ::Real)
    if a == 0
        return one(λ)
    end
    result = one(λ)
    for k in 0:(a - 1)
        result *= (p * λ - k) / (a - k)
    end
    return result
end

function d_product_term(a::Int, p::Int, λ::Real)
    if a == 0
        return 0.0
    end
    sum = 0.0
    for j in 0:(a - 1)
        term = 1.0
        for k in 0:(a - 1)
            if k != j
                term *= (p * λ - k) / (a - k)
            end
        end
        sum += p / (a - j) * term
    end
    return sum
end

function shapeFunctions_TetN(ξ::T, η::T, ζ::T, p::Int) where T<:Real
    nodes = generate_tet_nodes(p)
    λ1 = 1 - ξ - η - ζ
    λ2 = ξ
    λ3 = η
    λ4 = ζ

    N = T[]
    dN_dξ = T[]
    dN_dη = T[]
    dN_dζ = T[]

    for (a, b, c, d) in nodes
        # Value
        P1 = product_term(a, p, λ1)
        P2 = product_term(b, p, λ2)
        P3 = product_term(c, p, λ3)
        P4 = product_term(d, p, λ4)

        Nval = P1 * P2 * P3 * P4
        push!(N, Nval)

        # Derivatives
        dP1 = d_product_term(a, p, λ1)
        dP2 = d_product_term(b, p, λ2)
        dP3 = d_product_term(c, p, λ3)
        dP4 = d_product_term(d, p, λ4)

        dλ1_dξ = -1.0
        dλ2_dξ = 1.0
        dλ3_dξ = 0.0
        dλ4_dξ = 0.0

        dλ1_dη = -1.0
        dλ2_dη = 0.0
        dλ3_dη = 1.0
        dλ4_dη = 0.0

        dλ1_dζ = -1.0
        dλ2_dζ = 0.0
        dλ3_dζ = 0.0
        dλ4_dζ = 1.0

        dNξ = (dP1 * dλ1_dξ) * P2 * P3 * P4 +
              P1 * (dP2 * dλ2_dξ) * P3 * P4 +
              P1 * P2 * (dP3 * dλ3_dξ) * P4 +
              P1 * P2 * P3 * (dP4 * dλ4_dξ)

        dNη = (dP1 * dλ1_dη) * P2 * P3 * P4 +
              P1 * (dP2 * dλ2_dη) * P3 * P4 +
              P1 * P2 * (dP3 * dλ3_dη) * P4 +
              P1 * P2 * P3 * (dP4 * dλ4_dη)

        dNζ = (dP1 * dλ1_dζ) * P2 * P3 * P4 +
              P1 * (dP2 * dλ2_dζ) * P3 * P4 +
              P1 * P2 * (dP3 * dλ3_dζ) * P4 +
              P1 * P2 * P3 * (dP4 * dλ4_dζ)

        push!(dN_dξ, dNξ)
        push!(dN_dη, dNη)
        push!(dN_dζ, dNζ)
    end

    return N, dN_dξ, dN_dη, dN_dζ
end

# Dispatch
function shape_functions(::Tet4, ξ::T, η::T, ζ::T) where T<:Real
    return shapeFunctions_TetN(ξ, η, ζ, 1)
end
function shape_functions(::Tet10, ξ::T, η::T, ζ::T) where T<:Real
    return shapeFunctions_TetN(ξ, η, ζ, 2)
end
function shape_functions(::Tet20, ξ::T, η::T, ζ::T) where T<:Real
    return shapeFunctions_TetN(ξ, η, ζ, 3)
end
function shape_functions(::Tet35, ξ::T, η::T, ζ::T) where T<:Real
    return shapeFunctions_TetN(ξ, η, ζ, 4)
end
function shape_functions(::Tet56, ξ::T, η::T, ζ::T) where T<:Real
    return shapeFunctions_TetN(ξ, η, ζ, 5)
end


# function generate_tet_nodes(p::Int)
#     p == 0 && return [(0, 0, 0, 0)]
#     nodes = NTuple{4, Int}[]
#     for d in 0:p
#         for c in 0:(p - d)
#             for b in 0:(p - d - c)
#                 a = p - d - c - b
#                 push!(nodes, (a, b, c, d))
#             end
#         end
#     end
#     return nodes
# end

# function product_term(a::Int, p::Int, λ::Real)
#     a == 0 && return 1.0
#     term = 1.0
#     for k in 0:(a-1)
#         term *= (p * λ - k) / (a - k)
#     end
#     return term
# end

# function d_product_term(a::Int, p::Int, λ::Real)
#     a == 0 && return 0.0
#     result = 0.0
#     for j in 0:(a-1)
#         term_val = 1.0
#         for k in 0:(a-1)
#             if k != j
#                 term_val *= (p * λ - k) / (a - k)
#             end
#         end
#         result += p * term_val / (a - j)
#     end
#     return result
# end

# function shapeFunctions_TetN(ξ::T, η::T, ζ::T, p::Int) where T<:Real
#     λ1 = 1 - ξ - η - ζ
#     λ2 = ξ
#     λ3 = η
#     λ4 = ζ

#     nodes = generate_tet_nodes(p)
#     n_nodes = length(nodes)
#     N = zeros(T, n_nodes)
#     dN_dξ = zeros(T, n_nodes)
#     dN_dη = zeros(T, n_nodes)
#     dN_dζ = zeros(T, n_nodes)

#     # Precompute basis functions and their derivatives
#     P = [Vector{T}(undef, 4) for _ in 1:4]
#     dP = [Vector{T}(undef, 4) for _ in 1:4]

#     for i in 1:4
#         for (idx, (a, b, c, d)) in enumerate(nodes)
#             if i == 1
#                 P[1][idx] = product_term(a, p, λ1)
#                 dP[1][idx] = d_product_term(a, p, λ1)
#             elseif i == 2
#                 P[2][idx] = product_term(b, p, λ2)
#                 dP[2][idx] = d_product_term(b, p, λ2)
#             elseif i == 3
#                 P[3][idx] = product_term(c, p, λ3)
#                 dP[3][idx] = d_product_term(c, p, λ3)
#             else
#                 P[4][idx] = product_term(d, p, λ4)
#                 dP[4][idx] = d_product_term(d, p, λ4)
#             end
#         end
#     end

#     # Compute shape functions and derivatives
#     for (idx, (a, b, c, d)) in enumerate(nodes)
#         P1 = P[1][idx]
#         P2 = P[2][idx]
#         P3 = P[3][idx]
#         P4 = P[4][idx]
        
#         dP1 = dP[1][idx]
#         dP2 = dP[2][idx]
#         dP3 = dP[3][idx]
#         dP4 = dP[4][idx]

#         # Shape function
#         N[idx] = P1 * P2 * P3 * P4

#         # Derivatives using product rule
#         dN_dξ_val = (dP1 * (-1) * P2 * P3 * P4) +
#                     (P1 * dP2 * 1 * P3 * P4) +
#                     (P1 * P2 * dP3 * 0 * P4) +
#                     (P1 * P2 * P3 * dP4 * 0)
        
#         dN_dη_val = (dP1 * (-1) * P2 * P3 * P4) +
#                     (P1 * dP2 * 0 * P3 * P4) +
#                     (P1 * P2 * dP3 * 1 * P4) +
#                     (P1 * P2 * P3 * dP4 * 0)
        
#         dN_dζ_val = (dP1 * (-1) * P2 * P3 * P4) +
#                     (P1 * dP2 * 0 * P3 * P4) +
#                     (P1 * P2 * dP3 * 0 * P4) +
#                     (P1 * P2 * P3 * dP4 * 1)
        
#         dN_dξ[idx] = dN_dξ_val
#         dN_dη[idx] = dN_dη_val
#         dN_dζ[idx] = dN_dζ_val
#     end

#     return N, dN_dξ, dN_dη, dN_dζ
# end

# # Dispatch functions remain the same
# function shape_functions(::Tet4, ξ::T, η::T, ζ::T) where T<:Real
#     return shapeFunctions_TetN(ξ, η, ζ, 1)
# end
# function shape_functions(::Tet10, ξ::T, η::T, ζ::T) where T<:Real
#     return shapeFunctions_TetN(ξ, η, ζ, 2)
# end
# function shape_functions(::Tet20, ξ::T, η::T, ζ::T) where T<:Real
#     return shapeFunctions_TetN(ξ, η, ζ, 3)
# end
# function shape_functions(::Tet35, ξ::T, η::T, ζ::T) where T<:Real
#     return shapeFunctions_TetN(ξ, η, ζ, 4)
# end
# function shape_functions(::Tet56, ξ::T, η::T, ζ::T) where T<:Real
#     return shapeFunctions_TetN(ξ, η, ζ, 5)
# end