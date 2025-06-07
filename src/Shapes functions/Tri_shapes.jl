# =============================================
# Triangular Element Structs
# =============================================


# =============================================
# Shape Function Implementations
# =============================================

# Helper functions for barycentric coordinates
ζ(ξ, η) = 1 - ξ - η

# # =============================================
# # Tri3 Implementation (Linear)
# # =============================================

function shape_functions(::Tri3, ξ::T, η::T) where T<:Real
    N = [ζ(ξ,η), ξ, η]
    ∇N = [-1 -1; 1 0; 0 1]
    return N, ∇N[:,1], ∇N[:,2]
end

# # =============================================
# # Tri6 Implementation (Quadratic)
# # =============================================

function shape_functions(::Tri6, ξ::T, η::T) where T<:Real
    z = ζ(ξ,η)
    N = [
        z*(2z - 1),
        ξ*(2ξ - 1),
        η*(2η - 1),
        4z*ξ,
        4ξ*η,
        4η*z
    ]
    ∇N = [
        (-3 + 4 * ξ + 4 * η)  (-3 + 4 * ξ + 4 * η)
        (4ξ - 1)       0.0
        0.0            (4η - 1)
        (4z - 4ξ)      (-4ξ)
        (4η)           (4ξ)
        (-4η)          (4z - 4η)
    ]
    return N, ∇N[:,1], ∇N[:,2]
end

# =============================================
# Tri10 Implementation (Cubic)
# =============================================

function shape_functions(::Tri10, ξ::T, η::T) where T<:Real
    z = ζ(ξ, η)
    N = [
        z*(3z - 1)*(3z - 2)/2,        # Vertex 1
        ξ*(3ξ - 1)*(3ξ - 2)/2,        # Vertex 2
        η*(3η - 1)*(3η - 2)/2,        # Vertex 3
        9//2*z*ξ*(3z - 1),            # Edge 1-2 node 1
        9//2*z*ξ*(3ξ - 1),            # Edge 1-2 node 2
        9//2*ξ*η*(3ξ - 1),            # Edge 2-3 node 1
        9//2*ξ*η*(3η - 1),            # Edge 2-3 node 2
        9//2*η*z*(3η - 1),            # Edge 3-1 node 1
        9//2*η*z*(3z - 1),            # Edge 3-1 node 2
        27*z*ξ*η                       # Internal node
    ]
    
    # Derivatives of each shape function (∂N/∂ξ, ∂N/∂η)
    ∇N = [
        # Vertex 1: N1
        (-27z^2 + 18z - 2)/2  (-27z^2 + 18z - 2)/2
        # Vertex 2: N2
        (27ξ^2 - 18ξ + 2)/2   0
        # Vertex 3: N3
        0                     (27η^2 - 18η + 2)/2
        # Edge 1-2 node 1: N4
        (9//2)*(3z^2 - z + ξ - 6z*ξ)  (9//2)*ξ*(1 - 6z)
        # Edge 1-2 node 2: N5
        (9//2)*(-3ξ^2 + ξ + z*(6ξ - 1))  -9//2*ξ*(3ξ - 1)
        # Edge 2-3 node 1: N6
        (9//2)*η*(6ξ - 1)               (9//2)*ξ*(3ξ - 1)
        # Edge 2-3 node 2: N7
        (9//2)*η*(3η - 1)               (9//2)*ξ*(6η - 1)
        # Edge 3-1 node 1: N8
        -9//2*η*(3η - 1)                (9//2)*(-3η^2 + η + z*(6η - 1))
        # Edge 3-1 node 2: N9
        (9//2)*η*(1 - 6z)               (9//2)*(3z^2 - z + η - 6z*η)
        # Internal node: N10
        27η*(z - ξ)                     27ξ*(z - η)
    ]
    return N, ∇N[:,1], ∇N[:,2]
end

# =============================================
# Tri15 Implementation (Quintic)
# =============================================

function shape_functions(::Tri15, ξ::T, η::T) where T<:Real
    z = ξ + η - 1
    N = Vector{T}(undef, 15)
    
    # Définition des fonctions de forme
    N[1]  = ((ξ + η - 1) * (2z + 1) * (4z + 1) * (4ξ + 4η - 1)) / 3
    N[2]  = (ξ * (2ξ - 1) * (4ξ - 3) * (4ξ - 1)) / 3
    N[3]  = (η * (2η - 1) * (4η - 3) * (4η - 1)) / 3
    N[4]  = -(16 * ξ * z * (2z + 1) * (4z + 1)) / 3
    N[5]  =  4 * ξ * z * (4ξ - 1) * (4z + 1)
    N[6]  = 16 * ξ * -z * (2ξ - 1) * (4ξ - 1) / 3
    N[7]  = (16 * η * ξ * (2ξ - 1) * (4ξ - 1)) / 3
    N[8]  = 4 * η * (4η - 1) * ξ * (4ξ - 1)
    N[9]  = (16 * η * ξ * (2η - 1) * (4η - 1)) / 3
    N[10] = -(16 * η * z * (2η - 1) * (4η - 1)) / 3
    N[11] = 4 * η * z * (4η - 1) * (4z + 1)
    N[12] = -(16 * η * z * (2z + 1) * (4z + 1)) / 3
    N[13] = 32 * η * ξ * z * (4z + 1)
    N[14] = -32 * η * ξ * z * (4ξ - 1)
    N[15] = -32 * η * ξ * z * (4η - 1)
  
    ∇N = zeros(15, 2)
    # Derivative of N1
    dN1_dz = ((2z + 1)*(4z + 3)*(4z + 1) + z*2*(4z + 3)*(4z + 1) + z*(2z + 1)*(4*(4z + 3) + 4*(4z + 1))) / 3
    ∇N[1, 1] = dN1_dz * (1)
    ∇N[1, 2] = dN1_dz * (1)
    
    # Derivative of N2
    ∇N[2, 1] = ((2ξ - 1)*(4ξ - 3)*(4ξ - 1) + ξ*2*(4ξ - 3)*(4ξ - 1) + ξ*(2ξ - 1)*(4*(4ξ - 3) + 4*(4ξ - 1))) / 3
    ∇N[2, 2] = 0
    
    # Derivative of N3
    ∇N[3, 2] = ((2η - 1)*(4η - 3)*(4η - 1) + η*2*(4η - 3)*(4η - 1) + η*(2η - 1)*(4*(4η - 3) + 4*(4η - 1))) / 3
    ∇N[3, 1] = 0
    

    # Derivative of N4
    ∇N[4, 1] = -(16/3) * (z*(2z + 1)*(4z + 1) + ξ*( (2z + 1)*(4z + 1) + z*(2*(4z + 1) + (2z + 1)*4) ))
    ∇N[4, 2] = -(16/3) * ξ * ( (2z + 1)*(4z + 1) + z*(2*(4z + 1) + (2z + 1)*4) )
    
    # Derivative of N5
    ∇N[5, 1] = 4* ( z*(4ξ - 1)*(4z + 1) + ξ*(4ξ - 1)*(4z + 1) + 4*ξ*z*((4ξ - 1)+(4z + 1)) )
    ∇N[5, 2] = 4 * ξ * ( (4ξ - 1)*(4z + 1) + z*(4ξ - 1)*4 )
    
    # Derivative of N6
    ∇N[6, 1] = -(16/3) * (z*(2ξ - 1)*(4ξ - 1) + ξ*( (2ξ - 1)*(4ξ - 1) + z*(2*(4ξ - 1) + (2ξ - 1)*4) ))
    ∇N[6, 2] = -(16/3) * ξ * ( (2ξ - 1)*(4ξ - 1) )
        

    # Derivative of N7
    ∇N[7, 1] = (16/3) * (η*(2ξ - 1)*(4ξ - 1) + ξ*η*(2*(4ξ - 1) + (2ξ - 1)*4))
    ∇N[7, 2] = (16/3) * ξ*(2ξ - 1)*(4ξ - 1)
    
    # Derivative of N8
    ∇N[8, 1] = 4 * η*(4η - 1) * ((4ξ - 1) + 4ξ)
    ∇N[8, 2] = 4 * (ξ*(4ξ - 1)*(4η - 1) + ξ*η*4*(4ξ - 1))
    
    # Derivative of N9
    ∇N[9, 1] = (16/3) * η*(2η - 1)*(4η - 1)
    ∇N[9, 2] = (16/3) * (ξ*(2η - 1)*(4η - 1) + ξ*η*(2*(4η - 1) + (2η - 1)*4))
    
    # Derivative of N10
    ∇N[10, 1] = -(16/3) * η * ( (2η - 1)*(4η - 1) )
    ∇N[10, 2] = -(16/3) * ( z*(2η - 1)*(4η - 1) + η*( (2η - 1)*(4η - 1) + z*(2*(4η - 1) + (2η - 1)*4) ) )
    
    # Derivative of N11
    ∇N[11, 1] = 4 * η * (4η - 1)*((4z + 1) + z*4 )
    ∇N[11, 2] = 4 * ( z*(4η - 1)*(4z + 1) + η*( (4η - 1)*(4z + 1) + z*(4*(4z + 1) + (4η - 1)*4) ) )
    
    # Derivative of N12
    ∇N[12, 1] = -(16/3) * η * ( (2z+1)*(4z+1) + z*(2*(4z+1) + (2z+1)*4) )
    ∇N[12, 2] = -(16/3) * ( z*(2z + 1)*(4z + 1) + η*( (2z + 1)*(4z + 1) + z*(2*(4z + 1) + (2z + 1)*4) ) )
    
    # Derivative of N13
    ∇N[13, 1] = 32 * η* ( z*(4z + 1) + ξ*( z*4 + (4z + 1)*1 ) )
    ∇N[13, 2] = 32 * ( ξ*z*(4z + 1) + ξ*η*( z*4 + (4z + 1)*1 ) )
    
    # Derivative of N14
    ∇N[14, 1] = -32 *  η* (z*(4ξ - 1) + ξ*( z*4 + (4ξ - 1)*1 ) )
    ∇N[14, 2] = -32 * ξ * ( z*(4ξ - 1) + η*(  (4ξ - 1)*1 ) )
    
    # Derivative of N15
    ∇N[15, 1] = -32 * η * ( z*(4η - 1) + ξ*(4η - 1)*1 ) 
    ∇N[15, 2] = -32 * ( ξ*z*(4η - 1) + ξ*η*4*z + ξ*η*(4η - 1)  )
        
    
    return N, ∇N[:,1], ∇N[:,2]
end


# =============================================
# Tri21 Implementation (Quadratic)
# =============================================

ξ_val_21 = [
    (0, 0),  # Node 1
    (1.0, 0.0),  # Node 2
    (0.0, 1.0),  # Node 3
    (0.2, 0.0),  # Node 4
    (0.4, 0.0),  # Node 5
    (0.6, 0.0),  # Node 6
    (0.8, 0.0),  # Node 7
    (0.8, 0.2),  # Node 8
    (0.6, 0.4),  # Node 9
    (0.4, 0.6),  # Node 10
    (0.2, 0.8),  # Node 11
    (0.0, 0.8),  # Node 12
    (0.0, 0.6),  # Node 13
    (0.0, 0.4),  # Node 14
    (0.0, 0.2),  # Node 15
    (0.2, 0.2),  # Node 16
    (0.4, 0.2),  # Node 17
    (0.2, 0.4),  # Node 18
    (0.6, 0.2),  # Node 19
    (0.2, 0.6),  # Node 20
    (0.4, 0.4)   # Node 21
]


# function shape_functions(::Tri21, ξ::T, η::T) where T<:Real
#     h = 1e-8
#     function bary_poly(i, j, k, ξ, η)
#         if i < 0 || j < 0 || k < 0
#             return 0.0
#         end
#         L1 = 1.0 - ξ - η
#         L2 = ξ
#         L3 = η
#         return factorial(5) / (factorial(i)*factorial(j)*factorial(k)) * (L1^i) * (L2^j) * (L3^k)
#     end

#     idxs = [
#         (5, 0, 0), (4, 1, 0), (3, 2, 0), (2, 3, 0), (1, 4, 0), (0, 5, 0),
#         (4, 0, 1), (3, 1, 1), (2, 2, 1), (1, 3, 1), (0, 4, 1),
#         (3, 0, 2), (2, 1, 2), (1, 2, 2), (0, 3, 2),
#         (2, 0, 3), (1, 1, 3), (0, 2, 3),
#         (1, 0, 4), (0, 1, 4),
#         (0, 0, 5)
#     ]

#     N = [bary_poly(i, j, k, ξ, η) for (i, j, k) in idxs]

#     Nξ_plus = [bary_poly(i, j, k, ξ+h, η) for (i, j, k) in idxs]
#     Nξ_minus = [bary_poly(i, j, k, ξ-h, η) for (i, j, k) in idxs]
#     dN_dξ = [(Nξ_plus[i] - Nξ_minus[i]) / (2h) for i in 1:21]

#     Nη_plus = [bary_poly(i, j, k, ξ, η+h) for (i, j, k) in idxs]
#     Nη_minus = [bary_poly(i, j, k, ξ, η-h) for (i, j, k) in idxs]
#     dN_dη = [(Nη_plus[i] - Nη_minus[i]) / (2h) for i in 1:21]

#     return N, dN_dξ, dN_dη
# end

# function shapeFunctions_Tri21_dual(ξ::Real, η::Real)
#     input = [ξ, η]
#     f(ξη) = shapeFunctions_Tri21(ξη[1], ξη[2])[1]  # only N, no grad
#     val, grad = ForwardDiff.value.(f(input)), ForwardDiff.gradient(f, input)
#     return val, grad[1:21:end], grad[2:21:end]  # reshape gradients per variable
# end
function shape_functions(::Tri21, ξ::T, η::T) where T<:Real
    h = 1e-8
    function bary_poly(i, j, k, ξ, η)
        if i < 0 || j < 0 || k < 0
            return 0.0
        end
        L1 = 1.0 - ξ - η
        L2 = ξ
        L3 = η
        return factorial(5) / (factorial(i)*factorial(j)*factorial(k)) * (L1^i) * (L2^j) * (L3^k)
    end

    idxs = [
        (5, 0, 0), (4, 1, 0), (3, 2, 0), (2, 3, 0), (1, 4, 0), (0, 5, 0),
        (4, 0, 1), (3, 1, 1), (2, 2, 1), (1, 3, 1), (0, 4, 1),
        (3, 0, 2), (2, 1, 2), (1, 2, 2), (0, 3, 2),
        (2, 0, 3), (1, 1, 3), (0, 2, 3),
        (1, 0, 4), (0, 1, 4),
        (0, 0, 5)
    ]

    N = [bary_poly(i, j, k, ξ, η) for (i, j, k) in idxs]

    Nξ_plus = [bary_poly(i, j, k, ξ+h, η) for (i, j, k) in idxs]
    Nξ_minus = [bary_poly(i, j, k, ξ-h, η) for (i, j, k) in idxs]
    dN_dξ = [(Nξ_plus[i] - Nξ_minus[i]) / (2h) for i in 1:21]

    Nη_plus = [bary_poly(i, j, k, ξ, η+h) for (i, j, k) in idxs]
    Nη_minus = [bary_poly(i, j, k, ξ, η-h) for (i, j, k) in idxs]
    dN_dη = [(Nη_plus[i] - Nη_minus[i]) / (2h) for i in 1:21]

    return N, dN_dξ, dN_dη
end

function shapeFunctions_Tri21_dual(ξ::Real, η::Real)
    function all_N(uv)
        L1 = 1 - uv[1] - uv[2]
        L2 = uv[1]
        L3 = uv[2]
        idxs = [
            (5, 0, 0), (4, 1, 0), (3, 2, 0), (2, 3, 0), (1, 4, 0), (0, 5, 0),
            (4, 0, 1), (3, 1, 1), (2, 2, 1), (1, 3, 1), (0, 4, 1),
            (3, 0, 2), (2, 1, 2), (1, 2, 2), (0, 3, 2),
            (2, 0, 3), (1, 1, 3), (0, 2, 3),
            (1, 0, 4), (0, 1, 4),
            (0, 0, 5)
        ]
        N = similar(uv, 21)
        for i in 1:21
            a, b, c = idxs[i]
            N[i] = factorial(5) / (factorial(a)*factorial(b)*factorial(c)) * L1^a * L2^b * L3^c
        end
        return N
    end

    uv = ForwardDiff.Dual.(fill(0.0, 2), (1.0, 0.0), (0.0, 1.0))
    uv[1] = ξ
    uv[2] = η
    N = all_N(uv)
    dN_dξ = ForwardDiff.derivative.(x -> all_N([x, η]), Ref(ξ))
    dN_dη = ForwardDiff.derivative.(y -> all_N([ξ, y]), Ref(η))
    return N, dN_dξ, dN_dη
end


# function lagrange_shape_functions_Tri(nodes::Vector{T}, ξ::T) where T<:Real
#     n = length(nodes)
#     L = zeros(T, n)
#     dL = zeros(T, n)

#     for i in 1:n
#         Li = one(T)
#         dLi = zero(T)
#         for j in 1:n
#             if j != i
#                 Li *= (ξ - nodes[j]) / (nodes[i] - nodes[j])
#                 term = one(T) / (nodes[i] - nodes[j])
#                 prod = one(T)
#                 for k in 1:n
#                     if k != i && k != j
#                         prod *= (ξ - nodes[k]) / (nodes[i] - nodes[k])
#                     end
#                 end
#                 dLi += term * prod
#             end
#         end
#         L[i] = Li
#         dL[i] = dLi
#     end
#     return L, dL
# end

# function multinomial(i, j, k)
#     return factorial(i + j + k) ÷ (factorial(i) * factorial(j) * factorial(k))
# end

# function gmsh_ordering_tri(p::Int)::Vector{Int}
#     idx_map = Dict{Tuple{Int, Int, Int}, Int}()
#     node_index = 1
#     for k in 0:p
#         for j in 0:(p - k)
#             i = p - j - k
#             idx_map[(i, j, k)] = node_index
#             node_index += 1
#         end
#     end

#     ordering = Int[]

#     # corner nodes
#     push!(ordering, idx_map[(p, 0, 0)])  # vertex 1 (L1)
#     push!(ordering, idx_map[(0, p, 0)])  # vertex 2 (L2)
#     push!(ordering, idx_map[(0, 0, p)])  # vertex 3 (L3)

#     # edge nodes
#     if p > 1
#         for i in 1:p-1
#             push!(ordering, idx_map[(p - i, i, 0)])  # edge 1: (1-2)
#             push!(ordering, idx_map[(0, p - i, i)])  # edge 2: (2-3)
#             push!(ordering, idx_map[(i, 0, p - i)])  # edge 3: (3-1)
#         end
#     end

#     # face interior nodes
#     if p > 2
#         for k in 1:p-2
#             for j in 1:p - k - 1
#                 i = p - j - k
#                 push!(ordering, idx_map[(i, j, k)])
#             end
#         end
#     end

#     return ordering
# end

# function shapeFunctions_TriN(ξ::T, η::T, p::Int) where T<:Real
#     L1 = 1 - ξ - η
#     L2 = ξ
#     L3 = η

#     N = T[]
#     dN_dξ = T[]
#     dN_dη = T[]
#     coords = Tuple{T,T}[]

#     for k in 0:p
#         for j in 0:(p - k)
#             i = p - j - k
#             coeff = multinomial(i, j, k)
#             Ni = coeff * L1^i * L2^j * L3^k

#             dNi_dL1 = coeff * i * L1^(i-1) * L2^j * L3^k
#             dNi_dL2 = coeff * j * L1^i * L2^(j-1) * L3^k
#             dNi_dL3 = coeff * k * L1^i * L2^j * L3^(k-1)

#             dNi_dξ = -dNi_dL1 + dNi_dL2
#             dNi_dη = -dNi_dL1 + dNi_dL3

#             push!(N, Ni)
#             push!(dN_dξ, dNi_dξ)
#             push!(dN_dη, dNi_dη)
#             push!(coords, (L2, L3))
#         end
#     end

#     reorder = gmsh_ordering_tri(p)

#     # println("--- TRI Node Coordinates (ξ, η) ---")
#     # for (i, (x, y)) in enumerate(coords)
#     #     println(rpad(i, 2), ": (", x, ", ", y, ")")
#     # end

#     # println("--- Gmsh TRI Ordering ---")
#     # for (i, idx) in enumerate(reorder)
#     #     println("Gmsh Index ", i+1, " → Lex Idx ", idx, ": ", coords[idx])
#     # end

#     return N[reorder], dN_dξ[reorder], dN_dη[reorder]
# end



# Dispatch functions for specific TRI types
# function shape_functions(::Tri3, ξ::T, η::T) where T<:Real
#     return shapeFunctions_TriN(ξ, η, 2)
# end

# function shape_functions(::Tri6, ξ::T, η::T) where T<:Real
#     return shapeFunctions_TriN(ξ, η, 3)
# end

# function shape_functions(::Tri10, ξ::T, η::T) where T<:Real
#     return shapeFunctions_TriN(ξ, η, 4)
# end

# function shape_functions(::Tri15, ξ::T, η::T) where T<:Real
#     return shapeFunctions_TriN(ξ, η, 5)
# end

# function shape_functions(::Tri21, ξ::T, η::T) where T<:Real
#     N, dN_dξ, dN_dη = shapeFunctions_TriN_with_derivatives(ξ, η, 5)
#     return N, dN_dξ, dN_dη
# end



