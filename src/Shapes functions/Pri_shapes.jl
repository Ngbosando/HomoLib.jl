# =============================================
# Prismatic Element Structs
# =============================================

# =============================================
# Pri6 Implementation (Linear)
# =============================================


function shape_functions(::Pri6, ξ::T, η::T, ζ::T) where T<:Real
    # Shape functions for linear prism (6 nodes)
    # Bottom triangle (ζ = -1)
    N = T[
        (1 - ξ - η)*(1 - ζ)/2,   # Node 1
        ξ*(1 - ζ)/2,             # Node 2
        η*(1 - ζ)/2,             # Node 3
        (1 - ξ - η)*(1 + ζ)/2,   # Node 4
        ξ*(1 + ζ)/2,             # Node 5
        η*(1 + ζ)/2              # Node 6
    ]
    
    ∇N = Matrix{T}(undef, 6, 3)
    # ∂N/∂ξ, ∂N/∂η, ∂N/∂ζ derivatives
 
    ∇N[1,:] = [-(1 - ζ)*0.5, -(1 - ζ)*0.5, -(1 - ξ - η)*0.5]
    ∇N[2,:] = [(1 - ζ)*0.5, 0.0, -ξ*0.5]
    ∇N[3,:] = [0.0, (1 - ζ)*0.5, -η*0.5]
    ∇N[4,:] = [-(1 + ζ)*0.5, -(1 + ζ)*0.5, (1 - ξ - η)*0.5]
    ∇N[5,:] = [(1 + ζ)*0.5, 0.0, ξ*0.5]
    ∇N[6,:] = [0.0, (1 + ζ)*0.5, η*0.5]
    
    return N, ∇N[:,1], ∇N[:,2], ∇N[:,3]
end

# =============================================
# Pri18 Implementation (Quadratic)
# =============================================

function shape_functions(::Pri18, ξ::T, η::T, ζ::T) where T<:Real
    node_coords = [
        (0.0, 0.0, -1.0), (1.0, 0.0, -1.0), (0.0, 1.0, -1.0),
        (0.5, 0.0, -1.0), (0.5, 0.5, -1.0), (0.0, 0.5, -1.0),
        (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0),
        (0.5, 0.0, 0.0), (0.5, 0.5, 0.0), (0.0, 0.5, 0.0),
        (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (0.0, 1.0, 1.0),
        (0.5, 0.0, 1.0), (0.5, 0.5, 1.0), (0.0, 0.5, 1.0)
    ]
    
    N = Vector(undef, 18)
    ∇N = Matrix(undef, 18, 3)
    
    function tri_basis(ξ::T, η::T)
        L1 = 1 - ξ - η
        L2 = ξ
        L3 = η
        return [
            L1*(2L1 - 1),
            L2*(2L2 - 1),
            L3*(2L3 - 1),
            4L1*L2,
            4L2*L3,
            4L3*L1
        ]
    end
    
    function tri_deriv(ξ::T, η::T)
        dN_dξ = [
            -4*(1 - ξ - η) + 1,
            4ξ - 1,
            0.0,
            4*(1 - 2ξ - η),
            -4η,
            4η
        ]
        dN_dη = [
            -4*(1 - ξ - η) + 1,
            0.0,
            4η - 1,
            -4ξ,
            4ξ,
            4*(1 - ξ - 2η)
        ]
        return dN_dξ, dN_dη
    end
    
    line_basis = [
        ζ*(ζ - 1)/2,
        1 - ζ^2,
        ζ*(ζ + 1)/2
    ]
    
    line_deriv = [
        (2ζ - 1)/2,
        -2ζ,
        (2ζ + 1)/2
    ]
    
    for (i, (ξi, ηi, ζi)) in enumerate(node_coords)
        tri_idx = ((i - 1) % 6) + 1
        layer = (i - 1) ÷ 6 + 1
        
        tri_vals = tri_basis(ξ, η)
        dTri_dξ, dTri_dη = tri_deriv(ξ, η)
        
        line_val = line_basis[layer]
        dLine_dζ = line_deriv[layer]
        
        N[i] = tri_vals[tri_idx] * line_val
        
        ∇N[i,1] = dTri_dξ[tri_idx] * line_val
        ∇N[i,2] = dTri_dη[tri_idx] * line_val
        ∇N[i,3] = tri_vals[tri_idx] * dLine_dζ
    end
    
    return N, ∇N[:,1], ∇N[:,2], ∇N[:,3]
end


# ---------------------------------------------
# Lagrange basis in 1D
# function lagrange_shape_functions(nodes::Vector{T}, ξ::T) where T<:Real
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

# ---------------------------------------------
# GMSH-style ordering
function gmsh_ordering_prism(p::Int)::Vector{Int}
    return collect(1:(p*(p+1)÷2) * p)  # Placeholder, user should provide actual gmsh ordering
end

# ---------------------------------------------
# Triangular Bernstein basis
function triangle_shape_derivatives(p::Int, ξ::T, η::T) where T
    L1 = 1 - ξ - η
    L2 = ξ
    L3 = η
    n = div((p + 1) * (p + 2), 2)
    N     = Vector{T}(undef, n)
    dN_dξ = similar(N)
    dN_dη = similar(N)

    idx = 1
    @inbounds for i in 0:p
        for j in 0:(p - i)
            k = p - i - j
            c = factorial(p) / (factorial(i) * factorial(j) * factorial(k))
            lk = L1^k; li = L2^i; lj = L3^j
            N[idx] = c * lk * li * lj
            dN_dξ[idx] = c * (k * L1^(k - 1) * li * lj + i * lk * L2^(i - 1) * lj)
            dN_dη[idx] = c * (k * L1^(k - 1) * li * lj + j * lk * li * L3^(j - 1))
            idx += 1
        end
    end

    return N, dN_dξ, dN_dη
end

# ---------------------------------------------
# Main shape function builder
function shapeFunctions_PrismN(ξ::T, η::T, ζ::T, p::Int) where T<:Real
    nodes_1d = collect(range(-1.0, 1.0; length=p+1))
    Lζ, dLζ = lagrange_shape_functions(nodes_1d, ζ)

    N_tri, dN_ξ, dN_η = triangle_shape_derivatives(p, ξ, η)

   n_tri = length(N_tri)
    total = (p + 1) * n_tri
    N     = Vector{T}(undef, total)
    dN_dξ = similar(N)
    dN_dη = similar(N)
    dN_dζ = similar(N)

    idx = 1
    @inbounds for k in 1:(p + 1)
        lk = Lζ[k]; dlk = dLζ[k]
        for i in 1:n_tri
            nt = N_tri[i]
            N[idx]     = nt * lk
            dN_dξ[idx] = dN_ξ[i] * lk
            dN_dη[idx] = dN_η[i] * lk
            dN_dζ[idx] = nt * dlk
            idx += 1
        end
    end

    reorder = gmsh_ordering_prism(p + 1)
    return N[reorder], dN_dξ[reorder], dN_dη[reorder], dN_dζ[reorder]
end

# ---------------------------------------------
# Dispatches
# shape_functions(::Pri6, ξ::T, η::T, ζ::T) where T<:Real = shapeFunctions_PrismN(ξ, η, ζ, 1)
# shape_functions(::Pri18, ξ::T, η::T, ζ::T) where T<:Real = shapeFunctions_PrismN(ξ, η, ζ, 2)
shape_functions(::Pri36, ξ::T, η::T, ζ::T) where T<:Real = shapeFunctions_PrismN(ξ, η, ζ, 3)
shape_functions(::Pri56, ξ::T, η::T, ζ::T) where T<:Real = shapeFunctions_PrismN(ξ, η, ζ, 4)
shape_functions(::Pri78, ξ::T, η::T, ζ::T) where T<:Real = shapeFunctions_PrismN(ξ, η, ζ, 5)

# ---------------------------------------------
# Reference coordinates generator
function reference_coords_prism(p::Int)
    coords = []
    for k in 0:p
        ζ = k / p
        for i in 0:p, j in 0:(p - i)
            ξ = i / p
            η = j / p
            push!(coords, (ξ, η, ζ))
        end
    end
    return coords
end
