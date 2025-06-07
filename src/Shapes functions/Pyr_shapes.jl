

# =============================================
# Pyramid Shape Functions (3D Coordinates)
# =============================================



# Pyramid local coordinates: 
# ξ, η ∈ [-1,1] (base quadrilateral directions)
# ζ ∈ [0,1] (height from base to apex)

# Helper function for apex blending
ζ_blend(ζ) = 1 - ζ

# =============================================
# Pyr5 Implementation (Linear)
# =============================================

function shape_functions(::Pyr5, ξ::T, η::T, ζ::T) where T<:Real
    α = ζ_blend(ζ)
    N = [
        (1 - ξ)*(1 - η)*α/4,  # Base node 1 (-1,-1,0)
        (1 + ξ)*(1 - η)*α/4,  # Base node 2 (1,-1,0)
        (1 + ξ)*(1 + η)*α/4,  # Base node 3 (1,1,0)
        (1 - ξ)*(1 + η)*α/4,  # Base node 4 (-1,1,0)
        ζ                       # Apex node 5 (0,0,1)
    ]
    
    ∇N = [
        # ∂N/∂ξ        ∂N/∂η         ∂N/∂ζ
        -(1 - η)*α/4,  -(1 - ξ)*α/4,  -(1 - ξ)*(1 - η)/4,  # N1
         (1 - η)*α/4,  -(1 + ξ)*α/4,  -(1 + ξ)*(1 - η)/4,  # N2
         (1 + η)*α/4,   (1 + ξ)*α/4,  -(1 + ξ)*(1 + η)/4,  # N3
        -(1 + η)*α/4,   (1 - ξ)*α/4,  -(1 - ξ)*(1 + η)/4,  # N4
        0.0,            0.0,            1.0                 # N5
    ]
    return N, ∇N[:,1], ∇N[:,2], ∇N[:,3]
end

# =============================================
# Pyr14 Implementation (Quadratic)
# =============================================

function shape_functions(::Pyr14, ξ::T, η::T, ζ::T) where T<:Real
    ζb = 1 - ζ  # Blending term (1 - ζ)
    
    N = Vector{T}(undef, 14)
    
   
    N = zeros(14)
    
    # Base corner nodes (1-4)
    N[1] = 0.25*(1 - ξ)*(1 - η)*(1 - ζ) - 0.25*(1 - ξ^2)*(1 - η)*(1 - ζ) - 0.25*(1 - ξ)*(1 - η^2)*(1 - ζ) - 0.25*(1 - ξ^2)*(1 - η^2)*(1 - ζ)
    N[2] = 0.25*(1 + ξ)*(1 - η)*(1 - ζ) - 0.25*(1 - ξ^2)*(1 - η)*(1 - ζ) - 0.25*(1 + ξ)*(1 - η^2)*(1 - ζ) - 0.25*(1 - ξ^2)*(1 - η^2)*(1 - ζ)
    N[3] = 0.25*(1 + ξ)*(1 + η)*(1 - ζ) - 0.25*(1 - ξ^2)*(1 + η)*(1 - ζ) - 0.25*(1 + ξ)*(1 - η^2)*(1 - ζ) - 0.25*(1 - ξ^2)*(1 - η^2)*(1 - ζ)
    N[4] = 0.25*(1 - ξ)*(1 + η)*(1 - ζ) - 0.25*(1 - ξ^2)*(1 + η)*(1 - ζ) - 0.25*(1 - ξ)*(1 - η^2)*(1 - ζ) - 0.25*(1 - ξ^2)*(1 - η^2)*(1 - ζ)
    
    # Apex node (5)
    N[5] = ζ
    
    # Base mid-edge nodes (6-9)
    N[6] = 0.5*(1 - ξ^2)*(1 - η)*(1 - ζ)
    N[7] = 0.5*(1 + ξ)*(1 - η^2)*(1 - ζ)
    N[8] = 0.5*(1 - ξ^2)*(1 + η)*(1 - ζ)
    N[9] = 0.5*(1 - ξ)*(1 - η^2)*(1 - ζ)
    
    # ηertical mid-edge nodes (10-13)
    N[10] = 4.0*ζ*(1 - ζ)*(1 + ξ)*(1 + η)
    N[11] = 4.0*ζ*(1 - ζ)*(1 - ξ)*(1 + η)
    N[12] = 4.0*ζ*(1 - ζ)*(1 - ξ)*(1 - η)
    N[13] = 4.0*ζ*(1 - ζ)*(1 + ξ)*(1 - η)
    
    # Base center node (14)
    N[14] = (1 - ξ^2)*(1 - η^2)*(1 - ζ)
    


    ∇N = zeros(T, 14, 3)
    
    ζp = 1 + ξ + η + ζ
    ∇N[1,:] = [
        ((1-η)*(ζb*(1 + ξ + η + ζ) - (1-ξ)*ζb) - (1-ξ)*(1-η)*ζb)/4,
        ((1-ξ)*(ζb*(1 + ξ + η + ζ) - (1-η)*ζb) - (1-ξ)*(1-η)*ζb)/4,
        (1-ξ)*(1-η)*(ζp - ζb)/4
    ]
    
    ζp = 1 - ξ + η + ζ
    ∇N[2,:] = [
        ((-1+η)*(ζb*(1 - ξ + η + ζ) + (1+ξ)*(1-η)*ζb))/4,
        ((1+ξ)*(ζb*(1 - ξ + η + ζ) - (1-η)*ζb))/4,
        (1+ξ)*(1-η)*(ζp - ζb)/4
    ]
    
    ζp = 1 - ξ - η + ζ
    ∇N[3,:] = [
        ((-1-η)*(ζb*(1 - ξ - η + ζ) + (1+ξ)*(1+η)*ζb)/4),
        ((-1-ξ)*(ζb*(1 - ξ - η + ζ) + (1+ξ)*(1+η)*ζb)/4),
        (1+ξ)*(1+η)*(ζp - ζb)/4
    ]
    
    ζp = 1 + ξ - η + ζ
    ∇N[4,:] = [
        ((1+η)*(ζb*(1 + ξ - η + ζ) - (1-ξ)*(1+η)*ζb)/4),
        ((-1+ξ)*(ζb*(1 + ξ - η + ζ) - (1+η)*ζb)/4),
        (1-ξ)*(1+η)*(ζp - ζb)/4
    ]
    
    ∇N[5,:] = [0.0, 0.0, 4ζ - 1]
    
    ∇N[6,:] = [
        -ξ*(1-η)*ζb,
        -(1-ξ^2)*ζb/2,
        -(1-ξ^2)*(1-η)/2
    ]
    
    ∇N[7,:] = [
        (1-η^2)*ζb/2,
        -η*(1+ξ)*ζb,
        -(1+ξ)*(1-η^2)/2
    ]
    
    ∇N[8,:] = [
        -ξ*(1+η)*ζb,
        (1-ξ^2)*ζb/2,
        -(1-ξ^2)*(1+η)/2
    ]
    
    ∇N[9,:] = [
        -(1-η^2)*ζb/2,
        -η*(1-ξ)*ζb,
        -(1-ξ)*(1-η^2)/2
    ]
    
    ∇N[10,:] = [
        -(1-η)*ζ*ζb/2,
        -(1-ξ)*ζ*ζb/2,
        (1-ξ)*(1-η)*(1-2ζ)/2
    ]
    
    ∇N[11,:] = [
        (1-η)*ζ*ζb/2,
        -(1+ξ)*ζ*ζb/2,
        (1+ξ)*(1-η)*(1-2ζ)/2
    ]
    
    ∇N[12,:] = [
        (1+η)*ζ*ζb/2,
        (1+ξ)*ζ*ζb/2,
        (1+ξ)*(1+η)*(1-2ζ)/2
    ]
    
    ∇N[13,:] = [
        -(1+η)*ζ*ζb/2,
        (1-ξ)*ζ*ζb/2,
        (1-ξ)*(1+η)*(1-2ζ)/2
    ]
    
    ∇N[14,:] = [
        -2ξ*(1-η^2)*ζb,
        -2η*(1-ξ^2)*ζb,
        -(1-ξ^2)*(1-η^2)
    ]

    return N, ∇N[:,1], ∇N[:,2], ∇N[:,3]
end




# function shapeFunctions_PyramidN(ξ::Real, η::Real, ζ::Real, p::Int)
#     # Step 1: Map to unit pyramid coordinates
#     x̂ = (ξ + 1) / 2
#     ŷ = (η + 1) / 2
#     ẑ = (ζ + 1) / 2

#     # Step 2: Generate base 1D nodes
#     nodes = collect(range(0.0, 1.0, length = p + 1))

#     # Step 3: Build Lagrange basis for 1D
#     function lagrange_basis(i, nodes, x)
#         prod([(x - nodes[j]) / (nodes[i] - nodes[j]) for j in 1:length(nodes) if j != i])
#     end
#     function lagrange_deriv(i, nodes, x)
#         sum([
#             prod([(x - nodes[k]) / (nodes[i] - nodes[k]) for k in 1:length(nodes) if k != i && k != j]) /
#             (nodes[i] - nodes[j]) for j in 1:length(nodes) if j != i
#         ])
#     end

#     # Step 4: Pyramid node coordinates in [0,1]^3
#     coords = []
#     for k in 0:p
#         h = k / p
#         nxy = p - k
#         if nxy == 0
#             push!(coords, (0.5, 0.5, h))
#             continue
#         end
#         for j in 0:nxy
#             for i in 0:nxy
#                 push!(coords, (i / nxy, j / nxy, h))
#             end
#         end
#     end

#     # Step 5: Shape functions
#     N = Float64[]
#     dN_dξ = Float64[]
#     dN_dη = Float64[]
#     dN_dζ = Float64[]

#     for (x, y, z) in coords
#         i = findfirst(isequal(x), nodes)
#         j = findfirst(isequal(y), nodes)
#         k = findfirst(isequal(z), nodes)

#         if i === nothing || j === nothing || k === nothing
#             error("Invalid local coordinate ($x, $y, $z) for order $p")
#         end

#         ℓx = lagrange_basis(i, nodes, x̂)
#         ℓy = lagrange_basis(j, nodes, ŷ)
#         ℓz = lagrange_basis(k, nodes, ẑ)

#         dℓx = lagrange_deriv(i, nodes, x̂) * 0.5
#         dℓy = lagrange_deriv(j, nodes, ŷ) * 0.5
#         dℓz = lagrange_deriv(k, nodes, ẑ) * 0.5

#         push!(N, ℓx * ℓy * ℓz)
#         push!(dN_dξ, dℓx * ℓy * ℓz)
#         push!(dN_dη, ℓx * dℓy * ℓz)
#         push!(dN_dζ, ℓx * ℓy * dℓz)
#     end

#     return N, dN_dξ, dN_dη, dN_dζ
# end

# function shape_functions(::Pyr5, ξ::T, η::T, ζ::T) where T<:Real
#     return shapeFunctions_PyramidN(ξ, η, ζ, 1)
# end

# function shape_functions(::Pyr14, ξ::T, η::T, ζ::T) where T<:Real
#     return shapeFunctions_PyramidN(ξ, η, ζ, 2)
# end

# function shape_functions(::Pyr29, ξ::T, η::T, ζ::T) where T<:Real
#     return shapeFunctions_PyramidN(ξ, η, ζ, 3)
# end

# function shape_functions(::Pyr50, ξ::T, η::T, ζ::T) where T<:Real
#     return shapeFunctions_PyramidN(ξ, η, ζ, 4)
# end

# function shape_functions(::Pyr77, ξ::T, η::T, ζ::T) where T<:Real
#     return shapeFunctions_PyramidN(ξ, η, ζ, 5)
# end
