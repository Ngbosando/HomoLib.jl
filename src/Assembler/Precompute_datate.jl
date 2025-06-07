# =============================================
# Geometric jacobian computation
# =============================================

function compute_jacobian(elem_connect, nodes, dN_parametric)
    dim = size(nodes, 2)
    n_nodes = length(elem_connect)
    num_gauss = length(dN_parametric.weights)

    invJ_detJ = Vector{Tuple{Float64, AbstractMatrix{<:Real}}}(undef, num_gauss)

    for e in 1:num_gauss
        dN_ξ = dN_parametric.shape_ξ[e]
        dN_η = dim >= 2 ? dN_parametric.shape_η[e] : zeros(length(dN_ξ))
        dN_ζ = dim == 3 ? dN_parametric.shape_ζ[e] : zeros(length(dN_ξ))

        J = if dim == 1
            sum(nodes[elem_connect[i], 1] * dN_ξ[i] for i in 1:n_nodes)
        elseif dim == 2
            J = @SMatrix [0.0 0.0; 0.0 0.0]
            for i in 1:n_nodes
                node_idx = elem_connect[i]
                x = nodes[node_idx, 1]
                y = nodes[node_idx, 2]
                J += @SMatrix [dN_ξ[i] * x  dN_ξ[i] * y;
                               dN_η[i] * x  dN_η[i] * y]
            end
            J
        else  # dim == 3
            J = zeros(3, 3)
            for i in 1:n_nodes
                x, y, z = nodes[elem_connect[i], 1:3]
                J[1,1] += dN_ξ[i] * x
                J[1,2] += dN_η[i] * x
                J[1,3] += dN_ζ[i] * x
                J[2,1] += dN_ξ[i] * y
                J[2,2] += dN_η[i] * y
                J[2,3] += dN_ζ[i] * y
                J[3,1] += dN_ξ[i] * z
                J[3,2] += dN_η[i] * z
                J[3,3] += dN_ζ[i] * z
            end
            J
        end

        detJ = det(J)
        invJ = inv(J)
        invJ_detJ[e] = (detJ, invJ)
    end
    return invJ_detJ
end


# =============================================
# Precompute all derivatives data
# =============================================

function build_B_matrices(nodes::Matrix{Float64},
    connectivity::Matrix{Int},
    material::Material,
    gauss_data,
    jacobian_data)
    dim = material.dim
    n_elem = size(connectivity, 1)
    n_gauss = length(gauss_data.weights)
    n_nodes_per_elem = size(connectivity, 2)

    B_data = Vector{Dict{Symbol, Vector{Matrix{Float64}}}}(undef, n_elem)
    
    Threads.@threads for e in 1:n_elem
        elem_jacobian = jacobian_data[e] 
        B_dict = Dict{Symbol, Vector{Matrix{Float64}}}()

        # Generic B-matrix construction
        for B_type in material.B_types
            B = if B_type == :strain
                build_strain_B(elem_jacobian, gauss_data, dim, n_nodes_per_elem)
            elseif B_type ∈ (:gradient, :electric_field, :pressure, :temperature_gradient)
                build_gradient_B(elem_jacobian, gauss_data, dim, n_nodes_per_elem)
            else
                error("Unsupported B-type: $B_type")
            end
            if material.symmetry == :out_of_plane
                # Append a row of zeros to each B matrix to capture ε₃₃
                lb = length(B)
                for i in 1: lb
                    B[i] = vcat(B[i], zeros(1, size(B[i], 2)))  # 1×(2n)
                end
            end
           
            B_dict[B_type] = B
        end

        B_data[e] = B_dict
    end
    return B_data
end
# =============================================
# B construction for strain 
# =============================================

function build_strain_B(elem_jac, gauss_data, dim, n_nodes)
    n_gauss = length(gauss_data.weights)
    strain_comp = dim == 1 ? 1 : dim == 2 ? 3 : 6
    B = [zeros(strain_comp, dim*n_nodes) for _ in 1:n_gauss]
    fill_strain_B!(B, elem_jac, gauss_data, dim, n_nodes)
    return B
end
# =============================================
# B construction for gradient
# =============================================

function build_gradient_B(elem_jac, gauss_data, dim, n_nodes)
    n_gauss = length(gauss_data.weights)
    B = [zeros(dim, n_nodes) for _ in 1:n_gauss]
 
    fill_gradient_B!(B, elem_jac, gauss_data, dim, n_nodes)
    return B
end

# =============================================
#  B computatopn for gradient 
# =============================================

function fill_gradient_B!(B, elem_jac, gauss_data, dim, n_nodes)
    # Prepare iterator based on dim
    if dim == 1
        iter = zip(gauss_data.shape_ξ)
    elseif dim == 2
        iter = zip(gauss_data.shape_ξ, gauss_data.shape_η)
    elseif dim == 3
        iter = zip(gauss_data.shape_ξ, gauss_data.shape_η, gauss_data.shape_ζ)
    else
        error("Unsupported dimension: $dim")
    end

    for (qp, dN_tuple) in enumerate(iter)
        _, invJ = elem_jac[qp]
        for i in 1:n_nodes
            ref_grad = dim == 1 ? SVector(dN_tuple[1][i]) :
                       dim == 2 ? SVector(dN_tuple[1][i], dN_tuple[2][i]) :
                       SVector(dN_tuple[1][i], dN_tuple[2][i], dN_tuple[3][i])
            ∇N_physical = invJ * ref_grad
            B[qp][:, i] = ∇N_physical[1:dim]
        end
    end
end

# =============================================
# B computatopn for strain 
# =============================================

function fill_strain_B!(B, elem_jac, gauss_data, dim, n_nodes)
    if dim == 1
        iter = zip(gauss_data.shape_ξ)
    elseif dim == 2
        iter = zip(gauss_data.shape_ξ, gauss_data.shape_η)
    elseif dim == 3
        iter = zip(gauss_data.shape_ξ, gauss_data.shape_η, gauss_data.shape_ζ)
    else
        error("Unsupported dimension: $dim")
    end

    for (qp, dN_tuple) in enumerate(iter)
        _, invJ = elem_jac[qp]
        for i in 1:n_nodes
            ref_grad = dim == 1 ? SVector(dN_tuple[1][i]) :
                       dim == 2 ? SVector(dN_tuple[1][i], dN_tuple[2][i]) :
                       SVector(dN_tuple[1][i], dN_tuple[2][i], dN_tuple[3][i])
            ∇N_physical = invJ * ref_grad

            if dim == 1
                B[qp][1, i] = ∇N_physical[1]
            elseif dim == 2
                B[qp][1, 2i-1] = ∇N_physical[1]    # ∂u/∂x
                B[qp][2, 2i  ] = ∇N_physical[2]    # ∂v/∂y
                B[qp][3, 2i-1] = ∇N_physical[2]    # ∂u/∂y
                B[qp][3, 2i  ] = ∇N_physical[1]    # ∂v/∂x
            else # dim == 3
                B[qp][1, 3i-2] = ∇N_physical[1]    # ∂u/∂x
                B[qp][2, 3i-1] = ∇N_physical[2]    # ∂v/∂y
                B[qp][3, 3i  ] = ∇N_physical[3]    # ∂w/∂z
                B[qp][4, 3i-1] = ∇N_physical[3]    # ∂v/∂z
                B[qp][4, 3i  ] = ∇N_physical[2]    # ∂w/∂y
                B[qp][5, 3i-2] = ∇N_physical[3]    # ∂u/∂z
                B[qp][5, 3i  ] = ∇N_physical[1]    # ∂w/∂x
                B[qp][6, 3i-2] = ∇N_physical[2]    # ∂u/∂y
                B[qp][6, 3i-1] = ∇N_physical[1]    # ∂v/∂x
            end
        end
    end
end


# =============================================
# Precompute all shapes data
# =============================================

function shape_data(element_type, int_order, dim)
    gauss_data = let order = int_order

        gauss_points, gauss_weights = integration_rule(element_type, int_order)
      
        shape_data = [(shape_functions(element_type, ξ...)..., w)
                       for (ξ, w) in zip(gauss_points, gauss_weights)]
        (    
            N       = [sd[1] for sd in shape_data], 
            shape_ξ = [sd[2] for sd in shape_data],  # dN_dξ derivatives
            shape_η = dim >= 2 ? [sd[3] for sd in shape_data] : nothing,  # dN_dη derivatives
            shape_ζ = dim == 3 ? [sd[4] for sd in shape_data] : nothing,  # dN_dζ derivatives
            weights = [sd[dim + 2] for sd in shape_data] 
           
        )
    end
    return gauss_data
end

# =============================================
# Precompute all jacobian data
# =============================================

function jacobian_data(connectivity, nodes, shape_data)
    Ne = size(connectivity, 1)
    jacobian_data = [compute_jacobian(connectivity[i, :], nodes, shape_data) for i in 1:Ne]
    return jacobian_data
end
