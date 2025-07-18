# =============================================
# Core Geometric Computations
# =============================================

# ----------------------
# Jacobian Computation
# ----------------------
    """
        compute_jacobian(elem_conn, nodes, shp, dim::Int; elem_dim)

        Compute Jacobian matrices and related geometric quantities for finite element analysis.

        # Arguments
        - "elem_conn": Element connectivity (indices of nodes)
        - "nodes": Node coordinates matrix (N×dim)
        - "shp": Shape function data structure
        - "dim": Physical dimension (1, 2, or 3)
        - "elem_dim": Parametric dimension of element

        # Returns
        Vector of tuples containing for each Gauss point:
        1. Absolute determinant of Jacobian
        2. Inverse of Jacobian matrix
        3. Jacobian matrix itself

        # Notes
        - Automatically corrects element orientation
        - Handles non-square Jacobians for shells/beams
        - Throws error for nearly singular Jacobians
    """
    function compute_jacobian(elem_conn, nodes, shp, dim::Int; elem_dim)
        n_nodes = length(elem_conn)
        ngp = length(shp.weights)
        out = Vector{Tuple{Float64, Matrix{Float64}, Matrix{Float64}}}(undef, ngp)
        
        # Create local copy of coordinates
        coords = nodes[elem_conn, :]
        
        # Pre-check orientation at first Gauss point for square Jacobians
        if dim == elem_dim
            dNξ = shp.shape_ξ[1]
            dNη = (elem_dim ≥ 2) ? shp.shape_η[1] : zeros(n_nodes)
            dNζ = (elem_dim == 3) ? shp.shape_ζ[1] : zeros(n_nodes)
            
            J_test = zeros(elem_dim, dim)
            for i in 1:n_nodes
                x, y, z = coords[i, 1], coords[i, 2], (dim == 3 ? coords[i, 3] : 0.0)
                
                # Parametric derivative ξ
                J_test[1, 1] += dNξ[i] * x
                J_test[1, 2] += dNξ[i] * y
                dim == 3 && (J_test[1, 3] += dNξ[i] * z)
                
                # Parametric derivative η
                if elem_dim ≥ 2
                    J_test[2, 1] += dNη[i] * x
                    J_test[2, 2] += dNη[i] * y
                    dim == 3 && (J_test[2, 3] += dNη[i] * z)
                end
                
                # Parametric derivative ζ
                if elem_dim == 3
                    J_test[3, 1] += dNζ[i] * x
                    J_test[3, 2] += dNζ[i] * y
                    J_test[3, 3] += dNζ[i] * z
                end
            end
        
            # Flip connectivity if Jacobian is negative
            if det(J_test) < 0
                reversed_conn = reverse(elem_conn)
                coords = nodes[reversed_conn, :]
            end
        end

        # Main Gauss point loop
        for q in 1:ngp
            dNξ = shp.shape_ξ[q]
            dNη = (elem_dim >= 2 && shp.shape_η !== nothing)  ? shp.shape_η[q] : zeros(n_nodes)
            dNζ = (elem_dim == 3 && shp.shape_ζ !== nothing) ? shp.shape_ζ[q] : zeros(n_nodes)

            # Initialize Jacobian matrix
            J = zeros(elem_dim, dim)
            
            # Build Jacobian with proper derivative separation
            for i in 1:n_nodes
                node_idx = elem_conn[i]
                x = nodes[node_idx, 1]
                if dim >=  2
                    y = nodes[node_idx, 2]
                else
                    y = nothing
                end
                if dim == 3
                    z = nodes[node_idx, 3]
                else 
                    z = nothing
                end
                # ξ-derivatives (row 1)
                J[1, 1] += dNξ[i] * x  # ∂x/∂ξ
                if y !== nothing
                    J[1, 2] += dNξ[i] * y  # ∂y/∂ξ
                end
                if z !== nothing
                    J[1, 3] += dNξ[i] * z  # ∂y/∂ξ
                end
               
                # η-derivatives (row 2)
                if elem_dim >= 2
                    J[2, 1] += dNη[i] * x  # ∂x/∂η
                    J[2, 2] += dNη[i] * y  # ∂y/∂η
                    if z !== nothing
                        J[2, 3] += dNη[i] * z  # ∂y/∂ξ
                    end  
                end
                
                # ζ-derivatives (row 3)
                if elem_dim == 3
                    J[3, 1] += dNζ[i] * x  # ∂x/∂ζ
                    J[3, 2] += dNζ[i] * y  # ∂y/∂ζ
                    J[3, 3] += dNζ[i] * z  # ∂z/∂ζ
                end
            end
           
            # Compute determinant and inverse
            if dim == elem_dim
                detJ = det(J)
                abs_detJ = abs(detJ)
                invJ = inv(J)
            else
                # Pseudo-determinant for non-square Jacobians
                abs_detJ = sqrt(abs(det(J * J')))
                invJ = pinv(J)  # Moore-Penrose pseudoinverse
            end
        
            abs_detJ < 1e-12 && error("Jacobian nearly singular at Gauss point $q")
            out[q] = (abs_detJ, invJ, J)
        end
        return out
    end

# ----------------------
# Shape Data Precomputation
# ----------------------
    """
        shape_data(element_type, int_order, dim)

        Precompute shape functions and derivatives at Gauss points.

        # Arguments
        - "element_type": Element type (e.g., Tri6, Quad8)
        - "int_order": Integration order
        - "dim": Physical dimension

        # Returns
        Named tuple containing:
        - "N": Shape functions at Gauss points
        - "shape_ξ": ξ-derivatives
        - "shape_η": η-derivatives (2D/3D)
        - "shape_ζ": ζ-derivatives (3D)
        - "weights": Integration weights

        # Notes
        - Automatically selects appropriate Gauss quadrature
        - Returns nothing for unused derivatives
    """
    function shape_data(element_type, int_order, dim)
        gauss_points, gauss_weights = integration_rule(element_type, int_order)
        shape_data = [(shape_functions(element_type, ξ...)..., w)
                    for (ξ, w) in zip(gauss_points, gauss_weights)]
        
        return (    
            N       = [sd[1] for sd in shape_data], 
            shape_ξ = [sd[2] for sd in shape_data],
            shape_η = dim >= 2 ? [sd[3] for sd in shape_data] : nothing,
            shape_ζ = dim == 3 ? [sd[4] for sd in shape_data] : nothing,
            weights = [sd[end] for sd in shape_data] 
        )
    end

# ----------------------
# Jacobian Data Precomputation
# ----------------------
    """
        jacobian_data(connectivity, nodes, shape_data, dim; elem_dim=dim)

        Precompute Jacobian data for all elements.

        # Arguments
        - "connectivity": Matrix of element connectivities
        - "nodes": Node coordinates
        - "shape_data": Precomputed shape function data
        - "dim": Physical dimension
        - "elem_dim": Parametric dimension

        # Returns
        Vector of Jacobian data (one entry per element)

        # Notes
        - Wrapper around "compute_jacobian" for multiple elements
        - Maintains consistent interface for all element types
    """
    function jacobian_data(connectivity, nodes, shape_data, dim; elem_dim=dim)
        Ne = size(connectivity, 1)
        return [compute_jacobian(connectivity[i, :], nodes, shape_data,dim; elem_dim) for i in 1:Ne]
    end

# =============================================
# B-Matrix Construction
# =============================================
    # ----------------------
    # Strain B-Matrix
    # ----------------------
        """
            build_strain_B(elem_jac, gauss_data, dim, n_nodes)

            Construct strain-displacement (B) matrices.

            # Arguments
            - "elem_jac": Precomputed Jacobian data
            - "gauss_data": Shape function data
            - "dim": Physical dimension
            - "n_nodes": Number of nodes per element

            # Returns
            Vector of B-matrices (one per Gauss point)

            # Notes
            - Handles 1D, 2D, and 3D cases
            - Returns properly ordered strain components
        """
        function build_strain_B(elem_jac, gauss_data, dim, n_nodes)
            n_gauss = length(gauss_data.weights)
            strain_comp = dim == 1 ? 1 : dim == 2 ? 3 : 6
            B = [zeros(strain_comp, dim*n_nodes) for _ in 1:n_gauss]
            fill_strain_B!(B, elem_jac, gauss_data, dim, n_nodes)
            return B
        end

    # ----------------------
    # Gradient B-Matrix
    # ----------------------

        """
            build_gradient_B(elem_jac, gauss_data, dim, n_nodes)

            Construct gradient operator matrices.

            # Arguments
            - "elem_jac": Precomputed Jacobian data
            - "gauss_data": Shape function data
            - "dim": Physical dimension
            - "n_nodes": Number of nodes per element

            # Returns
            Vector of gradient matrices (one per Gauss point)

            # Notes
            - Used for scalar field gradients
            - Consistent with strain B-matrix formulation
        """
        function build_gradient_B(elem_jac, gauss_data, dim, n_nodes)
            n_gauss = length(gauss_data.weights)
            B = [zeros(dim, n_nodes) for _ in 1:n_gauss]
            fill_gradient_B!(B, elem_jac, gauss_data, dim, n_nodes)
            return B
        end

    # ----------------------
    # Strain B-Matrix Computation
    # ----------------------
        """
            fill_strain_B!(B, elem_jac, gauss_data, dim, n_nodes)

            Fill preallocated strain B-matrices.

            # Arguments
            - "B": Preallocated matrix storage
            - "elem_jac": Jacobian data
            - "gauss_data": Shape function data
            - "dim": Physical dimension
            - "n_nodes": Nodes per element

            # Notes
            - In-place operation for performance
            - Implements proper Voigt ordering
        """
        function fill_strain_B!(B, elem_jac, gauss_data, dim, n_nodes)
            # Setup parametric derivative iterator
            iter = dim == 1 ? zip(gauss_data.shape_ξ) :
                dim == 2 ? zip(gauss_data.shape_ξ, gauss_data.shape_η) :
                zip(gauss_data.shape_ξ, gauss_data.shape_η, gauss_data.shape_ζ)

            for (qp, dN_tuple) in enumerate(iter)
                _, invJ = elem_jac[qp]
                for i in 1:n_nodes
                    # Compute physical gradient
                    ref_grad = dim == 1 ? SVector(dN_tuple[1][i]) :
                            dim == 2 ? SVector(dN_tuple[1][i], dN_tuple[2][i]) :
                            SVector(dN_tuple[1][i], dN_tuple[2][i], dN_tuple[3][i])
                    ∇N_physical = invJ * ref_grad

                    # Fill strain B-matrix components
                    if dim == 1
                        B[qp][1, i] = ∇N_physical[1]
                    elseif dim == 2
                        B[qp][1, 2i-1] = ∇N_physical[1]  # ∂u/∂x
                        B[qp][2, 2i  ] = ∇N_physical[2]  # ∂v/∂y
                        B[qp][3, 2i-1] = ∇N_physical[2]  # ∂u/∂y
                        B[qp][3, 2i  ] = ∇N_physical[1]  # ∂v/∂x
                    else  # dim == 3
                        B[qp][1, 3i-2] = ∇N_physical[1]  # ∂u/∂x
                        B[qp][2, 3i-1] = ∇N_physical[2]  # ∂v/∂y
                        B[qp][3, 3i  ] = ∇N_physical[3]  # ∂w/∂z
                        B[qp][4, 3i-1] = ∇N_physical[3]  # ∂v/∂z
                        B[qp][4, 3i  ] = ∇N_physical[2]  # ∂w/∂y
                        B[qp][5, 3i-2] = ∇N_physical[3]  # ∂u/∂z
                        B[qp][5, 3i  ] = ∇N_physical[1]  # ∂w/∂x
                        B[qp][6, 3i-2] = ∇N_physical[2]  # ∂u/∂y
                        B[qp][6, 3i-1] = ∇N_physical[1]  # ∂v/∂x
                    end
                end
            end
        end

    # ----------------------
    # Gradient B-Matrix Computation
    # ----------------------

        """
            fill_gradient_B!(B, elem_jac, gauss_data, dim, n_nodes)

            Fill preallocated gradient matrices.

            # Arguments
            - "B": Preallocated matrix storage
            - "elem_jac": Jacobian data
            - "gauss_data": Shape function data
            - "dim": Physical dimension
            - "n_nodes": Nodes per element

            # Notes
            - Optimized for scalar field problems
            - Reuses Jacobian computations
        """
        function fill_gradient_B!(B, elem_jac, gauss_data, dim, n_nodes)
            # Setup parametric derivative iterator
            iter = dim == 1 ? zip(gauss_data.shape_ξ) :
                dim == 2 ? zip(gauss_data.shape_ξ, gauss_data.shape_η) :
                zip(gauss_data.shape_ξ, gauss_data.shape_η, gauss_data.shape_ζ)

            for (qp, dN_tuple) in enumerate(iter)
                _, invJ = elem_jac[qp]
                for i in 1:n_nodes
                    # Compute physical gradient
                    ref_grad = dim == 1 ? SVector(dN_tuple[1][i]) :
                            dim == 2 ? SVector(dN_tuple[1][i], dN_tuple[2][i]) :
                            SVector(dN_tuple[1][i], dN_tuple[2][i], dN_tuple[3][i])
                    ∇N_physical = invJ * ref_grad
                    
                    # Store gradient components
                    B[qp][:, i] = ∇N_physical[1:dim]
                end
            end
        end

    # ----------------------
    # B-Matrix Assembly
    # ----------------------

        """
            build_B_matrices(nodes, connectivity, material, gauss_data, jacobian_data)

            Assemble all required B-matrices for a material.

            # Arguments
            - "nodes": Node coordinates
            - "connectivity": Element connectivity
            - "material": Material definition
            - "gauss_data": Shape function data
            - "jacobian_data": Precomputed Jacobians

            # Returns
            Vector of dictionaries containing B-matrices:
            - Keys: B-matrix types (:strain, :gradient, etc.)
            - Values: Vectors of matrices (per Gauss point)

            # Notes
            - Thread-parallel implementation
            - Handles multiple B-matrix types
            - Supports out-of-plane symmetry
        """
        function build_B_matrices(
            nodes::Union{Matrix{Float64}, Vector{Float64}},
            connectivity::Matrix{Int},
            material::Material,
            gauss_data,
            jacobian_data)

            dim = material.dim
            n_elem = size(connectivity, 1)
            n_nodes_per_elem = size(connectivity, 2)
            B_data = Vector{Dict{Symbol, Vector{Matrix{Float64}}}}(undef, n_elem)
            
            Threads.@threads for e in 1:n_elem
                elem_jacobian = jacobian_data[e] 
                B_dict = Dict{Symbol, Vector{Matrix{Float64}}}()

                # Build B-matrices for each required type
                for B_type in material.B_types
                    B = if B_type == :strain
                        build_strain_B(elem_jacobian, gauss_data, dim, n_nodes_per_elem)
                    elseif B_type ∈ (:strain_rate, :gradient, :electric_field, 
                                    :pressure, :temperature_gradient, :fluid_velocity)
                        build_gradient_B(elem_jacobian, gauss_data, dim, n_nodes_per_elem)
                    else
                        error("Unsupported B-type: $B_type")
                    end

                    # Special handling for out-of-plane symmetry
                    if material.symmetry == :out_of_plane 
                        for i in eachindex(B)
                            B[i] = vcat(B[i], zeros(1, size(B[i], 2)))
                        end
                    end
                
                    B_dict[B_type] = B
                end

                B_data[e] = B_dict
            end
            return B_data
        end