include("Precompute_data.jl")
include("material_model.jl")
include("resolve_material_type.jl")
include("transform_boundary.jl")

# =============================================
#  Element Stiffness computation
# =============================================
    # Dispatch function remains the same
        function compute_element_stiffness(material::Material, B_dict::Dict, jacobian_data, gauss_data)
            T = resolve_material_type(material)
            return compute_element_stiffness(Val(T), material, B_dict, jacobian_data, gauss_data)
        end

    # ==================== ELASTIC ==================== #
        function compute_element_stiffness(::Val{:elastic}, material::Material, B_dict, jacobian_data, gauss_data)
            C = material.tensors[1]
            B = B_dict[:strain]
            n_dofs = size(B[1], 2)
           
            K = zeros(n_dofs, n_dofs)
            CB = zeros(size(C, 1), n_dofs)
            
            for qp in 1:length(gauss_data.weights)
                detJ, _ = jacobian_data[qp]
                scaling = abs(detJ) * gauss_data.weights[qp]
               
                B_qp = B[qp]
                mul!(CB, C, B_qp)
                mul!(K, B_qp', CB, scaling, 1.0)
               
            end
            
            return K
        end

    # ==================== PIEZOELECTRIC ==================== #
        function compute_element_stiffness(::Val{:piezoelectric}, material::Material, B_dict, jacobian_data, gauss_data)
            C, ϵ, e = material.tensors
            B_u = B_dict[:strain]
            B_ϕ = B_dict[:electric_field]
            
            n_dofs_u = size(B_u[1], 2)
            n_dofs_ϕ = size(B_ϕ[1], 2)

            # Initialize submatrices
            K_uu = zeros(n_dofs_u, n_dofs_u)
            K_ϕu = zeros(n_dofs_ϕ, n_dofs_u)
            K_ϕϕ = zeros(n_dofs_ϕ, n_dofs_ϕ)

            # Temporary buffers
            CB = zeros(size(C, 1), n_dofs_u)
            eB = zeros(size(e, 1), n_dofs_u)
            ϵB = zeros(size(ϵ, 1), n_dofs_ϕ)
            
            for qp in eachindex(gauss_data.weights)
                detJ, _ = jacobian_data[qp]
                scaling = abs(detJ) * gauss_data.weights[qp]
                B_u_qp = B_u[qp]
                B_ϕ_qp = B_ϕ[qp]
                
                # Mechanical stiffness: K_uu += B_u' * C * B_u
                mul!(CB, C, B_u_qp)
                mul!(K_uu, B_u_qp', CB, scaling, 1.0)
                
                # Coupling stiffness: K_ϕu += B_ϕ' * e * B_u
                mul!(eB, e, B_u_qp)
                mul!(K_ϕu, B_ϕ_qp', eB, scaling, 1.0)
                
                # Dielectric stiffness: K_ϕϕ += B_ϕ' * ϵ * B_ϕ
                mul!(ϵB, ϵ, B_ϕ_qp)
                mul!(K_ϕϕ, B_ϕ_qp', ϵB, scaling, 1.0)
            end

            # Form field-wise block matrix
            K_block = [K_uu  -K_ϕu';
                    K_ϕu  K_ϕϕ]

            # Get number of nodes from B-matrix dimensions
            dofs_per_field = [get_dofs_per_node(Material(material.type, material.dim, material.symmetry, material.properties, Vector{AbstractMatrix{Real}}(), [bt],nothing)) for bt in material.B_types]
                
            n_nodes = size(B_u[1], 2) ÷ dofs_per_field[1]  # Displacement DOFs per node

            # Permute to node-wise ordering
            K_element = permute_element_matrix(K_block, material, n_nodes)
            return K_element
        end

    # ==================== POROELASTIC  ==================== #
        function compute_element_stiffness(::Val{:poroelastic}, material::Material, B_dict, jacobian_data, gauss_data)
            C = material.tensors[1]
            B = B_dict[:strain]
            n_dofs = size(B[1], 2)
           
            K = zeros(n_dofs, n_dofs)
            CB = zeros(size(C, 1), n_dofs)
             
            for qp in 1:length(gauss_data.weights)
              
                detJ, _ = jacobian_data[qp]
                scaling =  abs(detJ) * gauss_data.weights[qp]
               
                B_qp = B[qp]
                mul!(CB, C, B_qp)
                mul!(K, B_qp', CB, scaling, 1.0)
               
            end
            
            return K
        end

    # ==================== THERMAL ==================== #
        function compute_element_stiffness(::Val{:thermal}, material::Material, B_dict, jacobian_data, gauss_data)
            κ = material.tensors[1]
            B = B_dict[:temperature_gradient]
            n_dofs = size(B[1], 2)
            
            K = zeros(n_dofs, n_dofs)
            κB = zeros(size(κ, 1), n_dofs)
            
            for qp in 1:length(gauss_data.weights)
                detJ, _ = jacobian_data[qp]
              
                scaling =  abs(detJ) * gauss_data.weights[qp]
                B_qp = B[qp]
                mul!(κB, κ, B_qp)
                mul!(K, B_qp', κB, scaling, 1.0)
            end
            
            return K
        end

    # ==================== VISCOELASTIC ==================== #
        function compute_element_stiffness(::Val{:viscoelastic}, material::Material, B_dict, jacobian_data, gauss_data)
             # Place_Holder need to be implement latter  
        end
    # ====================  Stokes Flow ==================== #    
        function compute_element_stiffness(::Val{:stokes}, material, B_dict, jac_data, gauss_data)
            μ = material.properties[:μ]
            B_v = B_dict[:velocity]
            B_p = B_dict[:pressure]
            
            n_dofs_v = size(B_v[1], 2)
            n_dofs_p = size(B_p[1], 2)

            K_vv = zeros(n_dofs_v, n_dofs_v)  # Viscous term
            K_vp = zeros(n_dofs_v, n_dofs_p)  # Divergence term

            for qp in eachindex(gauss_data.weights)
                detJ, _ = jac_data[qp]
                scaling = abs(detJ) * gauss_data.weights[qp]
                B_v_qp = B_v[qp]
                B_p_qp = B_p[qp]

                K_vv .+= μ * scaling * (B_v_qp' * B_v_qp)  # μ∫∇v:∇w
                K_vp .+= scaling * (B_v_qp' * B_p_qp)      # ∫p(∇·v)
            end

            return [K_vv  K_vp; 
                    K_vp' zeros(n_dofs_p, n_dofs_p)]  # Saddle-point system
        end       
    # =============================================
    # Helper Functions for  Element Stiffness  
    # =============================================    

        function permute_element_matrix(K_block, material::Material, n_nodes::Int)
            fields = material.B_types
            length(fields) == 1 && return K_block  # No permutation needed for single-field
            
            # Get DOFs per field (e.g., 2 for strain, 1 for electric)
            dofs_per_field = [get_dofs_per_node(Material(material.type, material.dim, material.symmetry, material.properties, Vector{AbstractMatrix{Real}}(), [bt],nothing)) for bt in fields]
            dofs_per_node = sum(dofs_per_field)
            total_dofs = n_nodes * dofs_per_node
            @assert size(K_block, 1) == total_dofs "Matrix size mismatch"

            # Build permutation vector p: maps node-wise index -> block-wise index
            p = zeros(Int, total_dofs)
            offset_field = 0
            for (j, bt) in enumerate(fields)
                dof_field = dofs_per_field[j]
                for node in 0:(n_nodes-1)
                    offset_in_node = sum(dofs_per_field[1:j-1])
                    for k in 1:dof_field
                        desired_idx = node * dofs_per_node + offset_in_node + k
                        block_idx = offset_field + node * dof_field + k
                        p[desired_idx] = block_idx
                    end
                end
                offset_field += n_nodes * dof_field
            end
            return K_block[p, p]  # Reorder matrix
        end
# =============================================
#  Element Mass Matrix computation
# =============================================
    # Dispatch function
        function compute_element_mass(material::Material, gauss_data, jacobian_data)
            T = resolve_material_type(material.type)
            return compute_element_mass(Val(T), material, gauss_data, jacobian_data)
        end

    # ==================== ELASTIC ==================== #
        function compute_element_mass(::Val{:elastic}, material::Material, gauss_data, jacobian_data)
            ρ = material.mass_properties[:ρ]  # Density
            dim = material.dim
            n_nodes = length(gauss_data.N[1])
            n_dofs = n_nodes * dim
            
            M = zeros(n_dofs, n_dofs)
            
            for qp in 1:length(gauss_data.weights)
                detJ, _ = jacobian_data[qp]
                scaling = abs(detJ) * gauss_data.weights[qp] * ρ
                N = gauss_data.N[qp]
                
                for i in 1:n_nodes, j in 1:n_nodes
                    Nij = N[i] * N[j] * scaling
                    for d in 1:dim
                        i_idx = (i-1)*dim + d
                        j_idx = (j-1)*dim + d
                        M[i_idx, j_idx] += Nij
                    end
                end
            end
            
            return M
        end

    # ==================== PIEZOELECTRIC ==================== #
        function compute_element_mass(::Val{:piezoelectric}, material::Material, gauss_data, jacobian_data)
            ρ = material.mass_properties[:ρ]  # Density
            dim = material.dim
            n_nodes = length(gauss_data.N[1])
            n_dofs_u = n_nodes * dim
            n_dofs_ϕ = n_nodes  # Electric potential DOFs
            
            # Only displacement DOFs have mass
            M_uu = zeros(n_dofs_u, n_dofs_u)
            
            for qp in 1:length(gauss_data.weights)
                detJ, _ = jacobian_data[qp]
                scaling = abs(detJ) * gauss_data.weights[qp] * ρ
                N = gauss_data.N[qp]
                
                for i in 1:n_nodes, j in 1:n_nodes
                    Nij = N[i] * N[j] * scaling
                    for d in 1:dim
                        i_idx = (i-1)*dim + d
                        j_idx = (j-1)*dim + d
                        M_uu[i_idx, j_idx] += Nij
                    end
                end
            end
            
            # Full mass matrix (displacement + electric)
            M = zeros(n_dofs_u + n_dofs_ϕ, n_dofs_u + n_dofs_ϕ)
            M[1:n_dofs_u, 1:n_dofs_u] = M_uu
            
            return M
        end

    # ==================== POROELASTIC ==================== #
        function compute_element_mass(::Val{:poroelastic}, material::Material, gauss_data, jacobian_data)
            ρ = material.mass_properties[:ρ]  # Density 
            dim = material.dim
            n_nodes = length(gauss_data.N[1])
            n_dofs_u = n_nodes * dim
            n_dofs_p = n_nodes  # Pressure DOFs
            
            # Only displacement DOFs have mass
            M_uu = zeros(n_dofs_u, n_dofs_u)
            
            for qp in 1:length(gauss_data.weights)
                detJ, _ = jacobian_data[qp]
                scaling = abs(detJ) * gauss_data.weights[qp] * ρ
                N = gauss_data.N[qp]
                
                for i in 1:n_nodes, j in 1:n_nodes
                    Nij = N[i] * N[j] * scaling
                    for d in 1:dim
                        i_idx = (i-1)*dim + d
                        j_idx = (j-1)*dim + d
                        M_uu[i_idx, j_idx] += Nij
                    end
                end
            end
            
            # Full mass matrix (displacement + pressure)
            M = zeros(n_dofs_u + n_dofs_p, n_dofs_u + n_dofs_p)
            M[1:n_dofs_u, 1:n_dofs_u] = M_uu
            
            return M
        end

    # ==================== THERMAL ==================== #
        function compute_element_mass(::Val{:thermal}, material::Material, gauss_data, jacobian_data)
            # Thermal mass matrix: ∫ ρ c N' N dΩ
            ρ = material.mass_properties[:ρ]  # Density
            c = material.mass_properties[:c]  # heat capacity
            n_nodes = length(gauss_data.N[1])
            
            M = zeros(n_nodes, n_nodes)
            
            for qp in 1:length(gauss_data.weights)
                detJ, _ = jacobian_data[qp]
                scaling = abs(detJ) * gauss_data.weights[qp] * ρ * c
                N = gauss_data.N[qp]
                
                for i in 1:n_nodes, j in 1:n_nodes
                    M[i, j] += N[i] * N[j] * scaling
                end
            end
            
            return M
        end

    # ==================== VISCOELASTIC ==================== #
        function compute_element_mass(::Val{:viscoelastic}, material::Material, gauss_data, jacobian_data)
            ρ = material.mass_properties[:ρ]  # Density
            dim = material.dim
            n_nodes = length(gauss_data.N[1])
            n_dofs = n_nodes * dim
            
            M = zeros(n_dofs, n_dofs)
            
            for qp in 1:length(gauss_data.weights)
                detJ, _ = jacobian_data[qp]
                scaling = abs(detJ) * gauss_data.weights[qp] * ρ
                N = gauss_data.N[qp]
                
                for i in 1:n_nodes, j in 1:n_nodes
                    Nij = N[i] * N[j] * scaling
                    for d in 1:dim
                        i_idx = (i-1)*dim + d
                        j_idx = (j-1)*dim + d
                        M[i_idx, j_idx] += Nij
                    end
                end
            end
            
            return M
        end        
# =============================================
#  Functions for Global Matrix computation
# =============================================
    # =============================================
    # Helper Functions for  Global dof Assembler 
    # =============================================
        function assemble_global_dofs(connectivity::Union{Vector{Int} ,Matrix{Int}}, material::Material)
            n_elements, nodes_per_element = size(connectivity)
            dofs_per_node = get_dofs_per_node(material)
            global_dofs = Matrix{Int}(undef, n_elements, nodes_per_element * dofs_per_node)
          
            field_offsets = cumsum([get_dofs_per_node(Material(material.type, material.dim, material.symmetry, material.properties, Vector{AbstractMatrix{Real}}(), [bt],nothing)) 
                            for bt in material.B_types])
            
            for e in 1:n_elements
                elem_dofs = Int[]
                for node in connectivity[e, :]
                    base_dof = (node - 1) * dofs_per_node
                    for (i, bt) in enumerate(material.B_types)
                        n_dofs = get_dofs_per_node(Material(material.type, material.dim, material.symmetry, material.properties, Vector{AbstractMatrix{Real}}(), [bt],nothing))
                        append!(elem_dofs, base_dof .+ (1:n_dofs) .+ (i > 1 ? field_offsets[i-1] : 0))
                    end
                end
            
                global_dofs[e, :] = elem_dofs
            end
            
            return global_dofs
        end

        function should_apply_internal_force(material::Material)::Bool
            material.symmetry == :out_of_plane || haskey(material.properties, :α)
        end
        function get_dofs_per_node(material::Material)
            dof_counts = Dict(
                :strain => material.dim,
                :gradient => 1,
                :electric_field => 1,
                :pressure => 1,
                :temperature_gradient => 1
            )
            
            sum(dof_counts[bt] for bt in material.B_types)
        end

    # =============================================
    #  Global Matrix Assembler
    # =============================================
        function assemble_global_matrix(
            connectivity::Matrix{Int},
            nodes::Union{Matrix{Float64}, Vector{Float64}},
            element_type::Symbol,
            integration_order::Int,
            material::Union{Material, Vector{Material}},
            dim::Int,
            connect_elem_phase::Union{Nothing, Vector{Int}},
            Geometric_Data;
            NodalForces = nothing,
            BoundaryFace = nothing,
            PointForces = nothing)

            # ==================== PHASE HANDLING ==================== #
                is_composite = material isa Vector{Material}
                if is_composite
                    # Validate composite material input
                    isnothing(connect_elem_phase) && error("connect_elem_phase required for composites")
                    length(connect_elem_phase) == size(connectivity, 1) || error("connect_elem_phase length mismatch")
                    
                    unique_phases = sort(unique(connect_elem_phase))
            
                    length(material) == length(unique_phases) || error("Material count mismatch with phases")
                    phase_to_idx = Dict(zip(unique_phases, 1:length(unique_phases)))
                    
                    # Validate material consistency
                    ref_dofs = get_dofs_per_node(material[1])
                    ref_B_types = material[1].B_types
                    for mat in material
                        get_dofs_per_node(mat) == ref_dofs || error("Mismatch in dofs_per_node")
                        mat.B_types == ref_B_types || error("Mismatch in B_types")
                    end
                    dofs_per_node = ref_dofs
                else
                    dofs_per_node = get_dofs_per_node(material)
                end

            # ==================== PROBLEM SETUP ==================== #
                total_dofs = size(nodes, 1) * dofs_per_node
                compute_F = !isnothing(NodalForces) || !isnothing(BoundaryFace) || !isnothing(PointForces) ||
                            any(should_apply_internal_force(m) for m in (material isa Vector ? material : [material]))
                compute_BoundaryFace = !isnothing(BoundaryFace)

                n_nodes_per_element = size(connectivity, 2)
                n_dofs_per_element = n_nodes_per_element * dofs_per_node
                nb_elm = size(connectivity, 1)
                entries_per_elem = n_dofs_per_element^2
                total_entries = nb_elm * entries_per_elem

            # ==================== BOUNDARY FACE HANDLING ==================== #
                if compute_BoundaryFace
                    jacobian_fcache = Dict{Symbol, Vector}()
                    global_dofsforces = Dict{Symbol, Matrix{Int}}()
                    
                    for (name, bf) in BoundaryFace
                    
                        gauss_fdata = shape_data(bf.element_type, bf.int_order, bf.dim-1)
                        jacobian_fcache[name] = jacobian_data(bf.element_border, bf.nodes, gauss_fdata, dim; elem_dim = bf.dim)
                        global_dofsforces[name] = assemble_global_dofs(bf.element_border, material isa Vector ? material[1] : material)
                    end
                else
                    
                    jacobian_fcache = nothing
                    global_dofsforces = nothing
                    gauss_fdata = nothing
                end

            # ==================== PRECOMPUTATION ==================== #
        
                gauss_data = Geometric_Data.gauss_data
                jacobian_cache = Geometric_Data.jacobian_cache
                B_dicts = Geometric_Data.B_dicts
                ref_mat = is_composite ? material[1] : material
                global_dofs = assemble_global_dofs(connectivity, ref_mat)

                # ==================== GLOBAL MATRIX  ==================== #
                    I = Vector{Int}(undef, total_entries)
                    J = Vector{Int}(undef, total_entries)
                    V = Vector{Float64}(undef, total_entries)

                # ==================== FORCE MATRIX ==================== #
                    if compute_F
                        if length(ref_mat.type) != 1 && ref_mat.symmetry == :out_of_plane
                            # For coupled out-of-plane systems, create separate force storage per field
                            local_F = Dict{Symbol, Vector{Vector{Float64}}}()
                            for bt in ref_mat.B_types
                                local_F[bt] = [zeros(total_dofs) for _ in 1:Threads.nthreads()]
                            end
                        else
                            # For all other cases, use single force storage
                            local_F = Dict{Symbol, Vector{Vector{Float64}}}()
                            local_F[:default] = [zeros(total_dofs) for _ in 1:Threads.nthreads()]
                        end
                    else
                        local_F = nothing
                    end
                
            # ==================== ELEMENT PROCESSING ==================== #
                Threads.@threads for e in 1:nb_elm
                    tid = Threads.threadid()
                    
                    # Material selection
                    mat = is_composite ? material[phase_to_idx[connect_elem_phase[e]]] : material
                    
                    # Element data
                    B_dict = B_dicts[e]
                    J_data = jacobian_cache[e]
                    elem_dof = global_dofs[e, :]
                    
                    # Stiffness computation
                    Kₑ = compute_element_stiffness(mat, B_dict, J_data, gauss_data)
                    
                    # Mass computation 
                    Mₑ = isnothing(mat.mass_properties) ? nothing : compute_element_mass(mat, J_data, gauss_data)
   
                    # Global matrix assembly
                    idx_start = (e-1) * entries_per_elem
                
                    for j in 1:n_dofs_per_element, i in 1:n_dofs_per_element
                        idx = idx_start + (j-1)*n_dofs_per_element + i
                        I[idx] = elem_dof[i]
                        J[idx] = elem_dof[j]
                        V[idx] = Kₑ[i, j]
                    end
                    
                    # Force computation
                    if compute_F
                        
                        if !isnothing(NodalForces)
                            for (face_name, fsrc) in NodalForces
                                
                                isnothing(fsrc.fᵥ) && isnothing(fsrc.fₜ) && continue
                                
                                face_conn = Int[]
                                jac_fcache = nothing
                                face_coords = nothing
                                
                                # Boundary face matching
                                if !isnothing(fsrc.fₜ) && compute_BoundaryFace && haskey(BoundaryFace, face_name)
                                    
                                    bf = BoundaryFace[face_name]
                                    bf_faces = global_dofsforces[face_name]
                                    bf_jacs = jacobian_fcache[face_name]
                                    
                                    for (i, face) in enumerate(eachrow(bf_faces))
                                        if all(dof -> dof in elem_dof, face)
                                            face_conn = collect(face)
                                            jac_fcache = bf_jacs[i]
                                            node_ids = unique(floor.(Int, (face_conn .- 1) ./ dofs_per_node) .+ 1)
                                            face_coords = nodes[node_ids, :]
                                            break
                                        end
                                    end
                                end
                             
                                isempty(face_conn) && continue
                                
                                fᵉ = compute_element_force(
                                    Val(resolve_material_type(mat)), mat,
                                    J_data, gauss_data,
                                    jac_fcache, gauss_fdata,
                                    fsrc, face_conn, elem_dof, B_dict,
                                    face_coords)
                                
                                local_F[:default][tid][face_conn] .+= fᵉ
                            end
                        else
                            
                            # When NodalForces is not provided, use default forces for special cases.
                            # For instance, if material is out_of_plane or requires an internal force,
                            # define a default force vector (here zeros of appropriate length).
                            jac_fcache = nothing
                            gauss_fdata_local = nothing
                            face_conn = nothing
                            face_coords = nothing
                            default_f = should_apply_internal_force(mat) ? (fᵥ = zeros(mat.dim), fₜ = nothing) : (fᵥ = nothing, fₜ = nothing)
                        
                            fᵉ = compute_element_force(
                                    Val(resolve_material_type(mat)), mat,
                                    J_data, gauss_data,
                                    jac_fcache, gauss_fdata_local,
                                    default_f, face_conn, elem_dof, B_dict,
                                    face_coords)
                                
                            # local_F[tid][elem_dof] .+= fᵉ
                            
                            if length(ref_mat.type) != 1 && ref_mat.symmetry == :out_of_plane
                                # For coupled out-of-plane systems: fᵉ is a tuple of field forces
                                for (i, bt) in enumerate(ref_mat.B_types)
                                    
                                    local_F[bt][tid][elem_dof] .+= fᵉ[i]
                                end
                            else
                                # For all other cases: fᵉ is a single vector
                                local_F[:default][tid][elem_dof] .+= fᵉ
                            end
                        end
                    end
                end
            

            # ==================== GLOBAL ASSEMBLY ==================== #
                global_K = sparse(I, J, V, total_dofs, total_dofs)
                if compute_F
                    if length(ref_mat.type) != 1 && ref_mat.symmetry == :out_of_plane
                        # For coupled out-of-plane systems: combine and return separate forces
                        global_F_fields = []
                        for bt in ref_mat.B_types
                            push!(global_F_fields, reduce(+, local_F[bt]))
                        end
                        return global_K, global_F_fields...
                         # Add point forces for out-of-plane systems (not supported yet)
                        if !isnothing(PointForces)
                            error("PointForces not yet supported for out-of-plane coupled systems")
                        end
                    return global_K, global_F_fields...
                    else
                        # For all other cases: return single force vector
                        global_F = reduce(+, local_F[:default])
                         # ===== ADD CONCENTRATED FORCES HERE ===== #
                    if !isnothing(PointForces)
                        F_point = zeros(total_dofs)
                        for (node, f_vec) in PointForces
                            # Calculate base DOF index for this node
                            base_dof = (node - 1) * dofs_per_node
                            
                            if length(f_vec) != dofs_per_node
                                error("Force vector length mismatch for node $node: Expected $dofs_per_node DOFs")
                            end
                            
                            # Add force vector to global force
                            F_point[base_dof+1:base_dof+dofs_per_node] .+= f_vec
                        end
                        global_F .+= F_point
                    end
                    # ===== END CONCENTRATED FORCES ===== #
                        return global_K, global_F
                    end
                else
                    return global_K
                end
        end

# =============================================
#  Element force computation
# =============================================
    # ============================================================
    # Helper Functions for Force Computation
    # ============================================================
        # ==================== SHAPES N UTILITY FUNCTION  ==================== #
            function build_N_matrix(N_face::Vector{Float64}, dim::Int)
                n_nodes = length(N_face)
                N = zeros(dim, dim * n_nodes)
                for i in 1:n_nodes
                    for d in 1:dim
                        N[d, dim*(i-1)+d] = N_face[i]
                    end
                end
                return N
            end
        # ==================== TRACTION TYPE DEFINITIONS ==================== #
            abstract type AbstractTraction end

            struct ConstantTraction <: AbstractTraction
                value::Vector{Float64}
            end

            struct FunctionTraction{F} <: AbstractTraction
                f::F
            end

        # ==================== TRACTION EVALUATION METHODS ==================== #
            function evaluate(t::Nothing, x::AbstractVector)
                return nothing
            end

            function evaluate(t::ConstantTraction, x::AbstractVector)
                return t.value
            end

            function evaluate(t::FunctionTraction, x::AbstractVector)
                try
                    return t.f(x...)
                catch
                    error("Traction function must accept as many args as spatial dimensions (e.g. (x, y))")
                end
            end

        # ==================== FORCE COMPUTATION UTILITIES ==================== #
            function force_computation(input)
                if input === nothing
                    return nothing
                elseif input isa AbstractVector
                    return ConstantTraction(input)
                elseif input isa Function
                    return FunctionTraction(input)
                else
                    error("Unsupported traction definition: $input")
                end
            end

            function compute_x_qp(N_face_matrix::Matrix{Float64}, face_coords::Matrix{Float64}; tol=1e-10)
                n_nodes, dim = size(face_coords)
                N_scalar = zeros(n_nodes)
                
                # Extract scalar shape functions
                for i in 1:n_nodes
                    N_scalar[i] = N_face_matrix[1, (i - 1)*dim + 1]  # First component only
                end

                # Compute coordinates at quadrature point
                x_qp = zeros(dim)
                for d in 1:dim
                    coords_d = face_coords[:, d]
                    if maximum(coords_d) - minimum(coords_d) < tol
                        x_qp[d] = coords_d[1]  # Coordinate is constant
                    else
                        for i in 1:n_nodes
                            x_qp[d] += N_scalar[i] * coords_d[i]
                        end
                    end
                end

                return x_qp
            end

        # ==================== VOLUME FORCE COMPUTATIONS ==================== #
            function compute_volume_force_vector(N, scaling, fᵥ, dim, n_nodes)
                f = zeros(eltype(N), n_nodes )
                Nn = size(N,1)
                @inbounds for i in 1:Nn
                    Ni = N[i] * scaling
                    base = (i-1)*dim
                    @simd for d in 1:dim
                        f[base + d] = Ni * fᵥ[d]
                    end
                end
                return f
            end

            function compute_volume_force_scalar(N, scaling, fᵥ)
               f = zeros(eltype(N), length(N))
                @inbounds @simd for i in eachindex(N)
                    f[i] = N[i] * fᵥ * scaling
                end
                return f
            end

        # ==================== SPECIALIZED FORCE COMPUTATIONS ==================== #

            function compute_out_of_plane_force(B, tensor, scaling)
               
                M = zeros(size(tensor, 2))
                M[end] = 1  # Out-of-plane component
                return -B' * tensor * M * scaling
            end

            function thermal_expansion_force!(B, material, scaling)

                α = material.properties[:α]
                C = material.tensors[1]
                ε_th = ones(size(C, 1))
                ε_th[end] = 0


                    # Accumulate force contribution
                   return B' * C * ε_th * scaling
            end

        # ==================== SURFACE FORCE COMPUTATIONS ==================== #
            function compute_surface_force_vector(
                N_face::Union{Vector{Float64}, Matrix{Float64}},
                scaling::Float64,
                fₜ::AbstractTraction,
                dim::Int,
                face_coords::Matrix{Float64},
                n)
              
                f = zeros(size(N_face, 1))

                x_qp = compute_x_qp(N_face, face_coords)
                fₜ_val = fₜ isa AbstractTraction ? evaluate(fₜ, x_qp) : zeros(dim)
                fₜ_val = - fₜ_val * n
                
                f = N_face' * fₜ_val * scaling
              
                return f
            end

            function compute_surface_force_scalar(N_face, scaling, fₜ, face_conn, elem_conn)
                f = zeros(size(N_face, 1))
                x_qp = compute_x_qp(N_face, face_coords)
                fₜ_val = fₜ isa AbstractTraction ? evaluate(fₜ, x_qp) : zeros(dim)
                f = N_face * fₜ_val * scaling
                return f
            end
            # Helper function for permutation vector
            function get_node_wise_permutation(material::Material, n_nodes::Int)
                fields = material.B_types
                dofs_per_field = [get_dofs_per_node(Material(material.type, material.dim, material.symmetry, material.properties, Vector{AbstractMatrix{Real}}(), [bt],nothing)) for bt in fields]
                dofs_per_node = sum(dofs_per_field)
                total_dofs = n_nodes * dofs_per_node

                # Build permutation vector: field-wise -> node-wise
                p = zeros(Int, total_dofs)
                offset_field = 0
                for (j, bt) in enumerate(fields)
                    dof_field = dofs_per_field[j]
                    for node in 0:(n_nodes-1)
                        offset_in_node = sum(dofs_per_field[1:j-1])
                        for k in 1:dof_field
                            field_wise_idx = offset_field + node * dof_field + k
                            node_wise_idx = node * dofs_per_node + offset_in_node + k
                            p[node_wise_idx] = field_wise_idx
                        end
                    end
                    offset_field += n_nodes * dof_field
                end
                return p
            end

    # ============================================================
    # Force Computation
    # ============================================================
        # ==================== Elastic ==================== #
            function compute_element_force(::Val{:elastic}, material::Material,
                jac_data, gauss_data,
                jac_fdata, gauss_fdata,
                forces,
                face_conn::Union{Vector{Int}, Nothing},
                elem_conn::Union{Vector{Int}, Nothing},
                B_dict,
                face_coords::Union{Matrix{Float64}, Nothing})

                # Extract element properties
                B = B_dict[:strain]
                n_dofs = size(B[1], 2)
                dim = material.dim            
                f = zeros(n_dofs)
                # Volume force
                if !isnothing(forces.fᵥ)
                   
                    for qp in 1:length(gauss_data.weights)
                        detJ, _ = jac_data[qp]
                        scaling =  abs(detJ) * gauss_data.weights[qp]
                        
                        if material.symmetry == :out_of_plane
                            
                            f += compute_out_of_plane_force(B[qp], material.tensors[1], scaling)                        

                        elseif haskey(material.properties, :α  )
                            f += thermal_expansion_force!(B[qp], material, scaling)
                        else
                        # newman
                            N = gauss_data.N[qp]
                            f += compute_volume_force_vector(N, scaling, forces.fᵥ, dim, n_dofs)
                        end
                    end
                end
                
                # Traction force
                if !isnothing(forces.fₜ)    
                    f =  zeros(n_dofs)
                    center = mean(face_coords, dims=1)
                    for qp in 1:length(gauss_fdata.weights)
                        N_face = gauss_fdata.N[qp]                           
                        detJ, _, J = jac_fdata[qp]    
                       
                        # === Compute the normal vector  ===
                        if dim == 2  # 2D case (1D edge)
                            t = J[1, :]
                            n = [t[2], -t[1]]
                            gauss_pos = N_face' * face_coords
                            outward_check = dot(n, gauss_pos' - center[:]) < 0
                            if !outward_check
                                n = -n 
                            end
                        elseif dim == 3  # 3D case (2D face)                           
                            t₁ = J[:, 1]
                            t₂ = J[:, 2]
                            n = cross(t₁, t₂)
                        end
                        
                        
                        # Normalize the normal vector (unit normal)
                        n = n / norm(n)
            
                        scaling = abs(detJ) * gauss_fdata.weights[qp]
                        N_face = build_N_matrix(N_face, dim)   

                        f += compute_surface_force_vector(N_face, scaling, forces.fₜ, dim, face_coords, n)
                    end
                
                end
                  
                return f
            end
        # ==================== THERMAL ==================== #
            function compute_element_force(::Val{:thermal}, material::Material,
                jac_data, gauss_data,
                jac_fdata, gauss_fdata,
                forces,
                face_conn::Union{Vector{Int}, Nothing},
                elem_conn::Vector{Int},
                B_dict,
                face_coords::Matrix{Float64})

                # Extract element properties
                B = B_dict[:temperature_gradient]
                n_dofs = size(B[1], 2)
                dim = material.dim

                # Initialize force vector
                f = zeros(n_dofs)

                # Volume heat source
                if !isnothing(forces.fᵥ)
                    for qp in 1:length(gauss_data.weights)
                        detJ, _ = jac_data[qp]
                        scaling = abs(detJ) * gauss_data.weights[qp]
                        
                        if material.symmetry == :out_of_plane
                            f += compute_out_of_plane_force(B[qp], material.tensors[1], scaling)
                        else
                            N = gauss_data.N[qp]
                            f += compute_volume_force_scalar(N, scaling, forces.fᵥ)
                        end
                    end
                end
                
                # Surface heat flux
                if !isnothing(forces.fₜ)
                    f_face = zeros(length(face_conn))
                    for qp in 1:length(gauss_fdata.weights)
                        N_face = gauss_fdata.N[qp]
                        detJ, _ = jac_fdata[qp]
                        scaling = abs(detJ) * gauss_fdata.weights[qp]
                        f_face += compute_surface_force_scalar(N_face, scaling, forces.fₜ, face_conn, elem_conn)
                    end
                    f = f_face
                end
                
                return f
            end

        # ==================== PIEZOELECTRIC ==================== #
            function compute_element_force(::Val{:piezoelectric}, material::Material,
                jac_data, gauss_data,
                jac_fdata, gauss_fdata,
                forces,
                face_conn::Union{Vector{Int}, Nothing},
                elem_conn::Union{Vector{Int}, Nothing},
                B_dict,
                face_coords::Union{Matrix{Float64}, Nothing})

                # Get element properties
                B_u = B_dict[:strain]
                B_ϕ = B_dict[:electric_field]
                n_dofs_u = size(B_u[1], 2)
                n_dofs_ϕ = size(B_ϕ[1], 2)
                dim = material.dim
                n_nodes = n_dofs_u ÷ dim

                # Volume force - mechanical
                if !isnothing(forces.fᵥ)
                      # Initialize separate force vectors
                    f_u = zeros(n_dofs_u + n_dofs_ϕ)  # Mechanical forces
                    f_ϕ = zeros(n_dofs_ϕ + n_dofs_u)  # Electrical forces
                    for qp in 1:length(gauss_data.weights)
                        detJ, _ = jac_data[qp]
                        scaling = abs(detJ) * gauss_data.weights[qp]
                        
                        if material.symmetry == :out_of_plane
                            C, ϵ, e = material.tensors
                            B_u_qp = B_u[qp]
                            B_e_qp = B_ϕ[qp] 
                           
                            F_u = compute_out_of_plane_force(B_u_qp, C, scaling)
                            G_u = compute_out_of_plane_force(B_u_qp, transpose(e), scaling)
                            F_ϕ = compute_out_of_plane_force(B_e_qp, e, scaling)
                            G_ϕ = compute_out_of_plane_force(B_e_qp, ϵ, scaling)

                            f_u += vcat(F_u, F_ϕ)
                            f_ϕ += vcat(-G_u, G_ϕ)
                            perm = get_node_wise_permutation(material, n_nodes)
                            f_u_perm = f_u[perm]
                            f_ϕ_perm = f_ϕ[perm]
                            fe = (f_u_perm, f_ϕ_perm)
                
                        elseif haskey(material.properties, :thermal_expansion)
                            # Place_Holder need to be implement latter    
                        else
                          
                            N = gauss_data.N[qp]
                            f_u += compute_volume_force_vector(N, scaling, forces.fᵥ, dim, n_dofs)
                            f_ϕ += compute_volume_force_vector(N, scaling, forces.fᵥ, dim, n_dofs)
                        end
                    end
                end
                        
                # Surface forces 
                if !isnothing(face_conn) && !isnothing(face_coords)
                      # Initialize separate force vectors
                    f_u = zeros(n_dofs_u)  # Mechanical forces
                    f_ϕ = zeros(n_dofs_ϕ)  # Electrical forces
                    # Mechanical traction
                    if !isnothing(forces.fₜ)
                        center = mean(face_coords, dims=1)
                        for qp in 1:length(gauss_fdata.weights)
                            N_face = gauss_fdata.N[qp] 
                            detJ, _, J = jac_fdata[qp]    
                       
                            # === Compute the normal vector  ===
                            if dim == 2  # 2D case (1D edge)
                                t = J[1, :]
                                n = [t[2], -t[1]]
                                gauss_pos = N_face' * face_coords
                                outward_check = dot(n, gauss_pos' - center[:]) < 0
                                if !outward_check
                                    n = -n 
                                end
                            elseif dim == 3  # 3D case (2D face)                           
                                t₁ = J[:, 1]
                                t₂ = J[:, 2]
                                n = cross(t₁, t₂)
                            end
                                          
                            # Normalize the normal vector (unit normal)
                            n = n / norm(n)                  
                            N_face = build_N_matrix(N_face, dim)
                            scaling =   abs(detJ) * gauss_fdata.weights[qp]           
                            f_u += compute_surface_force_vector(N_face, scaling, forces.fₜ, dim, face_coords,n)
                        end
                    end
                      
                    # Electrical surface charge
                    if !isnothing(forces.f_q)
                        center = mean(face_coords, dims=1)
                        for qp in 1:length(gauss_fdata.weights)
                        N_face = gauss_fdata.N[qp]  
                        detJ, _, J = jac_fdata[qp]    
                       
                        # === Compute the normal vector  ===
                        if dim == 2  # 2D case (1D edge)
                            t = J[1, :]
                            n = [t[2], -t[1]]
                            gauss_pos = N_face' * face_coords
                            outward_check = dot(n, gauss_pos' - center[:]) < 0
                            if !outward_check
                                n = -n 
                            end
                        elseif dim == 3  # 3D case (2D face)                           
                            t₁ = J[:, 1]
                            t₂ = J[:, 2]
                            n = cross(t₁, t₂)
                        end
                                      
                        # Normalize the normal vector (unit normal)
                        n = n / norm(n)
                        scaling =   abs(detJ) * gauss_fdata.weights[qp]
                        f_ϕ += compute_surface_force_vector(N_face, scaling, forces.f_q, dim, face_coords,n)
                        
                    end
                    end

                    fe = hcat(f_u, f_ϕ)
                    perm = get_node_wise_permutation(material, n_nodes)
                    fe = fe[perm]
                end

              
                return fe
            end
        # ==================== POROELASTIC ==================== #
            function compute_element_force(::Val{:poroelastic}, material::Material,
                jac_data, gauss_data,
                jac_fdata, gauss_fdata,
                forces,
                face_conn::Union{Vector{Int}, Nothing},
                elem_conn::Vector{Int},
                B_dict,
                face_coords::Matrix{Float64})

                # Extract element properties
                B = B_dict[:strain]
                n_dofs = size(B[1], 2)
                dim = material.dim            

                # Volume force
                if !isnothing(forces.fᵥ)
                    f = zeros(n_dofs)
                    for qp in 1:length(gauss_data.weights)
                        detJ, _ = jac_data[qp]
                        scaling =  abs(detJ) * gauss_data.weights[qp]
                        
                        if material.symmetry == :out_of_plane
                            
                            f += compute_out_of_plane_force(B[qp], material.tensors[1], scaling)                        

                        elseif haskey(material.properties, :α  )
                            f += thermal_expansion_force!(B[qp], material, scaling)
                        else
                        # newman
                            N = gauss_data.N[qp]
                            f += compute_volume_force_vector(N, scaling, forces.fᵥ, dim, n_dofs)
                        end
                    end
                end
                
                # Traction force
                if !isnothing(forces.fₜ)
                    
                    f =  zeros(length(face_conn))
                    center = mean(face_coords, dims=1)
                    for qp in 1:length(gauss_fdata.weights)
                        N_face = gauss_fdata.N[qp]                           
                        detJ, _, J = jac_fdata[qp]    
                       
                        # === Compute the normal vector  ===
                        if dim == 2  # 2D case (1D edge)
                            t = J[1, :]
                            n = [t[2], -t[1]]
                            gauss_pos = N_face' * face_coords
                            outward_check = dot(n, gauss_pos' - center[:]) < 0
                            if !outward_check
                                n = -n 
                            end
                        elseif dim == 3  # 3D case (2D face)                           
                            t₁ = J[:, 1]
                            t₂ = J[:, 2]
                            n = cross(t₁, t₂)
                        end
                        
                        
                        # Normalize the normal vector (unit normal)
                        n = n / norm(n)
            
                        scaling = abs(detJ) * gauss_fdata.weights[qp]
                        N_face = build_N_matrix(N_face, dim)   

                        f += compute_surface_force_vector(N_face, scaling, forces.fₜ, dim, face_coords, n)
                    end
                end
                
                return f
            end
        # ==================== VISCOELASTIC ==================== #
            function compute_element_force(::Val{:viscoelastic}, material::Material,
                    jac_data, gauss_data,
                    jac_fdata, gauss_fdata,
                    forces,
                    face_conn::Union{Vector{Int}, Nothing},
                    elem_conn::Union{Vector{Int}, Nothing},
                    B_dict,
                    face_coords::Union{Matrix{Float64}, Nothing})

                 # Place_Holder need to be implement latter    
            end