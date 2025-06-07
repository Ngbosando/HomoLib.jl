include("Precompute_data.jl")
include("material_model.jl")
include("resolve_material_type.jl")
# =============================================
# Element Stiffness Calculation
# =============================================

# Dispatch based on material type
function compute_element_stiffness(material::Material, B_dict::Dict, jacobian_data, gauss_data)
    T = resolve_material_type(material.type)
    return compute_element_stiffness(Val(T), material, B_dict, jacobian_data, gauss_data)
end
function compute_element_stiffness(::Val{:elastic}, material::Material, B_dict, jacobian_data, gauss_data)
    C = material.tensors[1]
    B = B_dict[:strain]
    n_gauss = length(gauss_data.weights)
    n_dofs = size(B[1], 2)
    K = zeros(n_dofs, n_dofs)
    

    for qp in 1:n_gauss
        detJ, _ = jacobian_data[qp]
        scaling = 0.5 * abs(detJ) * gauss_data.weights[qp]
        K .+= B[qp]' * C * B[qp] * scaling
    end
    return K
end
function compute_element_stiffness(::Val{:piezoelectric}, material::Material, B_dict, jacobian_data, gauss_data)
    C, e, ϵ = material.tensors
    B_u = B_dict[:strain]
    B_e = B_dict[:electric_field]
    
    n_dofs_u = size(B_u[1], 2)
    n_dofs_e = size(B_e[1], 2)

    K_uu = zeros(n_dofs_u, n_dofs_u)
    K_ue = zeros(n_dofs_u, n_dofs_e)
    K_ee = zeros(n_dofs_e, n_dofs_e)
    # K = zeros(n_dofs_u + n_dofs_e, n_dofs_u + n_dofs_e)
    for qp in 1:length(gauss_data.weights)
        detJ, _ = jacobian_data[qp]
        scaling = 0.5 * abs(detJ) * gauss_data.weights[qp]
    
        # Mechanical block
        K_uu .+= B_u[qp]' * C * B_u[qp] * scaling
       
        # Coupling blocks
        K_ue .+= B_u[qp]' * ϵ * B_e[qp] * scaling
        # K[n_dofs_u+1:end, 1:n_dofs_u] .+= B_e[qp]' * e' * B_u[qp] * scaling
      
        # Electrical block
        K_ee .+= B_e[qp]' * e * B_e[qp] * scaling

    end
  
    K = [K_uu K_ue; K_ue' K_ee]

    return K
end
function compute_element_stiffness(::Val{:poroelastic}, material::Material, B_dict, jacobian_data, gauss_data)
    C, α, inv_M = material.tensors
    B_u = B_dict[:strain]
    B_p = B_dict[:pressure]
    
    n_dofs_u = size(B_u[1], 2)
    n_dofs_p = size(B_p[1], 2)
    K = zeros(n_dofs_u + n_dofs_p, n_dofs_u + n_dofs_p)
    
    for qp in 1:length(gauss_data.weights)
        detJ, _ = jacobian_data[qp]
        scaling = 0.5 * abs(detJ) * gauss_data.weights[qp]
        
        # Mechanical block
        K[1:n_dofs_u, 1:n_dofs_u] .+= B_u[qp]' * C * B_u[qp] * scaling
        
        # Coupling blocks
        K[1:n_dofs_u, n_dofs_p+1:end] .+= B_u[qp]' * α * B_p[qp] * scaling
        K[n_dofs_p+1:end, 1:n_dofs_u] .+= B_p[qp]' * α' * B_u[qp] * scaling
        
        # Pressure block
        K[n_dofs_p+1:end, n_dofs_p+1:end] .+= B_p[qp]' * inv_M * B_p[qp] * scaling
    end
    return K
end
function compute_element_stiffness(::Val{:thermal}, material::Material, B_dict, jacobian_data, gauss_data)
    κ = material.tensors[1]  # Thermal conductivity matrix
    B = B_dict[:temperature_gradient]
    n_gauss = length(gauss_data.weights)
    n_nodes = size(B[1], 2)
    K = zeros(n_nodes, n_nodes)
    

    for qp in 1:n_gauss
        detJ, _ = jacobian_data[qp]
        
        scaling = 0.5 * abs(detJ) * gauss_data.weights[qp]
        K .+= B[qp]' * κ * B[qp] * scaling
    end
    return K
end
function compute_element_stiffness(::Val{:viscoelastic}, material::Material, B_dict, jacobian_data, gauss_data)
    C∞ = material.tensors[1]  # Long-term stiffness tensor
    B = B_dict[:strain]
    n_gauss = length(gauss_data.weights)
    n_dofs = size(B[1], 2)
    K = zeros(n_dofs, n_dofs)

    for qp in 1:n_gauss
        detJ, _ = jacobian_data[qp]
        scaling = 0.5 * abs(detJ) * gauss_data.weights[qp]
        K .+= B[qp]' * C∞ * B[qp] * scaling
    end
    return K
end

function compute_element_stiffness(::Val{:thermoelastic}, material::Material, B_dict, jacobian_data, gauss_data)
    C, β, κ = material.tensors
    B_u = B_dict[:strain]
    B_t = B_dict[:temperature_gradient]

    n_dofs_u = size(B_u[1], 2)
    n_dofs_t = size(B_t[1], 2)
    K_uu = zeros(n_dofs_u , n_dofs_u )
    K_ut = zeros(n_dofs_u + n_dofs_t, n_dofs_u + n_dofs_t)
    K_tt = zeros(n_dofs_u + n_dofs_t, n_dofs_u + n_dofs_t)

    for qp in 1:length(gauss_data.weights)
        detJ, _ = jacobian_data[qp]
        scaling = 0.5 * abs(detJ) * gauss_data.weights[qp]

        # Mechanical block
        K[1:n_dofs_u, 1:n_dofs_u] .+= B_u[qp]' * C * B_u[qp] * scaling

        # Coupling blocks (β is strain_dim × 1, expand if needed)
        β_mat = size(β, 2) == 1 ? β * ones(1, size(B_t[qp], 1)) : β
        K[1:n_dofs_u, n_dofs_u+1:end] .+= B_u[qp]' * β_mat * B_t[qp] * scaling
        K[n_dofs_u+1:end, 1:n_dofs_u] .+= B_t[qp]' * β_mat' * B_u[qp] * scaling

        # Thermal block
        K[n_dofs_u+1:end, n_dofs_u+1:end] .+= B_t[qp]' * κ * B_t[qp] * scaling
    end

    return K
end
function compute_element_stiffness(::Val{:thermo_piezoelectric}, material::Material, B_dict, jacobian_data, gauss_data)
    C, e, ϵ, β, κ = material.tensors
    B_u = B_dict[:strain]
    B_e = B_dict[:electric_field]
    B_t = B_dict[:temperature_gradient]

    nu, ne, nt = size(B_u[1],2), size(B_e[1],2), size(B_t[1],2)
    K = zeros(nu + ne + nt, nu + ne + nt)

    for qp in 1:length(gauss_data.weights)
        detJ, _ = jacobian_data[qp]
        w = 0.5 * abs(detJ) * gauss_data.weights[qp]

        u = 1:nu
        e = nu+1 : nu+ne
        t = nu+ne+1 : nu+ne+nt

        # Mechanical
        K[u, u] .+= B_u[qp]' * C * B_u[qp] * w

        # Electric
        K[e, e] .+= B_e[qp]' * ϵ * B_e[qp] * w

        # Thermal
        K[t, t] .+= B_t[qp]' * κ * B_t[qp] * w

        # Couplings
        K[u, e] .+= B_u[qp]' * e * B_e[qp] * w
        K[e, u] .+= B_e[qp]' * e' * B_u[qp] * w

        K[u, t] .+= B_u[qp]' * β * B_t[qp] * w
        K[t, u] .+= B_t[qp]' * β' * B_u[qp] * w
    end
    return K
end

# =============================================
# DOF Handling
# =============================================

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
# implies internal forces even without user input.
# =============================================

function should_apply_internal_force(material::Material)::Bool
    material.symmetry == :out_of_plane || haskey(material.properties, :α)
end

# =============================================
# Global dof assembling 
# =============================================


function assemble_global_dofs(connectivity::Matrix{Int}, material::Material)
    n_elements, nodes_per_element = size(connectivity)
    dofs_per_node = get_dofs_per_node(material)
    global_dofs = Matrix{Int}(undef, n_elements, nodes_per_element * dofs_per_node)
    
    field_offsets = cumsum([get_dofs_per_node(Material(material.type, material.dim, material.symmetry, material.properties, [], [bt])) 
                      for bt in material.B_types])
    
    for e in 1:n_elements
        elem_dofs = Int[]
        for node in connectivity[e, :]
            base_dof = (node - 1) * dofs_per_node
            for (i, bt) in enumerate(material.B_types)
                n_dofs = get_dofs_per_node(Material(material.type, material.dim, material.symmetry, material.properties, [], [bt]))
                append!(elem_dofs, base_dof .+ (1:n_dofs) .+ (i > 1 ? field_offsets[i-1] : 0))
            end
        end
    
        global_dofs[e, :] = elem_dofs
    end
    
    return global_dofs
end

# =============================================
#  Global Matrix Assembler
# =============================================

function assemble_global_matrix(
    connectivity::Matrix{Int},
    nodes::Matrix{Float64},
    element_type::Symbol,
    integration_order::Int,
    material::Union{Material, Vector{Material}},
    dim::Int,
    connect_elem_phase::Union{Nothing, Vector{Int}};
    BoundaryFace::Union{NamedTuple, Nothing} = nothing,
    NodalForces::Union{NamedTuple, Nothing} = nothing)

   # ========== Phase Handling ========== #
    is_composite = material isa Vector{Material}
    if is_composite
        connect_elem_phase === nothing && error("connect_elem_phase required for composite materials")
        length(connect_elem_phase) == size(connectivity, 1) || error("connect_elem_phase length mismatch")

        unique_phases = sort(unique(connect_elem_phase))
        length(material) == length(unique_phases) || error("Material count mismatch with phases")

        phase_to_idx = Dict(zip(unique_phases, 1:length(unique_phases)))

        # Validate all materials are consistent
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

    total_dofs = size(nodes, 1) * dofs_per_node
    compute_F = (
    !isnothing(NodalForces) || !isnothing(BoundaryFace) ||
    any(should_apply_internal_force(m) for m in (material isa Vector ? material : [material])))


    gauss_data = shape_data(element_type, integration_order, dim)
    jacobian_cache = jacobian_data(connectivity, nodes, gauss_data)
    ref_mat = is_composite ? material[1] : material
    B_dicts = build_B_matrices(nodes, connectivity, ref_mat, gauss_data, jacobian_cache)
    global_dofs = assemble_global_dofs(connectivity, ref_mat)

    local_K = [spzeros(total_dofs, total_dofs) for _ in 1:Threads.nthreads()]
    local_F = compute_F ? [zeros(total_dofs) for _ in 1:Threads.nthreads()] : nothing
    
    Threads.@threads for e in 1:size(connectivity, 1)
        tid = Threads.threadid()
        mat = is_composite ? material[phase_to_idx[connect_elem_phase[e]]] : material
    
        
        B_dict = B_dicts[e]
        J_data = jacobian_cache[e]
        elem_dof = global_dofs[e, :]

        Ke = compute_element_stiffness(mat, B_dict, J_data, gauss_data)
     
        for i in eachindex(elem_dof), j in eachindex(elem_dof)
            local_K[tid][elem_dof[i], elem_dof[j]] += Ke[i, j]
        end

        if compute_F
            # Dummy face if needed
            face_conn = BoundaryFace === nothing ? elem_dof[1:2] : BoundaryFace.lin_elements[e, :]
            mat_type = resolve_material_type(mat.type)
            NodalForces = should_apply_internal_force(mat) ? NodalForces = (fᵥ = [0.0 0.0],fₜ = nothing) : NodalForces
            fᵉ = compute_element_force(Val(mat_type), mat, J_data, gauss_data,
                           nothing, nothing,
                           NodalForces,
                           face_conn, connectivity[e, :], B_dict)
           
            for i in eachindex(elem_dof)
                local_F[tid][elem_dof[i]] += fᵉ[i]
            end
        end
    end

    global_K = reduce(+, local_K)
    global_F = compute_F ? reduce(+, local_F) : nothing
    return compute_F ? (global_K, global_F) : global_K
end

# =============================================
# Functions for force computation
# =============================================

# ============================================================
# Helper Functions for Force Computation
# ============================================================

function compute_volume_force_vector(N, scaling, fᵥ, dim, n_nodes)
    f = zeros(n_nodes * dim)
    for i in 1:n_nodes
        for d in 1:dim
            idx = (i-1)*dim + d
            f[idx] = N[i] * fᵥ[d] * scaling
        end
    end
    return f
end

function compute_volume_force_scalar(N, scaling, fᵥ)
    Nf = length(N)
    f = zeros(Nf)
    for i in 1:Nf
        f[i] = N[i] * fᵥ * scaling
    end
    return f
end

function compute_out_of_plane_force(B, tensor, scaling)
    M = zeros(size(tensor, 1))
    M[end] = 1  # Out-of-plane component
    return -B' * tensor * M * scaling
end

function compute_surface_force_vector(N_face, scaling, fₜ, dim, face_conn, elem_conn)
    n_face_nodes = length(face_conn)
    f = zeros(length(elem_conn) * dim)
    for i_face in 1:n_face_nodes
        i_elem = findfirst(isequal(face_conn[i_face]), elem_conn)
        if !isnothing(i_elem)
            for d in 1:dim
                idx = (i_elem-1)*dim + d
                f[idx] = N_face[i_face] * fₜ[d] * scaling
            end
        end
    end
    return f
end

function compute_surface_force_scalar(N_face, scaling, fₜ, face_conn, elem_conn)
    n_face_nodes = length(face_conn)
    f = zeros(length(elem_conn))
    for i_face in 1:n_face_nodes
        i_elem = findfirst(isequal(face_conn[i_face]), elem_conn)
        if !isnothing(i_elem)
            f[i_elem] = N_face[i_face] * fₜ * scaling
        end
    end
    return f
end
function thermal_expansion_force!(
    f::Vector{Float64},
    B::Vector{Matrix{Float64}},
    material::Material,
    jac_data,
    gauss_data,
    elem_id::Int)

    α = material.properties[:thermal_expansion]
    ΔT = get(material.properties, :delta_temperature, 1.0)
    C = material.tensors[1]
    dim = material.dim

    for qp in 1:length(gauss_data.weights)
        detJ, _ = jac_data[elem_id][qp]
        scaling = 0.5 * abs(detJ) * gauss_data.weights[qp]

        if isa(α, Number)
            # Isotropic: α * I (Voigt form)
            ε_th = fill(0.0, size(B[1], 1))
            if dim == 2
                ε_th[1] = α
                ε_th[2] = α
                ε_th[3] = 0.0
            elseif dim == 3
                ε_th[1] = α
                ε_th[2] = α
                ε_th[3] = α
                ε_th[4:end] .= 0.0
            else
                error("Unsupported dim = $dim")
            end
            ε_th .*= ΔT
        else
            # Anisotropic thermal strain
            ε_th = ΔT * α
        end

        f .= B[qp]' * C * ε_th * scaling
    end
end


# ============================================================
# Force Computation
# ============================================================

function compute_element_force(::Val{:elastic}, material::Material,
    jac_data, gauss_data,
    jac_fdata, gauss_fdata,
    forces,
    face_conn::Vector{Int},
    elem_conn,
    B_dict)

    # Extract element properties
    B = B_dict[:strain]
    n_dofs = size(B[1], 2)
    dim = material.dim
    n_nodes = n_dofs 
   
    f = zeros(n_dofs)

    # Volume force
    if !isnothing(forces.fᵥ)
        for qp in 1:length(gauss_data.weights)
            detJ, _ = jac_data[qp]
            scaling = 0.5 * abs(detJ) * gauss_data.weights[qp]
             
            if material.symmetry == :out_of_plane
                 
                f += compute_out_of_plane_force(B[qp], material.tensors[1], scaling)
            

            elseif haskey(material.properties, :thermal_expansion)
               f += thermal_expansion_force!(f, B, material, jac_data, gauss_data, elem_id)
            else
               
                N = gauss_data.N[qp]
                f += compute_volume_force_vector(N, scaling, forces.fᵥ, dim, n_nodes)
            end
        end
    end
    
    # Traction force
    if !isnothing(forces.fₜ)
       
        for qp in 1:length(gauss_fdata.weights)
            N_face = gauss_fdata.N[qp]
            detJ = jac_fdata[face_conn][qp][1]
            scaling = abs(detJ) * gauss_fdata.weights[qp]
            f += compute_surface_force_vector(N_face, scaling, forces.fₜ, dim, face_conn, elem_conn)
        end
    end
    
    return f
end

function compute_element_force(::Val{:thermal}, material::Material,
    jac_data, gauss_data,
    jac_fdata, gauss_fdata,
    forces,
    face_conn::Vector{Int},
    elem_conn,B_dict)

    # Extract element properties
    B = B_dict[:temperature_gradient]
    n_dofs = size(B[1], 2)
   
    f = zeros(n_dofs)
    
    # Volume heat source
    if !isnothing(forces.fᵥ)
        for qp in 1:length(gauss_data.weights)
            detJ, _ = jac_data[qp]
            scaling = 0.5 * abs(detJ) * gauss_data.weights[qp]
            
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
        for qp in 1:length(gauss_fdata.weights)
            N_face = gauss_fdata.N[qp]
            detJ = jac_fdata[face_conn][qp][1]
            scaling = abs(detJ) * gauss_fdata.weights[qp]
            f .+= compute_surface_force_scalar(N_face, scaling, forces.fₜ, face_conn, elem_conn)
        end
    end
    
    return f
end

function compute_element_force(::Val{:piezoelectric}, material::Material,
    jac_data, gauss_data,
    jac_fdata, gauss_fdata,
    forces,
    face_conn::Vector{Int},
    elem_conn,B_dict)

    # Extract element properties
    B_u = B_dict[:strain]
    B_e = B_dict[:electric_field]
    n_dofs_u = size(B_u, 2)
    n_dofs_e = size(B_e, 2)
    dim = material.dim
    n_nodes = n_dofs_u ÷ dim
   
    f = zeros(n_dofs_u + n_dofs_e)
    
    # Volume force - mechanical
    if !isnothing(forces.fᵥ)
        for qp in 1:length(gauss_data.weights)
            detJ, _ = jac_data[qp]
            scaling = 0.5 * abs(detJ) * gauss_data.weights[qp]
            
            if material.symmetry == :out_of_plane
                f_mech = compute_out_of_plane_force(B_u[qp], material.tensors[1], scaling)
                f .+= f_mech
            elseif haskey(material.properties, :thermal_expansion)
                f = thermal_expansion_force!(f, B, material, jac_data, gauss_data, elem_id)
            else
                N = gauss_data.N[qp]
                f_mech = compute_volume_force_vector(N, scaling, forces.fᵥ, dim, n_nodes)
                f[1:n_dofs_u] .+= f_mech
            end
        end
    end
    
    # Traction force - mechanical
    if !isnothing(forces.fₜ)
        for qp in 1:length(gauss_fdata.weights)
            N_face = gauss_fdata.N[qp]
            detJ = jac_fdata[face_conn][qp][1]
            scaling = abs(detJ) * gauss_fdata.weights[qp]
            f_mech = compute_surface_force_vector(N_face, scaling, forces.fₜ, dim, face_conn, elem_conn)
            f[1:n_dofs_u] .+= f_mech
        end
    end
    
    return f
end

function compute_element_force(::Val{:thermoelastic}, material::Material,
    jac_data, gauss_data,
    jac_fdata, gauss_fdata,
    forces,
    face_conn::Vector{Int},
    elem_conn,B_dict)

    # Extract element properties
    B_u = B_dict[:strain]
    B_t = B_dict[:temperature_gradient]
    n_dofs_u = size(B_u[1], 2)
    n_dofs_t = size(B_t[1], 2)
    dim = material.dim
    n_nodes = n_dofs_u ÷ dim
   
    f = zeros(n_dofs_u + n_dofs_t)
    
    # Volume force - mechanical
    if !isnothing(forces.fᵥ)
        for qp in 1:length(gauss_data.weights)
            detJ, _ = jac_data[qp]
            scaling = 0.5 * abs(detJ) * gauss_data.weights[qp]
            
            if material.symmetry == :out_of_plane
                f_mech = compute_out_of_plane_force(B_u[qp], material.tensors[1], scaling)
                f[1:n_dofs_u] += f_mech
            else
                N = gauss_data.N[qp]
                f_mech = compute_volume_force_vector(N, scaling, forces.fᵥ, dim, n_nodes)
                f[1:n_dofs_u] += f_mech
            end
        end
    end
    
    # Volume force - thermal
    if !isnothing(forces.fᵥ_thermal) && !iszero(forces.fᵥ_thermal)
        for qp in 1:length(gauss_data.weights)
            detJ, _ = jac_data[qp]
            scaling = 0.5 * abs(detJ) * gauss_data.weights[qp]
            
            if material.symmetry == :out_of_plane
                f_thermal = compute_out_of_plane_force(B_t[qp], material.tensors[3], scaling)
                f[n_dofs_u+1:end] += f_thermal
            else
                N = gauss_data.N[qp]
                f_thermal = compute_volume_force_scalar(N, scaling, forces.fᵥ_thermal)
                f[n_dofs_u+1:end] += f_thermal
            end
        end
    end
    
    # Traction force - mechanical
    if !isnothing(forces.fₜ)
        for qp in 1:length(gauss_fdata.weights)
            N_face = gauss_fdata.N[qp]
            detJ = jac_fdata[face_conn][qp][1]
            scaling = abs(detJ) * gauss_fdata.weights[qp]
            f_mech = compute_surface_force_vector(N_face, scaling, forces.fₜ, dim, face_conn, elem_conn)
            f[1:n_dofs_u] += f_mech
        end
    end
    
    # Surface heat flux
    if !isnothing(forces.fₜ_thermal) && !iszero(forces.fₜ_thermal)
        for qp in 1:length(gauss_fdata.weights)
            N_face = gauss_fdata.N[qp]
            detJ = jac_fdata[face_conn][qp][1]
            scaling = abs(detJ) * gauss_fdata.weights[qp]
            f_thermal = compute_surface_force_scalar(N_face, scaling, forces.fₜ_thermal, face_conn, elem_conn)
            f[n_dofs_u+1:end] += f_thermal
        end
    end
    
    return f
end

function compute_element_force(::Val{:poroelastic}, material::Material,
    jac_data, gauss_data,
    jac_fdata, gauss_fdata,
    forces,
    face_conn::Vector{Int},
    elem_conn,B_dict)

    # Extract element properties
    B_u = B_dict[:strain]
    n_dofs_u = size(B_u[1], 2)
    n_dofs_p = size(B_dict[:pressure][1], 2)
    dim = material.dim
    n_nodes = n_dofs_u ÷ dim
   
    f = zeros(n_dofs_u + n_dofs_p)
    
    # Volume force - mechanical
    if !isnothing(forces.fᵥ)
        for qp in 1:length(gauss_data.weights)
            detJ, _ = jac_data[qp]
            scaling = 0.5 * abs(detJ) * gauss_data.weights[qp]
            
            if material.symmetry == :out_of_plane
                f_mech = compute_out_of_plane_force(B_u[qp], material.tensors[1], scaling)
                f[1:n_dofs_u] += f_mech
            else
                N = gauss_data.N[qp]
                f_mech = compute_volume_force_vector(N, scaling, forces.fᵥ, dim, n_nodes)
                f[1:n_dofs_u] += f_mech
            end
        end
    end
    
    # Traction force - mechanical
    if !isnothing(forces.fₜ)
        for qp in 1:length(gauss_fdata.weights)
            N_face = gauss_fdata.N[qp]
            detJ = jac_fdata[face_conn][qp][1]
            scaling = abs(detJ) * gauss_fdata.weights[qp]
            f_mech = compute_surface_force_vector(N_face, scaling, forces.fₜ, dim, face_conn, elem_conn)
            f[1:n_dofs_u] += f_mech
        end
    end
    
    # Pressure flux (surface force)
    if !isnothing(forces.fₚ) && !iszero(forces.fₚ)
        for qp in 1:length(gauss_fdata.weights)
            N_face = gauss_fdata.N[qp]
            detJ = jac_fdata[face_conn][qp][1]
            scaling = abs(detJ) * gauss_fdata.weights[qp]
            f_pressure = compute_surface_force_scalar(N_face, scaling, forces.fₚ, face_conn, elem_conn)
            f[n_dofs_u+1:end] += f_pressure
        end
    end
    
    return f
end

function compute_element_force(::Val{:thermo_piezoelectric}, material::Material,
    jac_data, gauss_data,
    jac_fdata, gauss_fdata,
    forces,
    face_conn::Vector{Int},
    elem_conn,B_dict)

    # Extract element properties
    B_u = B_dict[:strain]
    n_dofs_u = size(B_u[1], 2)
    n_dofs_e = size(B_dict[:electric_field][1], 2)
    n_dofs_t = size(B_dict[:temperature_gradient][1], 2)
    dim = material.dim
    n_nodes = n_dofs_u ÷ dim
   
    f = zeros(n_dofs_u + n_dofs_e + n_dofs_t)
    
    # Volume force - mechanical
    if !isnothing(forces.fᵥ)
        for qp in 1:length(gauss_data.weights)
            detJ, _ = jac_data[qp]
            scaling = 0.5 * abs(detJ) * gauss_data.weights[qp]
            
            if material.symmetry == :out_of_plane
                f_mech = compute_out_of_plane_force(B_u[qp], material.tensors[1], scaling)
                f[1:n_dofs_u] += f_mech
            else
                N = gauss_data.N[qp]
                f_mech = compute_volume_force_vector(N, scaling, forces.fᵥ, dim, n_nodes)
                f[1:n_dofs_u] += f_mech
            end
        end
    end
    
    # Volume force - thermal
    if !isnothing(forces.fᵥ_thermal) && !iszero(forces.fᵥ_thermal)
        for qp in 1:length(gauss_data.weights)
            detJ, _ = jac_data[qp]
            scaling = 0.5 * abs(detJ) * gauss_data.weights[qp]
            
            if material.symmetry == :out_of_plane
                f_thermal = compute_out_of_plane_force(B_t[qp], material.tensors[5], scaling)
                f[n_dofs_u+n_dofs_e+1:end] += f_thermal
            else
                N = gauss_data.N[qp]
                f_thermal = compute_volume_force_scalar(N, scaling, forces.fᵥ_thermal)
                f[n_dofs_u+n_dofs_e+1:end] += f_thermal
            end
        end
    end
    
    # Traction force - mechanical
    if !isnothing(forces.fₜ)
        for qp in 1:length(gauss_fdata.weights)
            N_face = gauss_fdata.N[qp]
            detJ = jac_fdata[face_conn][qp][1]
            scaling = abs(detJ) * gauss_fdata.weights[qp]
            f_mech = compute_surface_force_vector(N_face, scaling, forces.fₜ, dim, face_conn, elem_conn)
            f[1:n_dofs_u] += f_mech
        end
    end
    
    # Surface heat flux
    if !isnothing(forces.fₜ_thermal) && !iszero(forces.fₜ_thermal)
        for qp in 1:length(gauss_fdata.weights)
            N_face = gauss_fdata.N[qp]
            detJ = jac_fdata[face_conn][qp][1]
            scaling = abs(detJ) * gauss_fdata.weights[qp]
            f_thermal = compute_surface_force_scalar(N_face, scaling, forces.fₜ_thermal, face_conn, elem_conn)
            f[n_dofs_u+n_dofs_e+1:end] += f_thermal
        end
    end
    
    return f
end

function compute_element_force(::Val{:viscoelastic}, material::Material,
    jac_data, gauss_data,
    jac_fdata, gauss_fdata,
    forces,
    face_conn::Vector{Int},
    elem_conn,B_dict)

    # Extract element properties
    B = B_dict[:strain]
    n_dofs = size(B[1], 2)
    dim = material.dim
    n_nodes = n_dofs ÷ dim
   
    f = zeros(n_dofs)
    
    # Volume force
    if !isnothing(forces.fᵥ)
        for qp in 1:length(gauss_data.weights)
            detJ, _ = jac_data[qp]
            scaling = 0.5 * abs(detJ) * gauss_data.weights[qp]
            
            if material.symmetry == :out_of_plane
                f += compute_out_of_plane_force(B[qp], material.tensors[1], scaling)
            else
                N = gauss_data.N[qp]
                f += compute_volume_force_vector(N, scaling, forces.fᵥ, dim, n_nodes)
            end
        end
    end
    
    # Traction force
    if !isnothing(forces.fₜ)
        for qp in 1:length(gauss_fdata.weights)
            N_face = gauss_fdata.N[qp]
            detJ = jac_fdata[face_conn][qp][1]
            scaling = abs(detJ) * gauss_fdata.weights[qp]
            f += compute_surface_force_vector(N_face, scaling, forces.fₜ, dim, face_conn, elem_conn)
        end
    end
    
    return f
end    


