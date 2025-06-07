
include("compute_flux.jl")
function compute_effective_property(
    materials::Union{Material, Vector{Material}},
    elements::Matrix{Int},
    nodes::Matrix{Float64},
    connect_elem_phase::Union{Vector{Int}, Nothing},
    solver_results::NamedTuple,
    element_type::Symbol,
    integration_order::Int,
    dim::Int)
    # ========== Phase Handling ========== #
    is_composite = materials isa Vector{Material}
    if is_composite
        connect_elem_phase === nothing && error("connect_elem_phase required for composite materials")
        length(connect_elem_phase) == size(elements, 1) || error("connect_elem_phase length mismatch")

        unique_phases = sort(unique(connect_elem_phase))
        length(materials) == length(unique_phases) || error("Material count mismatch with phases")

        phase_to_idx = Dict(zip(unique_phases, 1:length(unique_phases)))

        # Validate all materials are consistent
        ref_dofs = get_dofs_per_node(materials[1])
        ref_B_types = materials[1].B_types
        for mat in materials
            get_dofs_per_node(mat) == ref_dofs || error("Mismatch in dofs_per_node")
            mat.B_types == ref_B_types || error("Mismatch in B_types")
        end
        dofs_per_node = ref_dofs
    else
        dofs_per_node = get_dofs_per_node(materials)
    end

    # ========== Precomputation ========== #
    gauss_data = shape_data(element_type, integration_order, dim)
    jac_cache = jacobian_data(elements, nodes, gauss_data)
    ref_mat = is_composite ? materials[1] : materials
    B_dicts = build_B_matrices(nodes, elements, ref_mat, gauss_data, jac_cache)

    # ========== Init Storage ========== #
    result_template = init_global_storage(materials, dim)
    result_acc = [deepcopy(result_template) for _ in 1:Threads.nthreads()]
    vol_acc = zeros(Threads.nthreads())
    global_dofs = assemble_global_dofs(elements, ref_mat)

    # ========== Main Loop ========== #
    @threads for e in 1:size(elements, 1)
        tid = Threads.threadid()
        mat = is_composite ? materials[phase_to_idx[connect_elem_phase[e]]] : materials
        elem_nodes = global_dofs[e, :]
 
        elem_sol = init_element_solution(mat, elem_nodes, dim)
       
        elem_sol = get_element_solution!(elem_sol, mat, elem_nodes, solver_results)
        mat_type = resolve_material_type(mat.type)
        elem_result = compute_element_contribution(
            Val(mat_type), mat, B_dicts[e], e, jac_cache, gauss_data, elem_sol
        )

        accumulate_results!(result_acc[tid], elem_result)
        vol_acc[tid] += elem_result.volume
    end

    # ========== Finalization ========== #
    total_volume = sum(vol_acc)
    final_result = finalize_results(result_acc, total_volume)

    return final_result
end

#-----------------------------------------------------------------
# Material-Specific Implementations
#-----------------------------------------------------------------

function init_element_solution(mat::Material, elem_nodes, dim::Int)
    node_count = length(elem_nodes)
    lc = get_load_case_count(mat, dim)
    mat_type = resolve_material_type(mat.type)
    if mat_type == :thermal
        return zeros(node_count, lc)

    elseif mat_type == :elastic || mat_type == :viscoelastic
        return zeros(node_count, lc)

    elseif mat_type == :piezoelectric
        return (
            mech = zeros(node_count, lc),
            elec = zeros(node_count, dim)
        )

    elseif mat_type == :thermoelastic
        return (
            mech = zeros(node_count, lc),
            temp = zeros(node_count, dim)
        )

    elseif mat_type == :poroelastic
        return (
            mech = zeros(node_count, lc),
            pres = zeros(node_count)
        )

    else
        error("Unsupported material type: $(mat_type)")
    end
end

function get_element_solution!(elem_sol, mat::Material, elem_nodes, solver_results)
    Ns = length(solver_results.U)
    mat_type = resolve_material_type(mat.type)
    if mat_type in (:thermal, :elastic, :viscoelastic)
        for i in 1:Ns
            cache = solver_results.U[i]
            elem_sol[:, i] = cache[elem_nodes]
        end

    elseif mat_type == :piezoelectric
        for i in 1:Ns
            elem_sol.mech[:, i] = solver_results.U[i][elem_nodes, 1]
            elem_sol.elec[:, i] = solver_results.V[i][elem_nodes]
        end

    elseif mat_type == :thermoelastic
        for i in 1:Ns
            elem_sol.mech[:, i] = solver_results.U[i][elem_nodes, 1]
            elem_sol.temp[:, i] = solver_results.V[i][elem_nodes]
        end

    elseif mat_type == :poroelastic
        for i in 1:Ns
            elem_sol.mech[:, i] = solver_results.U[i][elem_nodes, 1]
            elem_sol.pres[:] = solver_results.V[i][elem_nodes]  # Pres assumed scalar per node
        end

    else
        error("Unsupported material type: $(mat_type)")
    end

    return elem_sol
end


function compute_element_contribution(::Val{:thermal}, mat::Material, B_dict, elem_id, jacobian_data, gauss_data, Uᵉ)
    κ = mat.tensors[1]
    B = B_dict[:temperature_gradient]
    K_eff = zero(κ)
    volume = 0.0

    for qp in eachindex(gauss_data.weights)
        detJ, _ = jacobian_data[elem_id][qp]
        scale =  0.5 * abs(detJ) * gauss_data.weights[qp]
        ∇T = B[qp] * Uᵉ

        q = κ * ∇T

        K_eff += q  * scale

        volume += scale 
    end

    return (K = K_eff, volume = volume)
end

function compute_element_contribution(::Val{:elastic}, mat::Material, B_dict, elem_id, jacobian_data, gauss_data, Uᵉ)
    C = mat.tensors[1]
    B = B_dict[:strain]
 
    C_eff = zero(C)
    volume = 0.0

    for qp in eachindex(gauss_data.weights)
        detJ, _ = jacobian_data[elem_id][qp]
        scale = 0.5 * abs(detJ) * gauss_data.weights[qp]
       
        if mat.symmetry == :out_of_plane
            M = zero(C); M[4, 4] = 1.0
            ε = B[qp] * Uᵉ + M
            
        else
       
            ε = B[qp] * Uᵉ
         
        end
        σ = C * ε
        C_eff += σ  * scale
        volume += scale 
    end

    return (C = C_eff, volume = volume)
end

function compute_element_contribution(::Val{:piezoelectric}, mat::Material, B_dict, elem_id, jacobian_data, gauss_data, Uᵉ)
    C, e, ϵ = mat.tensors
    B_u = B_dict[:strain]
    B_e = B_dict[:electric_field]

    C_eff = zero(C)
    e_eff = zero(e)
    ϵ_eff = zero(ϵ)
    volume = 0.0

    for qp in eachindex(gauss_data.weights)
        detJ, _ = jacobian_data[elem_id][qp]
        scale = abs(detJ) * gauss_data.weights[qp]

        ε = B_u[qp] * Uᵉ.mech
        E = B_e[qp] * Uᵉ.elec
        σ = C * ε - e' * E
        D = e * ε + ϵ * E

        C_eff += σ  * scale
        e_eff += D  * scale
        ϵ_eff += E  * scale
        volume += scale
    end

    return (C = C_eff, e = e_eff, ϵ = ϵ_eff, volume = volume)
end

function compute_element_contribution(::Val{:thermoelastic}, mat::Material, B_dict, elem_id, jacobian_data, gauss_data, Uᵉ)
    C, β, κ = mat.tensors
    B_u = B_dict[:strain]
    B_t = B_dict[:temperature_gradient]

    C_eff = zero(C)
    β_eff = zero(β)
    κ_eff = zero(κ)
    volume = 0.0

    for qp in eachindex(gauss_data.weights)
        detJ, _ = jacobian_data[elem_id][qp]
        scale = abs(detJ) * gauss_data.weights[qp]

        ε = B_u[qp] * Uᵉ.mech
        ∇T = B_t[qp] * Uᵉ.temp
        σ = C * ε - β * ∇T
        q = κ * ∇T

        C_eff += σ  * scale
        β_eff += σ  * scale
        κ_eff += q  * scale
        volume += scale
    end

    return (C = C_eff, β = β_eff, κ = κ_eff, volume = volume)
end

function compute_element_contribution(::Val{:poroelastic}, mat::Material, B_dict, elem_id, jacobian_data, gauss_data, Uᵉ)
    C, α, M = mat.tensors
    B_u = B_dict[:strain]
    B_p = B_dict[:pressure]

    C_eff = zero(C)
    α_eff = zero(α)
    M_eff = 0.0
    volume = 0.0

    for qp in eachindex(gauss_data.weights)
        detJ, _ = jacobian_data[elem_id][qp]
        scale = abs(detJ) * gauss_data.weights[qp]

        ε = B_u[qp] * Uᵉ.mech
        p = B_p[qp] * Uᵉ.pres
        σ = C * ε - α * p
        ζ = α' * ε + inv(M) * p

        C_eff += σ  * scale
        α_eff += σ  * scale
        M_eff += ζ  * scale
        volume += scale
    end

    return (C = C_eff, α = α_eff, M = M_eff, volume = volume)
end

function compute_element_contribution(::Val{:viscoelastic}, mat::Material, B_dict, elem_id, jacobian_data, gauss_data, Uᵉ)
    C, τ = mat.tensors  # τ unused here unless you model time
    B = B_dict[:strain]

    C_eff = zero(C)
    volume = 0.0

    for qp in eachindex(gauss_data.weights)
        detJ, _ = jacobian_data[elem_id][qp]
        scale = abs(detJ) * gauss_data.weights[qp]

        ε = B[qp] * Uᵉ
        σ = C * ε

        C_eff += σ  * scale
        volume += scale
    end

    return (C = C_eff, volume = volume)
end


#-----------------------------------------------------------------
# Helper Functions
#-----------------------------------------------------------------

function get_load_case_count(mat::Material, dim::Int)
     mat_type = resolve_material_type(mat.type)
    if mat_type == :thermal
        return dim  

    elseif mat_type in (:elastic, :thermoelastic, :piezoelectric, :poroelastic, :viscoelastic)
        if mat.symmetry == :isotropic
            return div(dim * (dim + 1), 2)  # 2D:3, 3D:6
        elseif mat.symmetry == :transversely_isotropic
            return dim == 2 ? 4 : 6
        elseif mat.symmetry == :orthotropic
            return dim == 2 ? 3 : 6
        elseif mat.symmetry == :out_of_plane
            return 4
        else
            error("Unsupported symmetry: $(mat.symmetry)")
        end

    else
        error("Unsupported material type: $(mat_type)")
    end
end
function ntens(dim)
    dim == 2 ? 3 : 6
end
# Corrected init_global_storage (remove volume from NamedTuple)
function init_global_storage(materials::Union{Material,Vector{Material}}, dim::Int)
    mat = materials isa Vector ? materials[1] : materials
    nstr = get_load_case_count(mat, dim)
    mat_type = resolve_material_type(mat.type)
    if mat_type == :thermal
        return (K = zeros(dim, dim),)

    elseif mat_type == :elastic
        return (C = zeros(nstr, nstr),)

    elseif mat_type == :piezoelectric
        return (
            C = zeros(nstr, nstr),
            e = zeros(dim, nstr),
            ϵ = zeros(dim, dim)
        )

    elseif mat_type == :thermoelastic
        return (
            C = zeros(nstr, nstr),
            β = zeros(nstr, dim),
            κ = zeros(dim, dim)
        )

    elseif mat_type == :poroelastic
        return (
            C = zeros(nstr, nstr),
            α = zeros(nstr),
            M = 0.0
        )

    else
        error("Unsupported material type: $(mat_type)")
    end
end
# Corrected accumulate_results! (skip volume)
function accumulate_results!(acc, elem_result)
    for field in propertynames(elem_result)
        if field != :volume && hasproperty(acc, field)
            acc_field = getproperty(acc, field)
            elem_field = getproperty(elem_result, field)
            acc_field .+= elem_field
        end
    end
   
end
# Corrected finalize_results (averages without volume)
function finalize_results(results, total_volume)
    final = deepcopy(results[1])
    for field in propertynames(final)
        final_field = getproperty(final, field)
        sum_field = sum(r -> getproperty(r, field), results)
        final_field .= sum_field ./ total_volume
    end
    return final
end