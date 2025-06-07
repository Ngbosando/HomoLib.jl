


function recover_field_values(
    elements::Matrix{Int},
    nodes::Matrix{Float64},
    material::Union{Material, Vector{Material}},
    Uresult::NamedTuple,
    connect_elem_phase::Union{Vector{Int}, Nothing},
    element_type::Symbol,
    integration_order::Int,
    dim::Int)
    # ========== Material Setup ========== #
     is_composite = material isa Vector{Material}
    if is_composite
        connect_elem_phase === nothing && error("connect_elem_phase required for composite materials")
        length(connect_elem_phase) == size(elements, 1) || error("connect_elem_phase length mismatch")

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


    # ========== Precomputation ========== #
    gauss_data = shape_data(element_type, integration_order, dim)
    jac_cache = jacobian_data(elements, nodes, gauss_data)
    B_dicts = build_B_matrices(nodes, elements, is_composite ? material[1] : material, gauss_data, jac_cache)

    # ========== Node Connectivity ========== #
    connect_nodes_triangles, count_tri_per_node = find_elem_connected_to_nodes(elements, nodes[:,1])
    num_nodes = size(nodes, 1)

    # ========== Initialize Storage ========== #

    store_field = init_field_storage(is_composite ? material[1] : material, dim,num_nodes)
 
    global_dofs = assemble_global_dofs(elements, material isa Vector ? material[1] : material)
 
    # ========== Main Node Loop ========== #

    for i in 1:num_nodes
        nb_elem = count_tri_per_node[i]
        connected_elements = connect_nodes_triangles[i]
      
        mat_ref = is_composite ? material[1] : material
        sum_field = init_empty_field_sum(mat_ref)
        sum_vol = 0.0


        for j in 1:nb_elem
            # Get material properties
            elem = connected_elements[j]
            current_mat = is_composite ? material[phase_to_idx[connect_elem_phase[elem]]] : material
            
            # Get element solution
            vec_assembly = global_dofs[elem, :]
            elem_sol = extract_element_solution_fields(current_mat, vec_assembly, Uresult)
          

            # Compute element contribution
            current_type = resolve_material_type(current_mat.type)
            elem_result = integrate_element_quantities(
                Val(current_type),
                current_mat,
                B_dicts[elem],
                elem,
                jac_cache,
                gauss_data,
                elem_sol
            )
            
          
            accumulate_results!(sum_field, elem_result)
            sum_vol += elem_result.volume
        end
       
            final_field = finalize_field_average(sum_field, sum_vol)
           
            for name in propertynames(store_field)
                target = getproperty(store_field, name)
                source = getproperty(final_field, name)
               
                # Handle both vector and scalar results
                if ndims(target) == 2
                    target[i, :] .= vec(source)
                else
                    target[i] = source
                end
            end
           
        
    end
    
   
    return store_field
end

function init_field_storage(mat::Material, dim::Int, num_nodes::Int)
    mat_type = resolve_material_type(mat.type)
    if mat_type == :thermal
        (flux = zeros(num_nodes, dim),)  # NamedTuple with flux field
    elseif mat_type == :elastic
        if mat.symmetry == :out_of_plane
        (stress = zeros(num_nodes, ntens(dim) + 1), strain = zeros(num_nodes, ntens(dim) + 1))
        else
            (stress = zeros(num_nodes, ntens(dim)), strain = zeros(num_nodes, ntens(dim)))
        end
    elseif mat_type == :piezoelectric
        (stress =stress = zeros(num_nodes, ntens(dim)), strain = zeros(num_nodes, ntens(dim)), electric = zeros(num_nodes, dim))
    elseif mat_type == :thermoelastic
        (stress =stress = zeros(num_nodes, ntens(dim)), strain = zeros(num_nodes, ntens(dim)), heatflux = zeros(num_nodes, dim))
    elseif mat_type == :poroelastic
        (stress =stress = zeros(num_nodes, ntens(dim)), strain = zeros(num_nodes, ntens(dim)), pressure = zeros(num_nodes, 1))
    else
        error("Unsupported material type: $(mat.type)")
    end
end
function init_empty_field_sum(mat::Material)
    mat_type = resolve_material_type(mat.type)
    if mat_type == :thermal
        (flux = zeros(mat.dim),)
    elseif mat_type == :elastic
        if mat.symmetry == :out_of_plane
            s = ntens(mat.dim) + 1
        else
            s = ntens(mat.dim)
        end
        (stress = zeros(s), strain = zeros(s))
    elseif mat_type == :piezoelectric
        (stress = zeros(ntens(mat.dim)), electric = zeros(mat.dim))
    elseif mat_type == :thermoelastic
        (stress = zeros(ntens(mat.dim)), heatflux = zeros(mat.dim))
    elseif mat_type == :poroelastic
        (stress = zeros(ntens(mat.dim)), pressure = zeros(1))
    elseif mat_type == :viscoelastic
        (stress = zeros(ntens(mat.dim)),)
    else
        error("Unsupported material type: $(mat.type)")
    end
end

function integrate_element_quantities(::Val{:thermal}, mat, B_dict, elem, jac_cache, gauss_data, Uᵉ)
    κ = mat.tensors[1]
    B = B_dict[:temperature_gradient]
    flux = zeros(mat.dim)
    volume = 0.0

    for qp in eachindex(gauss_data.weights)
        detJ, _ = jac_cache[elem][qp]
        scaling = 0.5 * abs(detJ) * gauss_data.weights[qp]
        ∇T = B[qp] * Uᵉ
        flux .+= - (κ * ∇T) * scaling
        volume += scaling
    end
    (flux = flux, volume = volume)
end

function integrate_element_quantities(::Val{:elastic}, mat, B_dict, elem, jac_cache, gauss_data, Uᵉ)
    C = mat.tensors[1]
    B = B_dict[:strain]
    stress = zeros(size(B[1], 1))
    strain = zeros(size(B[1], 1))
    volume = 0.0
    
    for qp in eachindex(gauss_data.weights)
        detJ, _ = jac_cache[elem][qp]
        scaling = 0.5 * abs(detJ) * gauss_data.weights[qp]
        ε = B[qp] * Uᵉ
        σ = C * ε
      
        stress .+= σ * scaling
        strain .+= ε * scaling
   
      
        volume += scaling
    end
    (stress = stress, strain = strain, volume = volume)
end

function integrate_element_quantities(::Val{:piezoelectric}, mat, B_dict, e, jac_cache, gauss_data, Uᵉ)
    C, e_tensor, ϵ = mat.tensors
    B_u = B_dict[:strain]
    B_ϕ = B_dict[:electric_field]
    
    elem_stress = zero(C)
    elem_electric = zero(ϵ)
    volume = 0.0

    for qp in eachindex(gauss_data.weights)
        detJ, _ = jac_cache[e][qp]
        scale = abs(detJ) * gauss_data.weights[qp]
        
        ε = B_u[qp] * Uᵉ.mech
        E = -B_ϕ[qp] * Uᵉ.elec  # Electric field is negative gradient
        
        σ = C * ε + e_tensor' * E
        D = e_tensor * ε - ϵ * E
        
        elem_stress += σ * scale
        elem_electric += D * scale
        volume += scale
    end
    
    (stress = elem_stress, electric = elem_electric, volume = volume)
end

function integrate_element_quantities(::Val{:thermoelastic}, mat, B_dict, e, jac_cache, gauss_data, Uᵉ)
    C, β, κ = mat.tensors
    B_u = B_dict[:strain]
    B_t = B_dict[:temperature_gradient]
    
    elem_stress = zero(C)
    elem_heatflux = zero(κ)
    volume = 0.0

    for qp in eachindex(gauss_data.weights)
        detJ, _ = jac_cache[e][qp]
        scale = abs(detJ) * gauss_data.weights[qp]
        
        ε = B_u[qp] * Uᵉ.mech
        ∇T = B_t[qp] * Uᵉ.temp
        
        σ = C * (ε - β * ∇T)
        q = -κ * ∇T
        
        elem_stress += σ * scale
        elem_heatflux += q * scale
        volume += scale
    end
    
    (stress = elem_stress, heatflux = elem_heatflux, volume = volume)
end

function integrate_element_quantities(::Val{:poroelastic}, mat, B_dict, e, jac_cache, gauss_data, Uᵉ)
    C, α, M = mat.tensors
    B_u = B_dict[:strain]
    B_p = B_dict[:pressure]
    
    elem_stress = zero(C)
    elem_pressure = 0.0
    volume = 0.0

    for qp in eachindex(gauss_data.weights)
        detJ, _ = jac_cache[e][qp]
        scale = abs(detJ) * gauss_data.weights[qp]
        
        ε = B_u[qp] * Uᵉ.mech
        p = B_p[qp] * Uᵉ.pres
        
        σ = C * ε - α * p
        elem_stress += σ * scale
        elem_pressure += p * scale
        volume += scale
    end
    
    (stress = elem_stress, pressure = elem_pressure, volume = volume)
end

function integrate_element_quantities(::Val{:viscoelastic}, mat, B_dict, e, jac_cache, gauss_data, Uᵉ)
    C, η = mat.tensors  # Elasticity and viscosity tensors
    B = B_dict[:strain]
    
    elem_stress = zero(C)
    volume = 0.0

    for qp in eachindex(gauss_data.weights)
        detJ, _ = jac_cache[e][qp]
        scale = abs(detJ) * gauss_data.weights[qp]
        
        ε = B[qp] * Uᵉ.mech
        ε̇ = B[qp] * Uᵉ.mech_rate  # Assuming rate field exists
        
        σ = C * ε + η * ε̇
        elem_stress += σ * scale
        volume += scale
    end
    
    (stress = elem_stress, volume = volume)
end

function extract_element_solution_fields(mat::Material, vec_assembly, Uresult)
    mat_type = resolve_material_type(mat.type)
    if mat_type == :thermal
        Uresult.T[vec_assembly]
    elseif mat_type == :elastic
        Uresult.U[vec_assembly]
    elseif mat_type == :piezoelectric
        (mech = Uresult.U[vec_assembly], potential = Uresult.ϕ[vec_assembly])
    elseif mat_type == :thermoelastic
        (mech = Uresult.U[vec_assembly], temp = Uresult.T[vec_assembly])
    else
        error("Unsupported material type: $(mat.type)")
    end
end
function find_elem_connected_to_nodes(elements, Nₓ)
    # Number of nodes is taken from the length of Nₓ.
    num_nodes = length(Nₓ)
    num_elems = size(elements, 1)

    # Create a vector of empty integer vectors.
    # Each entry corresponds to a node and will hold the indices of elements elementsed to that node.
    elements_nodes_elements = [Int[] for _ in 1:num_nodes]

    # Loop over each element.
    # Each row in 'elements' represents an element (triangle, quadrilateral, etc.)
    # with a variable number of node indices.
    for elem in 1:num_elems
        # Loop over all node indices in the current element.
        # This works regardless of how many nodes per element there are.
        for node in elements[elem, :]
            push!(elements_nodes_elements[node], elem)
        end
    end

    # Create a count vector that holds the number of elements elementsed to each node.
    count_elem_per_node = [length(elems) for elems in elements_nodes_elements]

    return elements_nodes_elements, count_elem_per_node
end
function finalize_field_average(sum_field::NamedTuple, total_volume::Real)
    total_volume == 0.0 && error("Cannot finalize: total volume is zero")
    return map(x -> x ./ total_volume, sum_field)
end


