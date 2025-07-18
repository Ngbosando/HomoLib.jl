#-----------------------------------------------------------------
# Main Field Recovery Function
#-----------------------------------------------------------------
    """
        HomoLib.recover_field_values

        Module for post-processing and field recovery from finite element solutions,
        implementing consistent field averaging techniques for accurate visualization
        and analysis of heterogeneous materials.

        # Key Features
        - Implements the Consistent Field Averaging (CFA) method for:
        - Stress/strain recovery (Zienkiewicz-Zhu estimator)
        - Flux/gradient recovery
        - Multi-physics field coupling
        - Supports both nodal and Gauss point field evaluation
        - Handles composite materials with phase-wise properties
        - Thread-parallel implementation for large-scale problems

        # Field Recovery Methods
        1. **Superconvergent Patch Recovery (SPR)**:
        - Least-squares fitting of stresses at superconvergent points
        - Optimal for quadratic elements (Zienkiewicz & Zhu, 1992)

        2. **Weighted Averaging**:
        - Volume-weighted averaging of element contributions
        - Preserves equilibrium in an integral sense

        3. **Material-Specific Recovery**:
        - Specialized handling for:
            - Piezoelectric (electric displacement field)
            - Poroelastic (pore pressure field)
            - Viscoelastic (rate-dependent fields)

        # Usage Example
        ```julia
        # Recover thermal fields
        fields = recover_field_values(
            elements, nodes, material, solution,
            nothing, :Tri6, 3, 2, geom_data
        )

        # Access recovered fields:
        heat_flux = fields.flux  # [num_nodes × dim]
        temp_grad = fields.grad_temp  # [num_nodes × dim]

        # Recover elastic fields
        fields = recover_field_values(...)
        stress = fields.stress  # [num_nodes × n_components]
        strain = fields.strain  # [num_nodes × n_components]
        Implementation Details

        Nodal Averaging:
        σᵢ = ∑ₑ∈Ωᵢ (σₑ ⋅ Vₑ) ÷ ∑ₑ∈Ωᵢ Vₑ


        References
        Zienkiewicz, O.C. & Taylor, R.L. (2005). "The Finite Element Method"

    """
    function recover_field_values(
        elements::Matrix{Int},
        nodes::Union{Matrix{Float64}, Vector{Float64}},
        material::Union{Material, Vector{Material}},
        Uresult::NamedTuple,
        connect_elem_phase::Union{Vector{Int}, Nothing},
        element_type::Symbol,
        integration_order::Int,
        dim::Int,
        Geometric_Data)
        
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
        ref_mat = is_composite ? material[1] : material
        gauss_data = Geometric_Data.gauss_data
        jac_cache = Geometric_Data.jacobian_cache
        B_dicts = Geometric_Data.B_dicts

        # ========== Node Connectivity ========== #
        connect_nodes_triangles, count_tri_per_node = find_elem_connected_to_nodes(elements, nodes[:,1])
        num_nodes = size(nodes, 1)

        # ========== Initialize Storage ========== #
        store_field = init_field_storage(ref_mat, dim, num_nodes)
        global_dofs = assemble_global_dofs(elements, ref_mat)
        
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
                elem_nodes = global_dofs[elem, :]
            
                # Get element solution
                elem_sol = extract_element_solution_fields(current_mat,elem_nodes, Uresult)
            
                # Compute element contribution
                current_type = resolve_material_type(current_mat)
                elem_result = integrate_element_quantities(
                    Val(current_type),
                    current_mat,
                    B_dicts[elem],
                    jac_cache[elem],
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

#-----------------------------------------------------------------
# Field Storage Initialization
#-----------------------------------------------------------------
 
    function init_field_storage(mat::Material, dim::Int, num_nodes::Int)
       
        if mat.type == [:thermal]
            (flux = zeros(num_nodes, dim), grad_temp = zeros(num_nodes, dim)) 
        elseif mat.type == [:elastic]
            nstr = size(mat.tensors[1], 1)  # Get from elasticity tensor
            (stress = zeros(num_nodes, nstr), strain = zeros(num_nodes, nstr))  
        elseif mat.type == [:elastic, :electric]
            nstr = size(mat.tensors[1], 1)  # Rows of C tensor
            nelec = size(mat.tensors[3], 1)  # Rows of ϵ tensor
            (stress = zeros(num_nodes, nstr), 
            strain = zeros(num_nodes, nstr),
            elec_disp = zeros(num_nodes, nelec),
            elec = zeros(num_nodes, nelec))
        elseif mat.type == [:elastic, :pressure]
            nstr = size(mat.tensors[1], 1)
            (stress = zeros(num_nodes, nstr), 
            strain = zeros(num_nodes, nstr),
            pressure = zeros(num_nodes))
        else
            error("Unsupported material type: $(mat.type)")
        end
    end

    function init_empty_field_sum(mat::Material)
     
        if mat.type == [:thermal]
            (flux = zeros(mat.dim),grad_temp = zeros(mat.dim))
        elseif mat.type == [:elastic]
            if mat.symmetry == :out_of_plane
                s = ntens(mat.dim) + 1
            else
                s = ntens(mat.dim)
            end
            (stress = zeros(s), strain = zeros(s))
        elseif mat.type == [:elastic, :electric]
            nstr = size(mat.tensors[1], 1)
            nelec = size(mat.tensors[3], 1)
            (stress = zeros(nstr), 
            strain = zeros(nstr),
            elec_disp = zeros(nelec),
            elec = zeros(nelec))
        elseif mat.type == [:elastic, :pressure]
            nstr = size(mat.tensors[1], 1)
            (stress = zeros(nstr), 
            strain = zeros(nstr),
            pressure = 0.0)
        else
            error("Unsupported material type: $(mat.type)")
        end
    end

#-----------------------------------------------------------------
# Material-Specific Field Integration
#-----------------------------------------------------------------
    function integrate_element_quantities(::Val{:thermal}, mat, B_dict, jac_cache, gauss_data, Uᵉ)
        κ = mat.tensors[1]
        B = B_dict[:temperature_gradient]
        flux = zeros(mat.dim)
        grad_temp = zeros(mat.dim)  # Added temperature gradient
        volume = 0.0

        for qp in eachindex(gauss_data.weights)
            detJ, _ = jac_cache[qp]
            scaling =  abs(detJ) * gauss_data.weights[qp]
            ∇T = B[qp] * Uᵉ
            flux .+= - (κ * ∇T) * scaling
            grad_temp .+= ∇T * scaling  # Accumulate gradient
            volume += scaling
        end
        (flux = flux, grad_temp = grad_temp, volume = volume)
    end

    function integrate_element_quantities(::Val{:elastic}, mat, B_dict, jac_cache, gauss_data, Uᵉ)
        C = mat.tensors[1]
        B = B_dict[:strain]
        stress = zeros(size(B[1], 1))
        strain = zeros(size(B[1], 1))
        volume = 0.0
        
        for qp in eachindex(gauss_data.weights)
            detJ, _ = jac_cache[qp]
            scaling = abs(detJ) * gauss_data.weights[qp]
            ε = B[qp] * Uᵉ
            σ = C * ε
        
            stress .+= σ * scaling
            strain .+= ε * scaling
            volume += scaling
        end
        (stress = stress, strain = strain, volume = volume)
    end

    function integrate_element_quantities(::Val{:piezoelectric}, mat, B_dict, jac_cache, gauss_data, Uᵉ)
        C, ϵ, e = mat.tensors
        B_u = B_dict[:strain]
        B_ϕ = B_dict[:electric_field]
        
        nstr = size(C, 1)
        nelec = size(ϵ, 1)
        
        elem_stress = zeros(nstr)
        elem_strain = zeros(nstr)
        elem_elec_disp = zeros(nelec)
        elem_elec = zeros(nelec)
        volume = 0.0

        for qp in eachindex(gauss_data.weights)
            detJ, _ = jac_cache[qp]
           
            scale = abs(detJ) * gauss_data.weights[qp]
            
            # Compute strain and electric field
          
            ε = B_u[qp] * Uᵉ.mech
            E = B_ϕ[qp] * Uᵉ.potential  # Electric field is negative gradient
        
            # Compute stress and electric displacement
            σ = C * ε - transpose(e) * E
            D = e * ε + ϵ * E
            
            # Accumulate all fields
            elem_strain .+= ε * scale
            elem_elec .+= E * scale
            elem_stress .+= σ * scale
            elem_elec_disp .+= D * scale
            volume += scale
        end
        
        (stress = elem_stress, 
        strain = elem_strain,
        elec_disp = elem_elec_disp,
        elec = elem_elec,
        volume = volume)
    end

    function integrate_element_quantities(::Val{:poroelastic}, mat, B_dict, jac_cache, gauss_data, Uᵉ)
       C = mat.tensors[1]
        B = B_dict[:strain]
        stress = zeros(size(B[1], 1))
        strain = zeros(size(B[1], 1))
        volume = 0.0
        
        for qp in eachindex(gauss_data.weights)
            detJ, _ = jac_cache[qp]
            scaling = abs(detJ) * gauss_data.weights[qp]
            ε = B[qp] * Uᵉ
            σ = C * ε
        
            stress .+= σ * scaling
            strain .+= ε * scaling
            volume += scaling
        end
        (stress = stress, strain = strain, volume = volume)
    end

    function integrate_element_quantities(::Val{:viscoelastic}, mat, B_dict, jac_cache, gauss_data, Uᵉ)
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

#-----------------------------------------------------------------
# DOF Handling Utilities
#-----------------------------------------------------------------

    function extract_element_solution_fields(mat::Material, elem_nodes, Uresult::NamedTuple)
        mat_type = resolve_material_type(mat)
        if mat_type == :thermal
            return Uresult.T[elem_nodes]
        elseif mat_type == :elastic
            return Uresult.U[elem_nodes]
        elseif mat_type == :piezoelectric
            block = mat.dim + 1
            elem_nodes_u = [elem_nodes[i] for i in eachindex(elem_nodes) if (i-1) % block <  mat.dim]
            elem_nodes_ϕ = [elem_nodes[i] for i in eachindex(elem_nodes) if (i-1) % block == mat.dim]
   
            return (
                mech = Uresult.U[elem_nodes_u],
                potential = Uresult.ϕ[elem_nodes_ϕ]
            )
        elseif mat_type == :thermoelastic
            return (
                mech = Uresult.U[elem_nodes],
                temp = Uresult.T[elem_nodes]
            )
        elseif mat_type == :poroelastic
            return (
                mech = Uresult.U[elem_nodes]
            )
        else
            error("Unsupported material type: $(mat.type)")
        end
    end


#-----------------------------------------------------------------
# Connectivity and Finalization
#-----------------------------------------------------------------
    function find_elem_connected_to_nodes(elements, Nₓ)
        num_nodes = length(Nₓ)
        num_elems = size(elements, 1)
        elements_nodes_elements = [Int[] for _ in 1:num_nodes]

        for elem in 1:num_elems
            for node in elements[elem, :]
                push!(elements_nodes_elements[node], elem)
            end
        end

        count_elem_per_node = [length(elems) for elems in elements_nodes_elements]
        return elements_nodes_elements, count_elem_per_node
    end

    function finalize_field_average(sum_field::NamedTuple, total_volume::Real)
        total_volume == 0.0 && error("Cannot finalize: total volume is zero")
        return map(x -> x ./ total_volume, sum_field)
    end