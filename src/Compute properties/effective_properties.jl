include("compute_flux.jl")
#-----------------------------------------------------------------
# Main Function
#-----------------------------------------------------------------
    function compute_effective_property(
        materials::Union{Material, Vector{Material}},
        elements::Matrix{Int},
        nodes::Matrix{Float64},
        connect_elem_phase::Union{Vector{Int}, Nothing},
        solver_results::NamedTuple,
        element_type::Symbol,
        integration_order::Int,
        dim::Int,
        Geometric_Data)
        
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
            gauss_data = Geometric_Data.gauss_data
            jac_cache = Geometric_Data.jacobian_cache
            B_dicts = Geometric_Data.B_dicts
            ref_mat = is_composite ? materials[1] : materials
        

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
            
            mat_type = resolve_material_type(mat)
            elem_result = compute_element_contribution(
                Val(mat_type), mat, B_dicts[e], e, jac_cache, gauss_data, elem_sol
            )

            accumulate_results!(result_acc[tid], elem_result)
            vol_acc[tid] += elem_result.volume
            
        end

        # ========== Finalization ========== #
            total_volume = sum(vol_acc)
            final_result = finalize_results(result_acc, total_volume)

        return final_result,total_volume
    end
#-----------------------------------------------------------------
# Solution Initialization and Retrieval
#-----------------------------------------------------------------
    function init_element_solution(mat::Material, elem_nodes, dim::Int)
        n_nodes = length(elem_nodes)
        

        if mat.type in ([:thermal], [:electric])
            lc = get_load_case_count(mat, dim)
            return zeros(n_nodes, lc)
        elseif mat.type == [:elastic] && !haskey(mat.properties, :α_p)
            if mat.symmetry == :out_of_plane && dim == 2
                nstr = 4  # ε11, ε22, ε12, ε33
            else
                nstr = get_load_case_count(mat, dim)
            end
            return zeros(n_nodes, nstr)
        # Handle coupled physics - piezoelectric case
        elseif mat.type == [:elastic, :electric]
            # Determine load case counts based on symmetry
            if mat.symmetry == :out_of_plane && dim == 2
                nstr = 4  # ε11, ε22, ε12, ε33
                nelec = 3  # E1, E2, E3
            else
                nstr = get_load_case_count(Material(material.type, mat.dim, mat.symmetry, mat.properties, mat.tensors, mat.B_types,mat.mass_properties),dim)
                nelec = dim
            end
            nodes = Int(n_nodes / (dim + 1))
         
            return (
                U_u = zeros(dim * nodes, nstr),  # Displacement for ε̄ cases
                V_u = zeros(dim * nodes, nelec),  # Displacement for Ē cases
                U_ϕ = zeros(nodes, nstr),  # Potential for ε̄ cases
                V_ϕ = zeros(nodes, nelec)  # Potential for Ē cases
            )
        elseif mat.type == [:elastic, :pressure] || haskey(mat.properties, :α_p)
            # Determine load case counts based on symmetry
            if mat.symmetry == :out_of_plane && dim == 2
                nstr = 4  # ε11, ε22, ε12, ε33
                np = 3  # E1, E2, E3
            elseif haskey(mat.properties, :α_p) && :pressure ∉ mat.type
                nstr = get_load_case_count(mat, dim) + 1
                return zeros(n_nodes, nstr)

            else
                nstr = get_load_case_count(Material(material.type, mat.dim, mat.symmetry, mat.properties, mat.tensors, mat.B_types,mat.mass_properties),dim)
                np = dim
                nodes = Int(n_nodes / (dim + 1))

                return (
                U_u = zeros(dim * nodes, nstr),  # Displacement 
                U_p = zeros(nodes, nstr),  # pressure 
            )   
            end
            
        
        # Add other coupled physics here as needed
        else
            error("Unsupported material combination: $(mat.type)")
        end
    end

    function get_element_solution!(elem_sol, mat::Material, elem_nodes, solver_results)
        # Handle single physics materials
        if length(mat.type) == 1
            mat_type = mat.type[1]
            if mat_type in (:thermal, :elastic, :viscoelastic)
                Ns = length(solver_results.U)
                for i in 1:Ns
                    elem_sol[:, i] = solver_results.U[i][elem_nodes]
                end
                return elem_sol
            end
        
        # Handle coupled physics - piezoelectric case
        elseif mat.type == [:elastic, :electric]
            nstr = size(elem_sol.U_u, 2)
            nelec = size(elem_sol.V_ϕ, 2)

            block = mat.dim + 1
            elem_nodes_u = [elem_nodes[i] for i in eachindex(elem_nodes) if (i-1) % block <  mat.dim]
            elem_nodes_ϕ = [elem_nodes[i] for i in eachindex(elem_nodes) if (i-1) % block == mat.dim]
            # Retrieve mechanical solutions
  
            for i in 1:nstr
                elem_sol.U_u[:, i] = solver_results.U_mech[i][elem_nodes_u]
                elem_sol.U_ϕ[:, i] = solver_results.V_elec[i][elem_nodes_ϕ]
            end
            
            # Retrieve electrical solutions
            for i in 1:nelec
                elem_sol.V_u[:, i] = solver_results.U_mech[nstr + i][elem_nodes_u]
                elem_sol.V_ϕ[:, i] = solver_results.V_elec[nstr + i][elem_nodes_ϕ]
            end
            return elem_sol
         elseif mat.type == [:elastic, :pressure] || haskey(mat.properties, :α_p)   
                Ns = length(solver_results.U)
                for i in 1:Ns
                    elem_sol[:, i] = solver_results.U[i][elem_nodes]
                end
                return elem_sol
        end
    end

#-----------------------------------------------------------------
# Result Handling Utilities
#-----------------------------------------------------------------
    function init_global_storage(materials::Union{Material,Vector{Material}}, dim::Int)
        mat = materials isa Vector ? materials[1] : materials
        
        # mat_type = resolve_material_type(mat.type)
        if mat.type == [:thermal]
            return (K = zeros(dim, dim),)

        elseif  mat.type == [:elastic] && !haskey(mat.properties, :α_p)
            nstr = get_load_case_count(mat, dim)
            return (C = zeros(nstr, nstr),)

        elseif  mat.type == [:elastic, :electric]
            nstr = get_load_case_count( Material([:elastic], mat.dim, mat.symmetry, mat.properties, mat.tensors, mat.B_types,mat.mass_properties), dim)
            if mat.symmetry == :out_of_plane
                 nelec = 3  # Assuming 3D electric field components
            else
                nelec = dim
            end
           
            return (
                C = zeros(nstr, nstr),
                e = zeros(nelec, nstr),
                ϵ = zeros(nelec, nelec)
            )

        elseif  mat.type == :thermoelastic
            return (
                C = zeros(nstr, nstr),
                β = zeros(nstr, dim),
                κ = zeros(dim, dim)
            )

         elseif mat.type == [:elastic, :pressure] || haskey(mat.properties, :α_p)
            nstr = get_load_case_count( Material([:elastic], mat.dim, mat.symmetry, mat.properties, mat.tensors, mat.B_types,mat.mass_properties), dim)
            return (
                C = zeros(nstr, nstr),
                Β = zeros(nstr)
            )

        else
            error("Unsupported material type: $( mat.type)")
        end
    end

    function accumulate_results!(acc, elem_result)
        for field in propertynames(elem_result)
            if field != :volume && hasproperty(acc, field)
                acc_field = getproperty(acc, field)
                elem_field = getproperty(elem_result, field)
                acc_field .+= elem_field
            end
        end
    end

    function finalize_results(results, total_volume)
        final = deepcopy(results[1])
        for field in propertynames(final)
            final_field = getproperty(final, field)
            sum_field = sum(r -> getproperty(r, field), results)
            final_field .= sum_field ./ total_volume
        end
        return final
    end

#-----------------------------------------------------------------
# Material Property Calculations
#-----------------------------------------------------------------
    function compute_element_contribution(::Val{:thermal}, mat::Material, B_dict, elem_id, jacobian_data, gauss_data, Uᵉ)
        κ = mat.tensors[1]
        B = B_dict[:temperature_gradient]
        K_eff = zero(κ)
        volume = 0.0

        for qp in eachindex(gauss_data.weights)
            detJ, _ = jacobian_data[elem_id][qp]
            scale =  abs(detJ) * gauss_data.weights[qp]
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
            scale = abs(detJ) * gauss_data.weights[qp]
        
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
        C, ϵ, e = mat.tensors
        B_u = B_dict[:strain]
        B_e = B_dict[:electric_field]
        Uᵤ = Uᵉ.U_u
        Uᵩ = Uᵉ.U_ϕ


        C_eff = zero(C)
        e_eff = zero(e)
        ϵ_eff = zero(ϵ)
        volume = 0.0

        for qp in eachindex(gauss_data.weights)
            detJ, _ = jacobian_data[elem_id][qp]
            scale =  abs(detJ) * gauss_data.weights[qp]
            B_u_qp = B_u[qp]
            B_ϕ_qp = B_e[qp]
            
            if mat.symmetry == :out_of_plane
                M = zero(C); M[4, 4] = 1.0
                ε = B_u_qp * Uᵤ + M
            else
                ε = B_u_qp * Uᵤ
            end
          
            E = -B_ϕ_qp * Uᵩ
            # Compute stress and electric displacement
            σ = C * ε - transpose(e) * E
            
            D = e * ε + ϵ * E

            E_ψ = B_ϕ_qp * Uᵉ.V_ϕ
            if mat.symmetry == :out_of_plane 
                D_k = zero(ϵ)
                D_k[3,3] = 1.0  
               
                E_total = D_k + E_ψ
            else
                E_total = E_ψ  
            end
            # Compute strain and electric field
            ε_ψ = B_u_qp * Uᵉ.V_u
          
            
            # Compute electric displacement

            D_ψ = e * ε_ψ + ϵ * E_total

            # Accumulate effective properties
            C_eff += σ  * scale
            e_eff += D  * scale
            ϵ_eff += D_ψ * scale

            volume += scale
        end

        return (C = C_eff, e = e_eff, ϵ = ϵ_eff, volume = volume)
    end

    function compute_element_contribution(::Val{:poroelastic}, mat::Material, B_dict, elem_id, jacobian_data, gauss_data, Uᵉ)
        C = mat.tensors[1]
        B = B_dict[:strain]
        Uᵤ = Uᵉ[:,1:end-1]
        Uₚ = Uᵉ[:,end]
     
        C_eff = zero(C)
        Β_eff = zeros(size(C,1))
   
        volume = 0.0

        for qp in eachindex(gauss_data.weights)
            detJ, _ = jacobian_data[elem_id][qp]
            scale = abs(detJ) * gauss_data.weights[qp]
        
            if mat.symmetry == :out_of_plane
                M = zero(C); M[4, 4] = 1.0
                ε = B[qp] * Uᵤ + M
            else
                ε = B[qp] * Uᵤ
            end
         
            
            Β_eff -= C * B[qp] * Uₚ  * scale
            σ = C * ε
            C_eff += σ  * scale
            volume += scale 
        end

                
        return (C = C_eff,Β = Β_eff, volume = volume)
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
        # mat_type = resolve_material_type(mat.type)
        if mat.type in ([:thermal], :electric, :pressure)
            if mat.symmetry == :isotropic
                return dim  
            elseif mat.symmetry == :out_of_plane && dim == 2 # Only in 2D
                return 3
            else
                error("Unsupported symmetry: $(mat.symmetry)")
            end

        elseif mat.type == [:elastic] 
            if mat.symmetry == :isotropic
                return div(dim * (dim + 1), 2)  # 2D:3, 3D:6
            elseif mat.symmetry == :transversely_isotropic
                return dim == 2 ? 4 : 6
            elseif mat.symmetry == :orthotropic
                return dim == 2 ? 3 : 6
            elseif mat.symmetry == :out_of_plane && dim == 2 # Only in 2D
                return 4
            else
                error("Unsupported symmetry: $(mat.symmetry)")
            end
        else
            error("Unsupported material type: $(mat.type)")
        end
    end

    function ntens(dim)
        dim == 2 ? 3 : 6
    end