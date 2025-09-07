include("compute_flux.jl") 

"""
ElementMatrices for thermal property calculations.
"""
mutable struct ThermalPropertyMatrices{T<:AbstractFloat}
    elem_sol::Matrix{T}         # Element solution vector 
    K_eff::Matrix{T}            # Effective conductivity contribution 
    volume::Base.RefValue{T}    # Element volume
    ∇T_qp::Matrix{T}            # temperature gradient at quadrature points
    q_qp::Matrix{T}             # heat flux at quadrature points
end

"""
ElementMatrices for elastic property calculations.
"""
mutable struct ElasticPropertyMatrices{T<:AbstractFloat}
    elem_sol::Matrix{T}         # Element solution vector 
    C_eff::Matrix{T}            # Effective stiffness contribution 
    volume::Base.RefValue{T}    # Element volume
    ε_qp::Matrix{T}             # strain at quadrature points
    σ_qp::Matrix{T}             # stress at quadrature points
end

"""
ElementMatrices for piezoelectric property calculations.
"""
mutable struct PiezoelectricPropertyMatrices{T<:AbstractFloat}
    C_eff::Matrix{T}
    e_eff::Matrix{T}
    ε_eff::Matrix{T}
    volume::Base.RefValue{T}

    # SOLUTION MATRICES 
    Uu_sol::Matrix{T}
    Uϕ_sol::Matrix{T}
    Vu_sol::Matrix{T}
    Vϕ_sol::Matrix{T}

    # QUADRATURE POINT Storage
    ε_qp_u::Matrix{T}
    ε_qp_v::Matrix{T} 
    E_qp_u::Matrix{T}
    E_qp_v::Matrix{T}
    σ_qp::Matrix{T}
    D_qp_u::Matrix{T}
    D_qp_v::Matrix{T}
end

"""
ElementMatrices for poroelastic property calculations.
"""
mutable struct PoroelasticPropertyMatrices{T<:AbstractFloat}
    # ACCUMULATORS 
    C_eff::Matrix{T}
    Β_eff::Vector{T} # Stores ∫σ_p dV, which becomes Β_eff
    volume::Base.RefValue{T}

    # SOLUTION VECTORS 
    Uu_sol::Matrix{T} # Solution matrix for mechanical load cases
    Up_sol::Vector{T} # Solution vector for the single pressure load case

    # QUADRATURE POINT Storage
    ε_qp_u::Matrix{T} # Strain matrix for mechanical cases
    σ_qp_u::Matrix{T} # Stress matrix for mechanical cases
    ε_qp_p::Vector{T} # Strain vector for pressure case
    σ_qp_p::Vector{T} # Stress vector for pressure case
end

# Create Property Matrices 

function create_property_matrices(mat::ThermalMaterial, NN::Int, n_load_cases::Int, T::Type=Float64)
    D = mat.dim
    return ThermalPropertyMatrices(
        zeros(T, NN, n_load_cases),
        zeros(T, D, n_load_cases),
        Ref(zero(T)),
        zeros(T, D, n_load_cases),
        zeros(T, D, n_load_cases)
    )
end

function create_property_matrices(mat::ElasticMaterial, NN::Int, n_load_cases::Int, T::Type=Float64)
    D, nstr = mat.dim, rows_voigt(mat.dim, mat.symmetry)
    return ElasticPropertyMatrices(
        zeros(T, D * NN, n_load_cases),
        zeros(T, nstr, n_load_cases),
        Ref(zero(T)),
        zeros(T, nstr, n_load_cases),
        zeros(T, nstr, n_load_cases)
    )
end

function create_property_matrices(mat::PiezoelectricMaterial, NN::Int, n_mech_lc::Int, n_elec_lc::Int, T::Type=Float64)
    D = mat.dim
    nstr = rows_voigt(D, mat.symmetry)
    rG = vec_rows(D, mat.symmetry)
   
    return PiezoelectricPropertyMatrices{T}(
        zeros(T, nstr, n_mech_lc),   # C_eff
        zeros(T, rG, n_mech_lc),     # e_eff
        zeros(T, rG, n_elec_lc),     # ε_eff
        Ref(zero(T)),                # volume

        zeros(T, D * NN, n_mech_lc), # Uu_sol
        zeros(T, NN, n_mech_lc),     # Uϕ_sol
        zeros(T, D * NN, n_elec_lc), # Vu_sol
        zeros(T, NN, n_elec_lc),     # Vϕ_sol

        zeros(T, nstr, n_mech_lc), # ε_qp_u 
        zeros(T, nstr, n_elec_lc), # ε_qp_v 
        zeros(T, rG, n_mech_lc),   # E_qp_u 
        zeros(T, rG, n_elec_lc),   # E_qp_v 
        zeros(T, nstr, n_mech_lc), # σ_qp 
        zeros(T, rG, n_mech_lc),   # D_qp_u 
        zeros(T, rG, n_elec_lc)    # D_qp_v 
    )
end

function create_property_matrices(mat::PoroelasticMaterial, NN::Int, n_mech_load_cases::Int, T::Type=Float64)
    D, nstr = mat.dim, rows_voigt(mat.dim, mat.symmetry)
    n_dofs_elem = D * NN
    
    return PoroelasticPropertyMatrices(
        zeros(T, nstr, n_mech_load_cases), # C_eff
        zeros(T, nstr),                   # Β_eff
        Ref(zero(T)),                     # volume
        zeros(T, n_dofs_elem, n_mech_load_cases), # Uu_sol
        zeros(T, n_dofs_elem),            # Up_sol
        zeros(T, nstr, n_mech_load_cases),# ε_qp_u 
        zeros(T, nstr, n_mech_load_cases),# σ_qp 
        zeros(T, nstr),                   # ε_qp_p 
        zeros(T, nstr)                    # σ_qp_p 
    )
end

# Element Contribution 

function compute_element_contribution!(
    prop_matrices::ThermalPropertyMatrices, mat::ThermalMaterial, B_e::BElem,
    shp::FESD{NN,NGP,D}, jacs::BJacs, e::Int) where {NN,NGP,D}

    fill!(prop_matrices.K_eff, 0.0)
    prop_matrices.volume[] = 0.0
    κ = mat.κ
    G_e = B_e.vector_gradient_operator

    @inbounds for qp in 1:NGP
        dV = abs(jacs[qp, e][1]) * shp.weights[qp]
        G_qp = @view G_e[:, :, qp]

        mul!(prop_matrices.∇T_qp, G_qp, prop_matrices.elem_sol) # ∇T = G * Tₑ
        mul!(prop_matrices.q_qp, κ, prop_matrices.∇T_qp)        # q = κ * ∇T

        # Accumulate effective conductivity: K_eff += q * dV
        prop_matrices.K_eff .+= prop_matrices.q_qp .* dV
        prop_matrices.volume[] += dV
    end
end

function compute_element_contribution!(
    prop_matrices::ElasticPropertyMatrices, mat::ElasticMaterial, B_e::BElem,
    shp::FESD{NN,NGP,D}, jacs::BJacs, e::Int) where {NN,NGP,D}

    fill!(prop_matrices.C_eff, 0.0)
    prop_matrices.volume[] = 0.0
    C = mat.C
    B_voigt_e = B_e.voigt_gradient_operator

    @inbounds for qp in 1:NGP
        dV = abs(jacs[qp, e][1]) * shp.weights[qp]
        B_qp = @view B_voigt_e[:, :, qp]

        # Calculate microscopic strain: ε = B * uₑ
        mul!(prop_matrices.ε_qp, B_qp, prop_matrices.elem_sol)

        if mat.symmetry == :out_of_plane
            # adds the unit strain component in the out-of-plane direction (ε₃₃ = 1).
            prop_matrices.ε_qp[end] += 1.0
        end

        # Calculate microscopic stress: σ = C * ε
        if haskey(mat.properties, :cracks)
         
            E = mat.properties[:E]
            ν = mat.properties[:ν]
            λ = (E * ν) / ((1 + ν) * (1 - 2ν))
            μ = E / (2(1 + ν))
            
            if compare_energy(prop_matrices.ε_qp, λ, μ; tol=1e-12)
                # Reduce stiffness for cracked elements
                C_reduced = 1e-8 * C
                mul!(prop_matrices.σ_qp, C_reduced, prop_matrices.ε_qp)
            else
                mul!(prop_matrices.σ_qp, C, prop_matrices.ε_qp)
            end
        else
            # Compute stress: σ = C * ε
            mul!(prop_matrices.σ_qp, C, prop_matrices.ε_qp)
        end

        # Accumulate effective stiffness: C_eff += σ * dV
        prop_matrices.C_eff .+= prop_matrices.σ_qp .* dV
        prop_matrices.volume[] += dV
    end
end

function compute_element_contribution!(
    prop_matrices::PiezoelectricPropertyMatrices,
    mat::PiezoelectricMaterial,
    B_e::BElem,
    shp::FESD{NN,NGP,D},
    jacs::BJacs,
    e::Int) where {NN, NGP, D}

    fill!(prop_matrices.C_eff, 0.0)
    fill!(prop_matrices.e_eff, 0.0)
    fill!(prop_matrices.ε_eff, 0.0)
    prop_matrices.volume[] = 0.0

    C, ϵ_tensor, e_tensor = mat.C, mat.ε, mat.e
    e_tensor_T = transpose(e_tensor)
    Bv_e = B_e.voigt_gradient_operator
    Gs_e = B_e.vector_gradient_operator
 
    @inbounds for qp in 1:NGP

        dV = abs(jacs[qp, e][1]) * shp.weights[qp]
        Bv = @view Bv_e[:, :, qp]
        Gs = @view Gs_e[:, :, qp]

        # Contributions from Mechanical Load Cases (Uu, Uϕ)

        # Microscopic strain: ε = Bᵤ * Uᵤ
        mul!(prop_matrices.ε_qp_u, Bv, prop_matrices.Uu_sol)
        if mat.symmetry == :out_of_plane
            # Add unit matrix M where M[4,4] = 1.0 
            prop_matrices.ε_qp_u[4, 4] += 1.0
        end

        # Microscopic electric field: E = -Bᵩ * Uᵩ
        mul!(prop_matrices.E_qp_u, Gs, prop_matrices.Uϕ_sol, -1.0, 0.0)

        # Stress: σ = C*ε - eᵀ*E
        mul!(prop_matrices.σ_qp, C, prop_matrices.ε_qp_u)
        mul!(prop_matrices.σ_qp, e_tensor_T, prop_matrices.E_qp_u, -1.0, 1.0)

        # Electric Displacement: D = e*ε + εˢ*E
        mul!(prop_matrices.D_qp_u, e_tensor, prop_matrices.ε_qp_u)
        mul!(prop_matrices.D_qp_u, ϵ_tensor, prop_matrices.E_qp_u, 1.0, 1.0)

        # Contributions from Electrical Load Cases (Vu, Vϕ)

        # Strain from electrical loading: ε_ψ = Bᵤ * Vᵤ
        mul!(prop_matrices.ε_qp_v, Bv, prop_matrices.Vu_sol)

        # Electric field from electrical loading: E_ψ = Bᵩ * Vᵩ
        mul!(prop_matrices.E_qp_v, Gs, prop_matrices.Vϕ_sol)
        if mat.symmetry == :out_of_plane
            # Add unit matrix D_k where D_k[3,3] = 1.0 
            prop_matrices.E_qp_v[3, 3] += 1.0
        end

        # Electric Displacement: D_ψ = e*ε_ψ + εˢ*E_total
        mul!(prop_matrices.D_qp_v, e_tensor, prop_matrices.ε_qp_v)
        mul!(prop_matrices.D_qp_v, ϵ_tensor, prop_matrices.E_qp_v, 1.0, 1.0)

        # Accumulate Integrated Results

        # C_eff += σ * dV
        prop_matrices.C_eff .+= prop_matrices.σ_qp .* dV

        # e_eff += D * dV
        prop_matrices.e_eff .+= prop_matrices.D_qp_u .* dV

        # ε_eff += D_ψ * dV
        prop_matrices.ε_eff .+= prop_matrices.D_qp_v .* dV

        # Accumulate total volume
        prop_matrices.volume[] += dV
    end
end

function compute_element_contribution!(
    prop_matrices::PoroelasticPropertyMatrices,
    mat::PoroelasticMaterial,
    B_e::BElem,
    shp::FESD{NN,NGP,D},
    jacs::BJacs,
    e::Int) where {NN, NGP, D}

    fill!(prop_matrices.C_eff, 0.0)
    fill!(prop_matrices.Β_eff, 0.0)
    prop_matrices.volume[] = 0.0

    C = mat.C
    B_voigt_e = B_e.voigt_gradient_operator

    @inbounds for qp in 1:NGP
        dV = abs(jacs[qp, e][1]) * shp.weights[qp]
        B_qp = @view B_voigt_e[:, :, qp]

        # Mechanical Loading Cases 
        # Strain: ε = B * Uᵤ
        mul!(prop_matrices.ε_qp_u, B_qp, prop_matrices.Uu_sol)
        if mat.symmetry == :out_of_plane
            # Add unit out-of-plane strain to all mechanical load cases
            prop_matrices.ε_qp_u[end] .+= 1.0
        end
        
        # Stress: σ = C * ε
        mul!(prop_matrices.σ_qp_u, C, prop_matrices.ε_qp_u)

        # Accumulate: C_eff += σ * dV
        prop_matrices.C_eff .+= prop_matrices.σ_qp_u .* dV

        # Pressure Loading Case 
        # Strain from pressure field: ε_p = B * Uₚ
        mul!(prop_matrices.ε_qp_p, B_qp, prop_matrices.Up_sol)
        
        # Accumulate: Β_eff -= C * ε_p * dV
        mul!(prop_matrices.σ_qp_p, C, prop_matrices.ε_qp_p)
        prop_matrices.Β_eff .-= prop_matrices.σ_qp_p .* dV
        
        # Accumulate Total Volume
        prop_matrices.volume[] += dV
    end
end

# Solution Retrieval and Storage Initialization

function get_element_solution!(prop_matrices::Union{ThermalPropertyMatrices, ElasticPropertyMatrices}, 
                               solver_results, elem_dofs::AbstractVector{Int})
    for i in 1:size(prop_matrices.elem_sol, 2)
        prop_matrices.elem_sol[:, i] .= solver_results.U[i][elem_dofs]
    end
end

function get_element_solution!(prop_matrices::PiezoelectricPropertyMatrices, 
                               solver_results, elem_dofs_u, elem_dofs_ϕ)
    n_mech_lc = size(prop_matrices.Uu_sol, 2)
    n_elec_lc = size(prop_matrices.Vu_sol, 2)
 
    for i in 1:n_mech_lc
        prop_matrices.Uu_sol[:, i] .= solver_results.U_mech[i][elem_dofs_u]
        prop_matrices.Uϕ_sol[:, i] .= solver_results.V_elec[i][elem_dofs_ϕ]
    end
    for i in 1:n_elec_lc
        prop_matrices.Vu_sol[:, i] .= solver_results.U_mech[n_mech_lc + i][elem_dofs_u]
        prop_matrices.Vϕ_sol[:, i] .= solver_results.V_elec[n_mech_lc + i][elem_dofs_ϕ]
    end
end

function get_element_solution!(prop_matrices::PoroelasticPropertyMatrices, 
    solver_results::NamedTuple, elem_dofs::AbstractVector{Int})

    n_mech_load_cases = length(solver_results.U) - 1
    
    for i in 1:n_mech_load_cases
        prop_matrices.Uu_sol[:, i] .= solver_results.U[i][elem_dofs]
    end

    prop_matrices.Up_sol[:] .= solver_results.U[end][elem_dofs]
end

# initialization functions

function init_thermal_storage(mat::ThermalMaterial, n_lc::Int)
    return (K = zeros(mat.dim, n_lc),)
end

function init_elastic_storage(mat::ElasticMaterial, n_lc::Int)
    return (C = zeros(rows_voigt(mat.dim, mat.symmetry), n_lc),)
end

function init_piezo_storage(mat::PiezoelectricMaterial, n_mech_lc::Int, n_elec_lc::Int)
    D = mat.dim
    return (
        C = zeros(rows_voigt(D, mat.symmetry), n_mech_lc),
        e = zeros(vec_rows(D, mat.symmetry), n_mech_lc),
        ϵ = zeros(vec_rows(D, mat.symmetry), n_elec_lc)
    )
end

function init_poroelastic_storage(mat::PoroelasticMaterial, n_mech_load_cases::Int)
    nstr = rows_voigt(mat.dim, mat.symmetry)
    return (
        C = zeros(nstr, n_mech_load_cases),
        Β = zeros(nstr) 
    )
end

# specializations for each material 

function compute_effective_property_thermal(
    materials::Vector{ThermalMaterial{M}},
    connectivity::Matrix{Int},
    material_indices_per_element::Vector{Int},
    solver_results::NamedTuple,
    NN::Int,
    geom::GeometricData) where {M}

    Ne = size(connectivity, 1)
    nthreads = Threads.nthreads()
    n_load_cases = length(solver_results.U)
    
    result_template = init_thermal_storage(materials[1], n_load_cases)
    result_acc = [deepcopy(result_template) for _ in 1:nthreads]
    vol_acc = zeros(nthreads)
    
    prop_matrices_per_thread = [[create_property_matrices(m, NN, n_load_cases) for m in materials] for _ in 1:nthreads]
    
    dofs_per_node = get_dofs_per_node(materials[1])
    single_field_dofs_buf = [Vector{Int}(undef, NN * dofs_per_node) for _ in 1:nthreads]


    @inbounds Threads.@threads for e in 1:Ne
        tid = Threads.threadid()
        mat_idx = material_indices_per_element[e]
        mat = materials[mat_idx]
        prop_matrices = prop_matrices_per_thread[tid][mat_idx]
        conn_e = @view connectivity[e, :]
        dofs_buf = single_field_dofs_buf[tid]
        
        dof_offsets = 1:dofs_per_node
        for a in 1:NN
            base = (conn_e[a] - 1) * dofs_per_node
            cols = (a - 1) * dofs_per_node + 1 : a * dofs_per_node
            dofs_buf[cols] .= base .+ dof_offsets
        end
        
        get_element_solution!(prop_matrices, solver_results, dofs_buf)
        
        compute_element_contribution!(
            prop_matrices, mat, geom.differential_operator_matrices[e],
            geom.shape_function_data, geom.jacobian_transformation_data, e
        )
        
        result_acc[tid].K .+= prop_matrices.K_eff
        vol_acc[tid] += prop_matrices.volume[]
    end
    
    # Finalize results
    total_volume = sum(vol_acc)
    final_K = sum(r.K for r in result_acc) ./ total_volume
    
    return (K = final_K,), total_volume
end

function compute_effective_property_elastic(
    materials::Vector{ElasticMaterial{M}},
    connectivity::Matrix{Int},
    material_indices_per_element::Vector{Int},
    solver_results::NamedTuple,
    NN::Int,
    geom::GeometricData) where {M}

    Ne = size(connectivity, 1)
    nthreads = Threads.nthreads()
    n_load_cases = length(solver_results.U)
    
    result_template = init_elastic_storage(materials[1], n_load_cases)
    result_acc = [deepcopy(result_template) for _ in 1:nthreads]
    vol_acc = zeros(nthreads)
    
    prop_matrices_per_thread = [[create_property_matrices(m, NN, n_load_cases) for m in materials] for _ in 1:nthreads]
    
    dofs_per_node = get_dofs_per_node(materials[1])
    single_field_dofs_buf = [Vector{Int}(undef, NN * dofs_per_node) for _ in 1:nthreads]
    
    @inbounds Threads.@threads for e in 1:Ne
        tid = Threads.threadid()
        mat_idx = material_indices_per_element[e]
        mat = materials[mat_idx]
        prop_matrices = prop_matrices_per_thread[tid][mat_idx]
        conn_e = @view connectivity[e, :]
        dofs_buf = single_field_dofs_buf[tid]
        
        dof_offsets = 1:dofs_per_node
        for a in 1:NN
            base = (conn_e[a] - 1) * dofs_per_node
            cols = (a - 1) * dofs_per_node + 1 : a * dofs_per_node
            dofs_buf[cols] .= base .+ dof_offsets
        end
        
        get_element_solution!(prop_matrices, solver_results, dofs_buf)
        
        compute_element_contribution!(
            prop_matrices, mat, geom.differential_operator_matrices[e],
            geom.shape_function_data, geom.jacobian_transformation_data, e
        )
        
        result_acc[tid].C .+= prop_matrices.C_eff
        vol_acc[tid] += prop_matrices.volume[]
    end
    
    total_volume = sum(vol_acc)
    final_C = sum(r.C for r in result_acc) ./ total_volume
    
    return (C = final_C,), total_volume
end

function compute_effective_property_piezo(
    materials::Vector{PiezoelectricMaterial{M1,M2,M3}},
    connectivity::Matrix{Int},
    material_indices_per_element::Vector{Int},
    solver_results::NamedTuple,
    NN::Int,
    dim::Int,
    geom::GeometricData) where {M1,M2,M3}
  
    Ne = size(connectivity, 1)
    nthreads = Threads.nthreads()
    n_mech_lc = size(materials[1].C,1)
    n_elec_lc = size(materials[1].ε,1) 
    
    result_template = init_piezo_storage(materials[1], n_mech_lc, n_elec_lc)
    result_acc = [deepcopy(result_template) for _ in 1:nthreads]
    vol_acc = zeros(nthreads)
    
    prop_matrices_per_thread = [[create_property_matrices(m, NN, n_mech_lc, n_elec_lc) for m in materials] for _ in 1:nthreads]
    
    dofs_per_node = get_dofs_per_node(materials[1])
    Nu = dofs_u_per_node(materials[1])
    mech_dofs_buf = [Vector{Int}(undef, Nu * NN) for _ in 1:nthreads]
    scalar_dofs_buf = [Vector{Int}(undef, NN) for _ in 1:nthreads]

    @inbounds Threads.@threads for e in 1:Ne
        tid = Threads.threadid()
        mat_idx = material_indices_per_element[e]
        mat = materials[mat_idx]
        prop_matrices = prop_matrices_per_thread[tid][mat_idx]
        conn_e = @view connectivity[e, :]

        local_mech_dofs!(mech_dofs_buf[tid], conn_e, dofs_per_node, dim)
        local_scal_dofs!(scalar_dofs_buf[tid], conn_e, dofs_per_node, dim, 1)
        
        get_element_solution!(prop_matrices, solver_results, mech_dofs_buf[tid], scalar_dofs_buf[tid])
            
        compute_element_contribution!(
            prop_matrices, mat, geom.differential_operator_matrices[e],
            geom.shape_function_data, geom.jacobian_transformation_data, e
        )
        
        result_acc[tid].C .+= prop_matrices.C_eff
        result_acc[tid].e .+= prop_matrices.e_eff
        result_acc[tid].ϵ .+= prop_matrices.ε_eff
        vol_acc[tid] += prop_matrices.volume[]
    end
    
    # Finalize results
    total_volume = sum(vol_acc)
    final_C = sum(r.C for r in result_acc) ./ total_volume
    final_e = sum(r.e for r in result_acc) ./ total_volume
    final_ε = sum(r.ϵ for r in result_acc) ./ total_volume
    
    return (C = final_C, e = final_e, ϵ = final_ε), total_volume
end

function compute_effective_property_poroelastic(
    materials::Vector{PoroelasticMaterial{M}},
    connectivity::Matrix{Int},
    material_indices_per_element::Vector{Int},
    solver_results::NamedTuple,
    NN::Int,
    geom::GeometricData) where {M}

    Ne = size(connectivity, 1)
    nthreads = Threads.nthreads()
    n_mech_load_cases = length(solver_results.U) - 1
    
    result_template = init_poroelastic_storage(materials[1], n_mech_load_cases)
    result_acc = [deepcopy(result_template) for _ in 1:nthreads]
    vol_acc = zeros(nthreads)
    
    prop_matrices_per_thread = [[create_property_matrices(m, NN, n_mech_load_cases) for m in materials] for _ in 1:nthreads]
    
    dofs_per_node = get_dofs_per_node(materials[1])
    dofs_buf_per_thread = [Vector{Int}(undef, NN * dofs_per_node) for _ in 1:nthreads]
    
    @inbounds Threads.@threads for e in 1:Ne
        tid = Threads.threadid()
        mat_idx = material_indices_per_element[e]
        mat = materials[mat_idx]
        prop_matrices = prop_matrices_per_thread[tid][mat_idx]
        conn_e = @view connectivity[e, :]
        dofs_buf = dofs_buf_per_thread[tid]

        for a in 1:NN
            base = (conn_e[a] - 1) * dofs_per_node
            cols = (a - 1) * dofs_per_node + 1 : a * dofs_per_node
            dofs_buf[cols] .= base .+ (1:dofs_per_node)
        end
        
        get_element_solution!(prop_matrices, solver_results, dofs_buf)

        compute_element_contribution!(
            prop_matrices, mat, geom.differential_operator_matrices[e],
            geom.shape_function_data, geom.jacobian_transformation_data, e
        )

        result_acc[tid].C .+= prop_matrices.C_eff
        result_acc[tid].Β .+= prop_matrices.Β_eff
        vol_acc[tid] += prop_matrices.volume[]
    end

    total_volume = sum(vol_acc)
    final_C = sum(r.C for r in result_acc) ./ total_volume
    final_B = sum(r.Β for r in result_acc) ./ total_volume
    
    return (C = final_C, Β = final_B), total_volume
end

# Main dispatcher function 

function compute_effective_property(
    materials::Vector{M},
    connectivity::Matrix{Int},
    connect_elem_phase::Vector{Int},
    solver_results::NamedTuple,
    dim::Int,
    geom::GeometricData) where {M<:AbstractMaterial}

    Ne, NN = size(connectivity)
    
    # Phase handling 
    unique_phases = sort(unique(connect_elem_phase))
    phase_to_material_idx = Dict{Int, Int}(phase => i for (i, phase) in enumerate(unique_phases))
    material_indices_per_element = [phase_to_material_idx[p] for p in connect_elem_phase]

    # dispatch based on material type
    if M <: ThermalMaterial
        return compute_effective_property_thermal(
            materials, connectivity, material_indices_per_element, solver_results, NN, geom
        )
    elseif M <: ElasticMaterial
        return compute_effective_property_elastic(
            materials, connectivity, material_indices_per_element, solver_results, NN, geom
        )
    elseif M <: PiezoelectricMaterial
        return compute_effective_property_piezo(
            materials, connectivity, material_indices_per_element, solver_results, NN, dim, geom
        )
    elseif M <: PoroelasticMaterial
        return compute_effective_property_poroelastic(
            materials, connectivity, material_indices_per_element, solver_results, NN, geom
        )
        
    else
        error("Unsupported material type: $(typeof(M))")
    end
end

# single material case
function compute_effective_property(
    material::M,
    connectivity::Matrix{Int},
    solver_results::NamedTuple,
    dim::Int,
    geom::GeometricData) where {M<:AbstractMaterial}
    
    Ne = size(connectivity, 1)
    connect_elem_phase = ones(Int, Ne)
    return compute_effective_property([material], connectivity, connect_elem_phase, solver_results, dim, geom)
end