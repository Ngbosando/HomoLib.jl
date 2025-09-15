#-----------------------------------------------------------------
# Field Recovery Function
#-----------------------------------------------------------------

"""
Storage for thermal field recovery results.
"""
struct ThermalFieldStorage{T<:AbstractFloat}
    flux::Matrix{T}         # [num_nodes × dim] - heat flux
    grad_temp::Matrix{T}    # [num_nodes × dim] - temperature gradient
end

"""
Storage for elastic field recovery results.
"""
struct ElasticFieldStorage{T<:AbstractFloat}
    stress::Matrix{T}       # [num_nodes × nstr] - stress components
    strain::Matrix{T}       # [num_nodes × nstr] - strain components
end

"""
Storage for piezoelectric field recovery results.
"""
struct PiezoelectricFieldStorage{T<:AbstractFloat}
    stress::Matrix{T}       # [num_nodes × nstr] - stress components
    strain::Matrix{T}       # [num_nodes × nstr] - strain components
    elec_disp::Matrix{T}    # [num_nodes × nelec] - electric displacement
    elec_field::Matrix{T}   # [num_nodes × nelec] - electric field
end

"""
Storage for poroelastic field recovery results.
"""
struct PoroelasticFieldStorage{T<:AbstractFloat}
    stress::Matrix{T}       # [num_nodes × nstr] - stress components
    strain::Matrix{T}       # [num_nodes × nstr] - strain components
    pressure::Vector{T}     # [num_nodes] - pore pressure
end


# Matrices for Field Recovery


"""
Working matrices for thermal field recovery at element level.
"""
mutable struct ThermalRecoveryMatrices{T<:AbstractFloat}
    elem_sol::Vector{T}     # Element solution vector
    flux_sum::Vector{T}     # Accumulated flux
    grad_sum::Vector{T}     # Accumulated temperature gradient
    volume_sum::Base.RefValue{T}
    
    # Quadrature point buffers
    ∇T_qp::Vector{T}
    q_qp::Vector{T}
end

"""
Working matrices for elastic field recovery at element level.
"""
mutable struct ElasticRecoveryMatrices{T<:AbstractFloat}
    elem_sol::Vector{T}     # Element solution vector
    stress_sum::Vector{T}   # Accumulated stress
    strain_sum::Vector{T}   # Accumulated strain
    volume_sum::Base.RefValue{T}
    
    # Quadrature point buffers
    ε_qp::Vector{T}
    σ_qp::Vector{T}
end

"""
Working matrices for piezoelectric field recovery at element level.
"""
mutable struct PiezoelectricRecoveryMatrices{T<:AbstractFloat}
    elem_sol_u::Vector{T}   # Mechanical solution
    elem_sol_ϕ::Vector{T}   # Electric potential solution
    stress_sum::Vector{T}   # Accumulated stress
    strain_sum::Vector{T}   # Accumulated strain
    elec_disp_sum::Vector{T} # Accumulated electric displacement
    elec_field_sum::Vector{T} # Accumulated electric field
    volume_sum::Base.RefValue{T}
    
    # Quadrature point buffers
    ε_qp::Vector{T}
    E_qp::Vector{T}
    σ_qp::Vector{T}
    D_qp::Vector{T}
end

"""
Working matrices for poroelastic field recovery at element level.
"""
mutable struct PoroelasticRecoveryMatrices{T<:AbstractFloat}
    elem_sol::Vector{T}     # Element solution vector
    stress_sum::Vector{T}   # Accumulated stress
    strain_sum::Vector{T}   # Accumulated strain
    pressure_sum::Base.RefValue{T} # Accumulated pressure
    volume_sum::Base.RefValue{T}
    
    # Quadrature point buffers
    ε_qp::Vector{T}
    σ_qp::Vector{T}
end

# Functions for Storage Creation


function create_field_storage(mat::ThermalMaterial, num_nodes::Int, T::Type=Float64)
    D = mat.dim
    return ThermalFieldStorage{T}(
        zeros(T, num_nodes, D),     # flux
        zeros(T, num_nodes, D)      # grad_temp
    )
end

function create_field_storage(mat::ElasticMaterial, num_nodes::Int, T::Type=Float64)
    nstr = rows_voigt(mat.dim, mat.symmetry)
    return ElasticFieldStorage{T}(
        zeros(T, num_nodes, nstr),  # stress
        zeros(T, num_nodes, nstr)   # strain
    )
end

function create_field_storage(mat::PiezoelectricMaterial, num_nodes::Int, T::Type=Float64)
    nstr = rows_voigt(mat.dim, mat.symmetry)
    nelec = vec_rows(mat.dim, mat.symmetry)
    return PiezoelectricFieldStorage{T}(
        zeros(T, num_nodes, nstr),  # stress
        zeros(T, num_nodes, nstr),  # strain
        zeros(T, num_nodes, nelec), # elec_disp
        zeros(T, num_nodes, nelec)  # elec_field
    )
end

function create_field_storage(mat::PoroelasticMaterial, num_nodes::Int, T::Type=Float64)
    nstr = rows_voigt(mat.dim, mat.symmetry)
    return PoroelasticFieldStorage{T}(
        zeros(T, num_nodes, nstr),  # stress
        zeros(T, num_nodes, nstr),  # strain
        zeros(T, num_nodes)         # pressure
    )
end

# Working Matrix Factories

function create_recovery_matrices(mat::ThermalMaterial, NN::Int, T::Type=Float64)
    D = mat.dim
    return ThermalRecoveryMatrices{T}(
        zeros(T, NN),
        zeros(T, D),
        zeros(T, D),
        Ref(zero(T)),
        zeros(T, D),
        zeros(T, D)
    )
end

function create_recovery_matrices(mat::ElasticMaterial, NN::Int, T::Type=Float64)
    nstr = rows_voigt(mat.dim, mat.symmetry)
    return ElasticRecoveryMatrices{T}(
        zeros(T, mat.dim * NN),
        zeros(T, nstr),
        zeros(T, nstr),
        Ref(zero(T)),
        zeros(T, nstr),
        zeros(T, nstr)
    )
end

function create_recovery_matrices(mat::PiezoelectricMaterial, NN::Int, T::Type=Float64)
    nstr = rows_voigt(mat.dim, mat.symmetry)
    nelec = vec_rows(mat.dim, mat.symmetry)
    return PiezoelectricRecoveryMatrices{T}(
        zeros(T, mat.dim * NN),
        zeros(T, NN),
        zeros(T, nstr),
        zeros(T, nstr),
        zeros(T, nelec),
        zeros(T, nelec),
        Ref(zero(T)),
        zeros(T, nstr),
        zeros(T, nelec),
        zeros(T, nstr),
        zeros(T, nelec)
    )
end

function create_recovery_matrices(mat::PoroelasticMaterial, NN::Int, T::Type=Float64)
    nstr = rows_voigt(mat.dim, mat.symmetry)
    return PoroelasticRecoveryMatrices{T}(
        zeros(T, mat.dim * NN),
        zeros(T, nstr),
        zeros(T, nstr),
        Ref(zero(T)),
        Ref(zero(T)),
        zeros(T, nstr),
        zeros(T, nstr)
    )
end


# Element Field Recovery


function recover_element_fields!(
    recovery_matrices::ThermalRecoveryMatrices,
    mat::ThermalMaterial,
    B_e::BElem,
    shp::FESD{NN,NGP,D},
    jacs::BJacs,
    e::Int) where {NN,NGP,D}
    
    fill!(recovery_matrices.flux_sum, 0.0)
    fill!(recovery_matrices.grad_sum, 0.0)
    recovery_matrices.volume_sum[] = 0.0
    
    κ = mat.κ
    G_e = B_e.vector_gradient_operator
    
    @inbounds for qp in 1:NGP
        dV = abs(jacs[qp, e][1]) * shp.weights[qp]
        G_qp = @view G_e[:, :, qp]
        
        # Compute temperature gradient: ∇T = G * T_e
        mul!(recovery_matrices.∇T_qp, G_qp, recovery_matrices.elem_sol)
        
        # Compute heat flux: q = -κ * ∇T
        mul!(recovery_matrices.q_qp, κ, recovery_matrices.∇T_qp, -1.0, 0.0)
        
        # Accumulate weighted results
        recovery_matrices.flux_sum .+= recovery_matrices.q_qp .* dV
        recovery_matrices.grad_sum .+= recovery_matrices.∇T_qp .* dV
        recovery_matrices.volume_sum[] += dV
    end
end

function recover_element_fields!(
    recovery_matrices::ElasticRecoveryMatrices,
    mat::ElasticMaterial,
    B_e::BElem,
    shp::FESD{NN,NGP,D},
    jacs::BJacs,
    e::Int) where {NN,NGP,D}
    
    fill!(recovery_matrices.stress_sum, 0.0)
    fill!(recovery_matrices.strain_sum, 0.0)
    recovery_matrices.volume_sum[] = 0.0
    
    C = mat.C
    B_voigt_e = B_e.voigt_gradient_operator
    
    @inbounds for qp in 1:NGP
        dV = abs(jacs[qp, e][1]) * shp.weights[qp]
        B_qp = @view B_voigt_e[:, :, qp]
        
        # Compute strain: ε = B * u_e
        mul!(recovery_matrices.ε_qp, B_qp, recovery_matrices.elem_sol)

        if mat.symmetry == :out_of_plane
            recovery_matrices.ε_qp[end] += 1.0
        end
         
        # Check for cracked material 
        if haskey(mat.properties, :cracks)
            E = mat.properties[:E]
            ν = mat.properties[:ν]
            λ = (E * ν) / ((1 + ν) * (1 - 2ν))
            μ = E / (2(1 + ν))
            
            if !compare_energy(recovery_matrices.ε_qp, λ, μ; tol=1e-12)
                # Reduce stiffness for cracked elements
               
                C_reduced = 1e-8 * C
                mul!(recovery_matrices.σ_qp, C_reduced, recovery_matrices.ε_qp)
            else
                mul!(recovery_matrices.σ_qp, C, recovery_matrices.ε_qp)
            end
        else
            # Compute stress: σ = C * ε
            mul!(recovery_matrices.σ_qp, C, recovery_matrices.ε_qp)
        end
        
        # Accumulate weighted results
        recovery_matrices.stress_sum .+= recovery_matrices.σ_qp .* dV
        recovery_matrices.strain_sum .+= recovery_matrices.ε_qp .* dV
        recovery_matrices.volume_sum[] += dV
    end
end

function recover_element_fields!(
    recovery_matrices::PiezoelectricRecoveryMatrices,
    mat::PiezoelectricMaterial,
    B_e::BElem,
    shp::FESD{NN,NGP,D},
    jacs::BJacs,
    e::Int) where {NN,NGP,D}
    

    fill!(recovery_matrices.stress_sum, 0.0)
    fill!(recovery_matrices.strain_sum, 0.0)
    fill!(recovery_matrices.elec_disp_sum, 0.0)
    fill!(recovery_matrices.elec_field_sum, 0.0)
    recovery_matrices.volume_sum[] = 0.0
    
    C, ε_tensor, e_tensor = mat.C, mat.ε, mat.e
    B_voigt_e = B_e.voigt_gradient_operator
    G_e = B_e.vector_gradient_operator
    
    @inbounds for qp in 1:NGP
        dV = abs(jacs[qp, e][1]) * shp.weights[qp]
        B_qp = @view B_voigt_e[:, :, qp]
        G_qp = @view G_e[:, :, qp]
        
        # Compute strain: ε = B * u_e
        mul!(recovery_matrices.ε_qp, B_qp, recovery_matrices.elem_sol_u)
        
        # Compute electric field: E = G * φ_e
        mul!(recovery_matrices.E_qp, G_qp, recovery_matrices.elem_sol_ϕ, 1.0, 0.0)
        
        # Compute stress: σ = C*ε - e^T*E
        mul!(recovery_matrices.σ_qp, C, recovery_matrices.ε_qp)
        mul!(recovery_matrices.σ_qp, transpose(e_tensor), recovery_matrices.E_qp, -1.0, 1.0)
        
        # Compute electric displacement: D = e*ε + ε*E
        mul!(recovery_matrices.D_qp, e_tensor, recovery_matrices.ε_qp)
        mul!(recovery_matrices.D_qp, ε_tensor, recovery_matrices.E_qp, 1.0, 1.0)
        
        # Accumulate weighted results
        recovery_matrices.stress_sum .+= recovery_matrices.σ_qp .* dV
        recovery_matrices.strain_sum .+= recovery_matrices.ε_qp .* dV
        recovery_matrices.elec_disp_sum .+= recovery_matrices.D_qp .* dV
        recovery_matrices.elec_field_sum .+= recovery_matrices.E_qp .* dV
        recovery_matrices.volume_sum[] += dV
    end
end

function recover_element_fields!(
    recovery_matrices::PoroelasticRecoveryMatrices,
    mat::PoroelasticMaterial,
    B_e::BElem,
    shp::FESD{NN,NGP,D},
    jacs::BJacs,
    e::Int) where {NN,NGP,D}
    
    fill!(recovery_matrices.stress_sum, 0.0)
    fill!(recovery_matrices.strain_sum, 0.0)
    recovery_matrices.pressure_sum[] = 0.0
    recovery_matrices.volume_sum[] = 0.0
    
    C = mat.C
    B_voigt_e = B_e.voigt_gradient_operator
    
    @inbounds for qp in 1:NGP
        dV = abs(jacs[qp, e][1]) * shp.weights[qp]
        B_qp = @view B_voigt_e[:, :, qp]
        N = shp.N[qp]
        
        # Compute strain: ε = B * u_e
        mul!(recovery_matrices.ε_qp, B_qp, recovery_matrices.elem_sol)
        
        # Compute stress: σ = C * ε
        mul!(recovery_matrices.σ_qp, C, recovery_matrices.ε_qp)
        
        # Compute pressure at quadrature point 
        pressure_qp = 0.0
        for a in 1:NN
            node_idx = (a-1)*get_dofs_per_node(mat) + mat.dim + 1
            if node_idx <= length(recovery_matrices.elem_sol)
                pressure_qp += N[a] * recovery_matrices.elem_sol[node_idx] #dont work need to correct it 
            end
        end
        
        # Accumulate weighted results
        recovery_matrices.stress_sum .+= recovery_matrices.σ_qp .* dV
        recovery_matrices.strain_sum .+= recovery_matrices.ε_qp .* dV
        recovery_matrices.pressure_sum[] += pressure_qp * dV
        recovery_matrices.volume_sum[] += dV
    end
end

# Solution Extraction

function extract_element_solution!(
    recovery_matrices::Union{ThermalRecoveryMatrices, ElasticRecoveryMatrices, PoroelasticRecoveryMatrices},
    Uresult::NamedTuple,
    elem_dofs::AbstractVector{Int})
    
    if hasfield(typeof(Uresult), :T)
        recovery_matrices.elem_sol .= Uresult.T[elem_dofs]
    else
        recovery_matrices.elem_sol .= Uresult.U[elem_dofs]
    end
end

function extract_element_solution!(
    recovery_matrices::PiezoelectricRecoveryMatrices,
    Uresult::NamedTuple,
    elem_dofs_u::AbstractVector{Int},
    elem_dofs_ϕ::AbstractVector{Int})
    
    recovery_matrices.elem_sol_u .= Uresult.U[elem_dofs_u]
    recovery_matrices.elem_sol_ϕ .= Uresult.ϕ[elem_dofs_ϕ]
end

# Main Recovery Function

function recover_field_values(
    materials::Vector{M},
    connectivity::Matrix{Int},
    nodes::Matrix{Float64},
    Uresult::NamedTuple,
    connect_elem_phase::Vector{Int},
    geom::GeometricData) where {M<:AbstractMaterial}
    
    Ne, NN = size(connectivity)
    num_nodes = size(nodes, 1)
    nthreads = Threads.nthreads()
    
    # Phase handling
    unique_phases = sort(unique(connect_elem_phase))
    phase_to_material_idx = Dict{Int, Int}(phase => i for (i, phase) in enumerate(unique_phases))
    material_indices_per_element = [phase_to_material_idx[p] for p in connect_elem_phase]
    
    # Create storage
    field_storage = create_field_storage(materials[1], num_nodes)
    
    # Thread-local working matrices
    recovery_matrices_per_thread = [[create_recovery_matrices(m, NN) for m in materials] for _ in 1:nthreads]
    
    # Node connectivity
    connect_nodes_elements, count_elem_per_node = find_elem_connected_to_nodes(connectivity, nodes[:,1])
    
    # DOF handling
    dofs_per_node = get_dofs_per_node(materials[1])
    dofs_buf_per_thread = [Vector{Int}(undef, NN * dofs_per_node) for _ in 1:nthreads]
    
    # For piezoelectric materials, need separate DOF buffers
    if M <: PiezoelectricMaterial
        mech_dofs_buf_per_thread = [Vector{Int}(undef, materials[1].dim * NN) for _ in 1:nthreads]
        elec_dofs_buf_per_thread = [Vector{Int}(undef, NN) for _ in 1:nthreads]
    end
    
    # Main node loop
    @inbounds Threads.@threads for i in 1:num_nodes
        tid = Threads.threadid()
        nb_elem = count_elem_per_node[i]
        connected_elements = connect_nodes_elements[i]
        
        # Initialize accumulators for this node
        if M <: ThermalMaterial
            flux_acc = zeros(materials[1].dim)
            grad_acc = zeros(materials[1].dim)
        elseif M <: ElasticMaterial
            nstr = rows_voigt(materials[1].dim, materials[1].symmetry)
            stress_acc = zeros(nstr)
            strain_acc = zeros(nstr)
        elseif M <: PiezoelectricMaterial
            nstr = rows_voigt(materials[1].dim, materials[1].symmetry)
            nelec = vec_rows(materials[1].dim, materials[1].symmetry)
            stress_acc = zeros(nstr)
            strain_acc = zeros(nstr)
            elec_disp_acc = zeros(nelec)
            elec_field_acc = zeros(nelec)
        elseif M <: PoroelasticMaterial
            nstr = rows_voigt(materials[1].dim, materials[1].symmetry)
            stress_acc = zeros(nstr)
            strain_acc = zeros(nstr)
            pressure_acc = 0.0
        end
        
        total_volume = 0.0
        
        for j in 1:nb_elem
            elem = connected_elements[j]
            mat_idx = material_indices_per_element[elem]
            mat = materials[mat_idx]
            recovery_matrices = recovery_matrices_per_thread[tid][mat_idx]
            
            conn_e = @view connectivity[elem, :]
            
            if M <: PiezoelectricMaterial
                # Handle piezoelectric DOFs
                local_mech_dofs!(mech_dofs_buf_per_thread[tid], conn_e, dofs_per_node, mat.dim)
                local_scal_dofs!(elec_dofs_buf_per_thread[tid], conn_e, dofs_per_node, mat.dim, 1)
                extract_element_solution!(recovery_matrices, Uresult, 
                                        mech_dofs_buf_per_thread[tid], elec_dofs_buf_per_thread[tid])
            else
                # Handle single-field materials
                dofs_buf = dofs_buf_per_thread[tid]
                for a in 1:NN
                    base = (conn_e[a] - 1) * dofs_per_node
                    cols = (a-1)*dofs_per_node+1 : a*dofs_per_node
                    dofs_buf[cols] .= base .+ (1:dofs_per_node)
                end
                extract_element_solution!(recovery_matrices, Uresult, dofs_buf)
            end
           
            
            # Recover fields for this element
            recover_element_fields!(
                recovery_matrices, mat, geom.differential_operator_matrices[elem],
                geom.shape_function_data, geom.jacobian_transformation_data, elem
            )
            
            # Accumulate results
            vol = recovery_matrices.volume_sum[]
            if vol > 0
                if M <: ThermalMaterial
                    flux_acc .+= recovery_matrices.flux_sum
                    grad_acc .+= recovery_matrices.grad_sum
                elseif M <: ElasticMaterial
                    stress_acc .+= recovery_matrices.stress_sum
                    strain_acc .+= recovery_matrices.strain_sum
                elseif M <: PiezoelectricMaterial
                    stress_acc .+= recovery_matrices.stress_sum
                    strain_acc .+= recovery_matrices.strain_sum
                    elec_disp_acc .+= recovery_matrices.elec_disp_sum
                    elec_field_acc .+= recovery_matrices.elec_field_sum
                elseif M <: PoroelasticMaterial
                    stress_acc .+= recovery_matrices.stress_sum
                    strain_acc .+= recovery_matrices.strain_sum
                    pressure_acc += recovery_matrices.pressure_sum[]
                end
                total_volume += vol
            end
        end
        
        # Finalize results for this node
        if total_volume > 0
            if M <: ThermalMaterial
                field_storage.flux[i, :] .= flux_acc ./ total_volume
                field_storage.grad_temp[i, :] .= grad_acc ./ total_volume
            elseif M <: ElasticMaterial
                field_storage.stress[i, :] .= stress_acc ./ total_volume
                field_storage.strain[i, :] .= strain_acc ./ total_volume
            elseif M <: PiezoelectricMaterial
                field_storage.stress[i, :] .= stress_acc ./ total_volume
                field_storage.strain[i, :] .= strain_acc ./ total_volume
                field_storage.elec_disp[i, :] .= elec_disp_acc ./ total_volume
                field_storage.elec_field[i, :] .= elec_field_acc ./ total_volume
            elseif M <: PoroelasticMaterial
                field_storage.stress[i, :] .= stress_acc ./ total_volume
                field_storage.strain[i, :] .= strain_acc ./ total_volume
                field_storage.pressure[i] = pressure_acc / total_volume
            end
        end
    end
    
    return field_storage
end

# Single material 
function recover_field_values(
    material::M,
    connectivity::Matrix{Int},
    nodes::Matrix{Float64},
    Uresult::NamedTuple,
    geom::GeometricData) where {M<:AbstractMaterial}
    
    Ne = size(connectivity, 1)
    connect_elem_phase = ones(Int, Ne)
    return recover_field_values([material], connectivity, nodes, Uresult, connect_elem_phase, geom)
end


# Utility Functions 

function find_elem_connected_to_nodes(connectivity::Matrix{Int}, node_coords::Vector{Float64})
    num_nodes = length(node_coords)
    num_elems = size(connectivity, 1)
    elements_nodes_elements = [Int[] for _ in 1:num_nodes]

    for elem in 1:num_elems
        for node in connectivity[elem, :]
            push!(elements_nodes_elements[node], elem)
        end
    end

    count_elem_per_node = [length(elems) for elems in elements_nodes_elements]
    return elements_nodes_elements, count_elem_per_node
end

