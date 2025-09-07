# =============================================
# Assembler/assemble.jl
# =============================================

include("Precompute_data.jl")
include("material_model.jl")
include("transform_boundary.jl")

const FESD{NN,NGP,D,T} = FiniteElementShapeData{NN,NGP,D,T}
const BElem{T} = BAtElem{<:DifferentialOperatorStorage{T}}
const BJacs    = AbstractMatrix{<:Tuple}

# DOF utilities 

dofs_u_per_node(mat::ElasticMaterial)         = mat.dim
dofs_u_per_node(mat::ThermalMaterial)         = 0
dofs_u_per_node(mat::PiezoelectricMaterial) = mat.dim
dofs_u_per_node(mat::PoroelasticMaterial)   = mat.dim

dofs_s_per_node(mat::ElasticMaterial)         = 0
dofs_s_per_node(mat::ThermalMaterial)         = 1
dofs_s_per_node(mat::PiezoelectricMaterial) = 1
dofs_s_per_node(mat::PoroelasticMaterial)   = 0

get_dofs_per_node(mat::AbstractMaterial) = dofs_u_per_node(mat) + dofs_s_per_node(mat)

function build_global_dofs(connectivity::Matrix{Int}, dofs_per_node::Int)
    Ne, NN = size(connectivity)
    gdofs = Matrix{Int}(undef, Ne, NN * dofs_per_node)
    dof_offsets = 1:dofs_per_node
    
    @inbounds for e in 1:Ne, a in 1:NN
        base = (connectivity[e, a] - 1) * dofs_per_node
        cols = (a-1)*dofs_per_node+1 : a*dofs_per_node
        gdofs[e, cols] .= base .+ dof_offsets
    end
    return gdofs
end

# DOFs field-wise (for bloc piezo)
function local_mech_dofs_fieldwise!(out::Vector{Int}, conn_row::AbstractVector{Int}, D::Int)
    NN = length(conn_row)
    @inbounds for a in 1:NN
        node = conn_row[a]
        base = (node-1)*D
        @simd for d in 1:D
            out[(a-1)*D + d] = base + d
        end
    end
    nothing
end

function local_scalar_dofs_fieldwise!(out::Vector{Int}, conn_row::AbstractVector{Int}, Nu::Int)
    NN = length(conn_row)
    @inbounds for a in 1:NN
        node = conn_row[a]
        out[a] = Nu + node
    end
    nothing
end

# DOFs / field
function local_mech_dofs!(out::Vector{Int}, conn_row, dofs_per_node::Int, D::Int)
    NN = length(conn_row)
    @inbounds for a in 1:NN
        base = (conn_row[a] - 1) * dofs_per_node
        idx_start = (a - 1) * D
        @simd for α in 1:D
            out[idx_start + α] = base + α
        end
    end
    nothing
end

function local_scal_dofs!(out::Vector{Int}, conn_row, dofs_per_node::Int, offset_u::Int, Ns::Int)
    NN = length(conn_row)
    @inbounds for a in 1:NN
        base = (conn_row[a] - 1) * dofs_per_node
        @simd for s in 1:Ns
            out[(a-1)*Ns + s] = base + offset_u + s
        end
    end
    nothing
end

function local_pressure_dofs!(dofs_buf::Vector{Int}, conn_e::AbstractVector{Int}, 
                             dofs_per_node::Int, dim::Int)
    # Extract pressure DOFs 
    for (i, node) in enumerate(conn_e)
        dofs_buf[i] = (node - 1) * dofs_per_node + dim + 1
    end
end

rows_voigt(D::Int, symmetry::Symbol) = (D==1 ? 1 : (D==2 ? 3 : 6)) + (symmetry === :out_of_plane)
vec_rows(D::Int, symmetry::Symbol)   = D + (symmetry === :out_of_plane ? 1 : 0)

nnz_K_per_elem(mat::ElasticMaterial, NN::Int)       = (mat.dim*NN)^2
nnz_K_per_elem(mat::ThermalMaterial, NN::Int)       = NN^2
nnz_K_per_elem(mat::PiezoelectricMaterial, NN::Int) = (mat.dim*NN + NN)^2
nnz_K_per_elem(mat::PoroelasticMaterial, NN::Int)       = (mat.dim*NN)^2

nnz_M_per_elem(mat::ElasticMaterial, NN::Int, flag::Bool)       = flag ? (mat.dim*NN)^2 : 0
nnz_M_per_elem(mat::ThermalMaterial, NN::Int, flag::Bool)       = flag ? NN^2 : 0
nnz_M_per_elem(mat::PiezoelectricMaterial, NN::Int, flag::Bool) = flag ? (mat.dim*NN)^2 : 0
nnz_M_per_elem(mat::PoroelasticMaterial, NN::Int, flag::Bool)       = flag ? (mat.dim*NN)^2 : 0

#  K/M structures 

mutable struct ElasticElementMatrices{T<:AbstractFloat}
    Kuu::Matrix{T}
    Muu::Union{Nothing,Matrix{T}}
    CB::Matrix{T}
end

mutable struct ThermalElementMatrices{T<:AbstractFloat}
    Kss::Matrix{T}
    Css::Union{Nothing,Matrix{T}}
    κG::Matrix{T}
end

mutable struct PiezoElementMatrices{T<:AbstractFloat}
    Kuu::Matrix{T}
    Kus::Matrix{T}
    Kss::Matrix{T}
    Muu::Union{Nothing,Matrix{T}}
    CB::Matrix{T}
    eTG::Matrix{T}
    εG::Matrix{T}
end

mutable struct PoroelasticElementMatrices{T<:AbstractFloat}
    Kuu::Matrix{T}
    Muu::Union{Nothing,Matrix{T}}
    CB::Matrix{T}
end

function allocate_element_matrices(mat::ElasticMaterial, NN::Int, T::DataType=Float64)::ElasticElementMatrices{T}
    D, nstr = mat.dim, rows_voigt(mat.dim, mat.symmetry)
    ElasticElementMatrices(zeros(T, D*NN, D*NN), nothing, zeros(T, nstr, D*NN))
end

function allocate_element_matrices(mat::ThermalMaterial, NN::Int, T::DataType=Float64)::ThermalElementMatrices{T}
    D = mat.dim
    rG = vec_rows(D, mat.symmetry)
    ThermalElementMatrices(zeros(T, NN, NN), nothing, zeros(T, rG, NN))
end

function allocate_element_matrices(mat::PiezoelectricMaterial, NN::Int, T::DataType=Float64)::PiezoElementMatrices{T}
    D  = mat.dim
    nstr = rows_voigt(D, mat.symmetry)
    rG = vec_rows(D, mat.symmetry)
    PiezoElementMatrices(zeros(T, D*NN, D*NN),
            zeros(T, D*NN, NN),
            zeros(T, NN, NN),
            nothing,
            zeros(T, nstr, D*NN),
            zeros(T, nstr, NN),
            zeros(T, rG,   NN))
end

function allocate_element_matrices(mat::PoroelasticMaterial, NN::Int, T::DataType=Float64)::PoroelasticElementMatrices{T}
    D, nstr = mat.dim, rows_voigt(mat.dim, mat.symmetry)
    PoroelasticElementMatrices(zeros(T, D*NN, D*NN), nothing, zeros(T, nstr, D*NN))
end

allocate_mass_matrix!(ElementMatrices::ElasticElementMatrices, D::Int, NN::Int) = ElementMatrices.Muu === nothing && (ElementMatrices.Muu = zeros(eltype(ElementMatrices.Kuu), D*NN, D*NN))
allocate_mass_matrix!(ElementMatrices::ThermalElementMatrices, NN::Int)         = ElementMatrices.Css === nothing && (ElementMatrices.Css = zeros(eltype(ElementMatrices.Kss), NN, NN))
allocate_mass_matrix!(ElementMatrices::PiezoElementMatrices, D::Int, NN::Int)   = ElementMatrices.Muu === nothing && (ElementMatrices.Muu = zeros(eltype(ElementMatrices.Kuu), D*NN, D*NN))

# Element stiffness Computation

function compute_element_matrices_elastic!(
    ElementMatrices::Union{ElasticElementMatrices, PoroelasticElementMatrices}, mat::Union{ElasticMaterial, PoroelasticMaterial}, B_e::BElem,
    shp::FESD{NN,NGP,D}, jacs::BJacs, e::Int; compute_mass::Bool=false) where {NN,NGP,D}
    
    fill!(ElementMatrices.Kuu, 0.0)
    C = mat.C
    
    if compute_mass
        allocate_mass_matrix!(ElementMatrices, D, NN)
        fill!(ElementMatrices.Muu, 0.0)
    end
    
    B_voigt_e = B_e.voigt_gradient_operator
    
    @inbounds for qp in 1:NGP
        B_qp = @view B_voigt_e[:, :, qp]       
        detJ= abs(jacs[qp, e][1]) 
        scaling = detJ * shp.weights[qp]
        mul!(ElementMatrices.CB, C, B_qp)
        mul!(ElementMatrices.Kuu, B_qp', ElementMatrices.CB, scaling, 1.0)
        
        if compute_mass
            N = shp.N[qp]
            m_factor = mat.mass_properties[:ρ] * scaling
            Muu = ElementMatrices.Muu
            
            @simd for jb in 1:NN
                j_base = (jb-1)*D
                for ib in 1:NN
                    i_base = (ib-1)*D
                    m = m_factor * N[ib] * N[jb]
                    for α in 1:D
                        Muu[i_base + α, j_base + α] += m
                    end
                end
            end
        end
    end
    return nothing
end

function compute_element_matrices_thermal!(
    ElementMatrices::ThermalElementMatrices, mat::ThermalMaterial, B_e::BElem,
    shp::FESD{NN,NGP,D}, jacs::BJacs, e::Int; compute_mass::Bool=false) where {NN,NGP,D}
    
    fill!(ElementMatrices.Kss, 0.0)
    κ = mat.κ

    if compute_mass
        allocate_mass_matrix!(ElementMatrices, NN)
        fill!(ElementMatrices.Css, 0.0)
    end
    
    G_e = B_e.vector_gradient_operator

    @inbounds for qp in 1:NGP
        detJ = abs(jacs[qp, e][1])
        scaling   = detJ * shp.weights[qp]
        G    = @view G_e[:, :, qp] 

        # K += G' * κ * G * scaling
        mul!(ElementMatrices.κG, κ, G)
        mul!(ElementMatrices.Kss, G', ElementMatrices.κG, scaling, 1.0)

        if compute_mass
            # Css += (ρ*c*dV) * N * N'
            N = shp.N[qp]
            mfac = mat.mass_properties[:ρ] * mat.c * scaling
            if N isa StridedVector{<:BlasFloat}
                BLAS.ger!(mfac, N, N, ElementMatrices.Css)
            else
                @inbounds for j in 1:NN
                    Nj = N[j]
                    @simd for i in 1:NN
                        ElementMatrices.Css[i,j] += mfac * N[i] * Nj
                    end
                end
            end
        end
    end
    return nothing
end

function compute_element_matrices_piezo!(
    ElementMatrices::PiezoElementMatrices, mat::PiezoelectricMaterial, B_e::BElem,
    shp::FESD{NN,NGP,D}, jacs::BJacs, e::Int; compute_mass::Bool=false) where {NN,NGP,D}
    
    fill!(ElementMatrices.Kuu, 0.0); fill!(ElementMatrices.Kus, 0.0); fill!(ElementMatrices.Kss, 0.0)
    C, ε, eT = mat.C, mat.ε, transpose(mat.e)

    if compute_mass
        allocate_mass_matrix!(ElementMatrices, D, NN)
        fill!(ElementMatrices.Muu, 0.0)
    end
    
    Bv_e = B_e.voigt_gradient_operator
    Gs_e = B_e.vector_gradient_operator

    @inbounds for qp in 1:NGP
        detJ = abs(jacs[qp, e][1])
        scaling = detJ * shp.weights[qp]
        
        Bv = @view Bv_e[:, :, qp]
        Gs = @view Gs_e[:, :, qp]
        
        # Stiffness sub-matrices 
        mul!(ElementMatrices.CB, C, Bv)
        mul!(ElementMatrices.Kuu, Bv', ElementMatrices.CB, scaling, 1.0)

        mul!(ElementMatrices.eTG, eT, Gs)
        mul!(ElementMatrices.Kus, Bv', ElementMatrices.eTG, scaling, 1.0)

        mul!(ElementMatrices.εG, ε, Gs)
        mul!(ElementMatrices.Kss, Gs', ElementMatrices.εG, scaling, 1.0)
        
        if compute_mass
            N = shp.N[qp]
            m_factor = mat.mass_properties[:ρ] * scaling
            Muu = ElementMatrices.Muu

            @simd for jb in 1:NN
                j_base = (jb-1)*D
                for ib in 1:NN
                    i_base = (ib-1)*D
                    m = m_factor * N[ib] * N[jb]
                    for α in 1:D
                        Muu[i_base + α, j_base + α] += m
                    end
                end
            end
        end
    end
    return nothing
end


# Assembly Logic 

function scatter_block!(
    I_global, J_global, V_global, ptr::Int,
    K_local::AbstractMatrix{T},
    g_rows::Vector{Int}, g_cols::Vector{Int},
    val_multiplier::T = one(T)) where {T}  

    n_rows, n_cols = size(K_local)
    @inbounds for j in 1:n_cols
        g_col_j = g_cols[j]
        base = ptr + (j - 1) * n_rows - 1
        @simd for i in 1:n_rows
            I_global[base + i] = g_rows[i]
            J_global[base + i] = g_col_j
            V_global[base + i] = K_local[i, j] * val_multiplier
        end
    end
    return n_rows * n_cols
end

function assemble_system_matrices!(
    I::NamedTuple{(:K, :M), Tuple{Vector{Int}, Vector{Int}}},
    J::NamedTuple{(:K, :M), Tuple{Vector{Int}, Vector{Int}}},
    V::NamedTuple{(:K, :M), Tuple{Vector{Float64}, Vector{Float64}}},
    connectivity::Matrix{Int}, nodes::Matrix{Float64},
    materials::Vector{<:AbstractMaterial}, spatial_dim::Int,
    geometric_data::GeometricData,
    material_indices_per_element::Vector{Int},
    compute_mass::Bool,
    compute_forces::Bool,
    boundary_data,  
    nodal_forces::Union{Nothing,Dict},
    boundary_faces::Union{Nothing,Dict},
    local_force_vectors::Dict{Symbol, Vector{Vector{Float64}}})

    Ne, NN = size(connectivity)
    reference_material = materials[1]
    dofs_per_node = get_dofs_per_node(reference_material)
    shape_functions, jacobians, B_operators = geometric_data.shape_function_data, geometric_data.jacobian_transformation_data, geometric_data.differential_operator_matrices
    D = spatial_dim
    Nu, Ns = dofs_u_per_node(reference_material), dofs_s_per_node(reference_material)

    elastic_matrices = nothing
    thermal_matrices = nothing 
    piezo_matrices = nothing
    poro_matrices = nothing

    # Local DOF 
    mech_dofs = Nu > 0 ? Vector{Int}(undef, spatial_dim * NN) : nothing
    scalar_dofs = Ns > 0 ? Vector{Int}(undef, NN) : nothing

    stiffness_ptr = 1
    mass_ptr = 1

    @views @inbounds for e in 1:Ne
        current_material = materials[material_indices_per_element[e]]
        element_connectivity = connectivity[e, :]

        if current_material isa ElasticMaterial
            elastic_matrices === nothing && (elastic_matrices = allocate_element_matrices(current_material, NN))
            compute_element_matrices_elastic!(elastic_matrices, current_material, B_operators[e], shape_functions, jacobians, e; compute_mass=compute_mass)
            local_mech_dofs!(mech_dofs, element_connectivity, dofs_per_node, D)

            # Stiffness matrix assembly
            stiffness_ptr += scatter_block!(I.K, J.K, V.K, stiffness_ptr, elastic_matrices.Kuu, mech_dofs, mech_dofs)

            # Mass matrix assembly
            if compute_mass
                allocate_mass_matrix!(elastic_matrices, spatial_dim, NN)
                mass_ptr += scatter_block!(I.M, J.M, V.M, mass_ptr, elastic_matrices.Muu, mech_dofs, mech_dofs)
            end

            # Force vector assembly
            if compute_forces    
                assemble_F_elem!(current_material, mech_dofs, nodes,
                                       jacobians, B_operators[e], shape_functions,
                                       dofs_per_node,
                                       local_force_vectors,
                                       nodal_forces, boundary_faces,
                                       boundary_data.jacobian_fcache, boundary_data.global_dofsforces, boundary_data.face_gauss_fdata,
                                       e)
            end

        elseif current_material isa ThermalMaterial
            thermal_matrices === nothing && (thermal_matrices = allocate_element_matrices(current_material, NN))
            compute_element_matrices_thermal!(thermal_matrices, current_material, B_operators[e], shape_functions, jacobians, e; compute_mass=compute_mass)
            local_scal_dofs!(scalar_dofs, element_connectivity, dofs_per_node, 0, 1)
            n_local = length(scalar_dofs)

            # Stiffness matrix assembly 
            stiffness_ptr += scatter_block!(I.K, J.K, V.K, stiffness_ptr, thermal_matrices.Kss, scalar_dofs, scalar_dofs)

            # Mass matrix assembly  
            if compute_mass
                for j in 1:n_local
                    global_col_j = scalar_dofs[j]
                    base = mass_ptr + (j-1)*n_local - 1
                    @simd for i in 1:n_local
                        I.M[base+i] = scalar_dofs[i]
                        J.M[base+i] = global_col_j
                        V.M[base+i] = (thermal_matrices.Css::Matrix{Float64})[i,j]
                    end
                end
                mass_ptr += n_local*n_local
            end

        elseif current_material isa PiezoelectricMaterial
            piezo_matrices === nothing && (piezo_matrices = allocate_element_matrices(current_material, NN))
            compute_element_matrices_piezo!(piezo_matrices, current_material, B_operators[e], shape_functions, jacobians, e; compute_mass=compute_mass)
            local_mech_dofs!(mech_dofs, element_connectivity, dofs_per_node, D)
            local_scal_dofs!(scalar_dofs, element_connectivity, dofs_per_node, D, 1)
            nu, ns = length(mech_dofs), length(scalar_dofs)

            # Kuu block - mechanical stiffness
            stiffness_ptr += scatter_block!(I.K, J.K, V.K, stiffness_ptr, piezo_matrices.Kuu, mech_dofs, mech_dofs)
            
            # Kus block - mechanical-electrical coupling
            stiffness_ptr += scatter_block!(I.K, J.K, V.K, stiffness_ptr, piezo_matrices.Kus, mech_dofs, scalar_dofs)
            
            # Ksu block - electrical-mechanical coupling (transpose of Kus)
            Ksu_transposed = transpose(piezo_matrices.Kus)
            stiffness_ptr += scatter_block!(I.K, J.K, V.K, stiffness_ptr, Ksu_transposed, scalar_dofs, mech_dofs, -1.0)
            
            # Kss block - electrical stiffness
            stiffness_ptr += scatter_block!(I.K, J.K, V.K, stiffness_ptr, piezo_matrices.Kss, scalar_dofs, scalar_dofs)
            
            # Mass matrix assembly (only mechanical)
            if compute_mass
                mass_ptr += scatter_block!(I.M, J.M, V.M, mass_ptr, piezo_matrices.Muu, mech_dofs, mech_dofs)
            end

            if compute_forces  

                combined_conn = Vector{Int}(undef, nu + ns)
                idx = 1
                @inbounds for node in 1:NN
                   
                    mstart = (node - 1) * D + 1
                    mend   = node * D
                    @inbounds for k in mstart:mend
                        combined_conn[idx] = mech_dofs[k]; idx += 1
                    end
                    
                    combined_conn[idx] = scalar_dofs[node]; idx += 1
                end
                assemble_F_elem!(current_material, combined_conn, nodes,
                                       jacobians, B_operators[e], shape_functions,
                                       dofs_per_node,
                                       local_force_vectors,
                                       nodal_forces, boundary_faces,
                                       boundary_data.jacobian_fcache, boundary_data.global_dofsforces, boundary_data.face_gauss_fdata,
                                       e)
            end
        else current_material isa PoroelasticMaterial
            if current_material.B_types == [:strain]
                elastic_matrices === nothing && (elastic_matrices = allocate_element_matrices(current_material, NN))
                compute_element_matrices_elastic!(elastic_matrices, current_material, B_operators[e], shape_functions, jacobians, e; compute_mass=compute_mass)
                local_mech_dofs!(mech_dofs, element_connectivity, dofs_per_node, D)

                # Stiffness matrix assembly
                stiffness_ptr += scatter_block!(I.K, J.K, V.K, stiffness_ptr, elastic_matrices.Kuu, mech_dofs, mech_dofs)

                # Mass matrix assembly
                if compute_mass
                    allocate_mass_matrix!(elastic_matrices, spatial_dim, NN)
                    mass_ptr += scatter_block!(I.M, J.M, V.M, mass_ptr, elastic_matrices.Muu, mech_dofs, mech_dofs)
                end

                # Force vector assembly
                if compute_forces    
                    assemble_F_elem!(current_material, mech_dofs, nodes,
                                        jacobians, B_operators[e], shape_functions,
                                        dofs_per_node,
                                        local_force_vectors,
                                        nodal_forces, boundary_faces,
                                        boundary_data.jacobian_fcache, boundary_data.global_dofsforces, boundary_data.face_gauss_fdata,
                                        e)
                end
            
            else
                # Place holder
            end
        end
    end
    return nothing
end

function assemble_KMF(
    connectivity::Matrix{Int}, nodes::Matrix{Float64},
    materials::Vector{<:AbstractMaterial}, dim::Int,
    connect_elem_phase::Vector{Int}, geom::GeometricData;
    compute_mass::Bool=false, NodalForces::Union{Nothing, Dict}=nothing,
    PointForces::Union{Nothing, Dict}=nothing,
    BoundaryFace=nothing)
    Ne, NN = size(connectivity)

    unique_phases = sort(unique(connect_elem_phase))
    if length(unique_phases) != length(materials)
        @warn "The number of unique phases ($(length(unique_phases))) does not match the number of materials provided ($(length(materials)))."
    end
    phase_to_material_idx = Dict{Int, Int}(phase => i for (i, phase) in enumerate(unique_phases))
    material_indices_per_element = [phase_to_material_idx[p] for p in connect_elem_phase]


    mat0 = materials[1]
    dofs_per_node = get_dofs_per_node(mat0)
    total_dofs = size(nodes, 1) * dofs_per_node
    
    nnzK = nnz_K_per_elem(mat0, NN) * Ne
    nnzM = compute_mass ? nnz_M_per_elem(mat0, NN, true) * Ne : 0
    compute_F = !isnothing(NodalForces) || !isnothing(BoundaryFace) || !isnothing(PointForces) ||
                any(should_apply_internal_force(m) for m in materials)
    compute_BoundaryFace = !isnothing(BoundaryFace)
    
    boundary_data = compute_BoundaryFace ? preprocess_boundary_faces(BoundaryFace, connectivity, materials) : 
                 (jacobian_fcache=nothing, global_dofsforces=nothing, face_element_map=nothing, face_gauss_fdata=nothing)

    I = (K=Vector{Int}(undef, nnzK), M=Vector{Int}(undef, nnzM))
    J = (K=Vector{Int}(undef, nnzK), M=Vector{Int}(undef, nnzM))
    V = (K=Vector{Float64}(undef, nnzK), M=Vector{Float64}(undef, nnzM))
    local_F = make_local_F(mat0, total_dofs, compute_F)

    assemble_system_matrices!(I, J, V, connectivity, nodes, materials, dim, geom, 
                      material_indices_per_element, compute_mass, compute_F, boundary_data,
                      NodalForces, BoundaryFace, local_F)

    K = sparse(I.K, J.K, V.K, total_dofs, total_dofs)
    M = compute_mass ? sparse(I.M, J.M, V.M, total_dofs, total_dofs) : spzeros(total_dofs, total_dofs)
    F = compute_F ? finalize_F(local_F, mat0, total_dofs, dofs_per_node, PointForces) : zeros(total_dofs)

    return (K=K, M=M, F=F)
end

# for single material case
function assemble_KMF(
    connectivity::Matrix{Int}, nodes::Matrix{Float64},
    material::AbstractMaterial, dim::Int, geom::GeometricData;
    compute_mass=false,
    NodalForces::Union{Nothing,Dict}=nothing,
    PointForces::Union{Nothing,Dict}=nothing,
    BoundaryFace=nothing)
    
    Ne = size(connectivity, 1)
    assemble_KMF(connectivity, nodes, [material], dim, ones(Int, Ne), geom;
                compute_mass=compute_mass, NodalForces=NodalForces, PointForces=PointForces, BoundaryFace=BoundaryFace)
end
 
# Force Computation

# Struct F 
abstract type AbstractForceMatrices end

abstract type AbstractFaceForceMatrices end

struct ElasticElementForceMatrices{DNN,D,T} <: AbstractForceMatrices
    f::MVector{DNN,T}
end

struct ThermalElementForceMatrices{NN,D,T} <: AbstractForceMatrices
    f::MVector{NN,T}
end

struct PiezoElementForceMatrices{DNN,NN,D,T} <: AbstractForceMatrices
    fu::MVector{DNN,T}
    fϕ::MVector{NN,T}
end

struct ElasticFaceMatrices{DNF,NF,D,T} <: AbstractFaceForceMatrices
    f::MVector{DNF,T}
    n::MVector{D,T}
    x::MVector{D,T}
    c::MVector{D,T}
end

struct ThermalFaceMatrices{NF,D,T} <: AbstractFaceForceMatrices
    f::MVector{NF,T}
    n::MVector{D,T}
    x::MVector{D,T}
    c::MVector{D,T}
end

struct PiezoFaceMatrices{DNF,NF,D,T} <: AbstractFaceForceMatrices
    fu::MVector{DNF,T}
    fϕ::MVector{NF,T}
    n::MVector{D,T}
    x::MVector{D,T}
    c::MVector{D,T}
end

function allocate_force_element_matrices(::Union{ElasticMaterial, PoroelasticMaterial}, ::Val{NN}, ::Val{D}, ::Type{T}=Float64) where {NN,D,T}
    DNN = D*NN
    ElasticElementForceMatrices{DNN,D,T}(MVector{DNN,T}(zeros(T, DNN)))
end

function allocate_force_element_matrices(::ThermalMaterial, ::Val{NN}, ::Val{D}, ::Type{T}=Float64) where {NN,D,T}
    ThermalElementForceMatrices{NN,D,T}(MVector{NN,T}(zeros(T, NN)))
end

function allocate_force_element_matrices(::PiezoelectricMaterial, ::Val{NN}, ::Val{D}, ::Type{T}=Float64) where {NN,D,T}
    DNN = D*NN
    PiezoElementForceMatrices{DNN,NN,D,T}(MVector{DNN,T}(zeros(T, DNN)), MVector{NN,T}(zeros(T, NN)))
end

function allocate_face_element_matrices(::Union{ElasticMaterial, PoroelasticMaterial}, ::Val{NF}, ::Val{D}, ::Type{T}=Float64) where {NF,D,T}
    DNF = D*NF
    ElasticFaceMatrices{DNF,NF,D,T}( MVector{DNF,T}(zeros(T, DNF)),
                               MVector{D,T}(zeros(T,D)), MVector{D,T}(zeros(T,D)), MVector{D,T}(zeros(T,D)) )
end

function allocate_face_element_matrices(::ThermalMaterial, ::Val{NF}, ::Val{D}, ::Type{T}=Float64) where {NF,D,T}
    ThermalFaceMatrices{NF,D,T}( MVector{NF,T}(zeros(T,NF)),
                           MVector{D,T}(zeros(T,D)), MVector{D,T}(zeros(T,D)), MVector{D,T}(zeros(T,D)) )
end

function allocate_face_element_matrices(::PiezoelectricMaterial, ::Val{NF}, ::Val{D}, ::Type{T}=Float64) where {NF,D,T}
    DNF = D*NF
    PiezoFaceMatrices{DNF,NF,D,T}( MVector{DNF,T}(zeros(T,DNF)), MVector{NF,T}(zeros(T,NF)),
                             MVector{D,T}(zeros(T,D)), MVector{D,T}(zeros(T,D)), MVector{D,T}(zeros(T,D)) )
end

# Struct Traction and force evaluation definition 

abstract type AbstractTraction end
struct ConstantTraction{T<:Real} <: AbstractTraction
    value::Vector{T}
end

struct FunctionTraction{F} <: AbstractTraction
    f::F
end

evaluate(::Nothing, ::AbstractVector) = nothing
evaluate(t::ConstantTraction, ::AbstractVector) = t.value
evaluate(t::FunctionTraction, x::AbstractVector) = t.f(x...)
force_computation(x) = x === nothing ? nothing :
                       x isa AbstractVector ? ConstantTraction(x) :
                       x isa Function ? FunctionTraction(x) :
                       error("Unsupported traction: $x")

function x_from_face!(x::AbstractVector{T}, N::SVector{NF,T}, X::AbstractMatrix{T}) where {NF,T}
    D = size(X,2)
    @inbounds for d in 1:D
        s = zero(T)
        for i in 1:NF; s += N[i] * X[i,d]; end
        x[d] = s
    end
    return x
end     

# Preprocessing boundary
function preprocess_boundary_faces(BoundaryFace,
                                   connectivity::Matrix{Int},
                                   materials::Vector{<:AbstractMaterial})
    jacobian_fcache = nothing
    global_dofsforces = nothing
    face_element_map = nothing
    gauss_fdata = nothing

    if !isnothing(BoundaryFace)
        jacobian_fcache = Dict{Symbol, AbstractMatrix{<:Tuple}}()
        global_dofsforces = Dict{Symbol, Matrix{Int}}()
        face_element_map = Dict{Symbol, Dict{Int, Int}}()
        for (name, bf) in BoundaryFace
            face_map = Dict{Int, Int}()
            for e in 1:size(connectivity, 1)
                for (idx, face) in enumerate(eachrow(bf.element_border))
                    if all(node -> node in connectivity[e, :], face)
                        face_map[e] = idx
                        break
                    end
                end
            end
            if !isempty(face_map)
                gauss_fdata = compute_shape_function_data(bf.element_type, bf.int_order, bf.dim-1, size(bf.element_border,2))
                ref_mat = materials[1]
                jacobian_fcache[name] = compute_element_jacobian_data(bf.element_border, bf.nodes, gauss_fdata, Val(bf.dim))
                global_dofsforces[name] = build_global_dofs(bf.element_border, get_dofs_per_node(ref_mat))
                face_element_map[name] = face_map
            end
        end
    end
    return (jacobian_fcache=jacobian_fcache,
            global_dofsforces=global_dofsforces,
            face_element_map=face_element_map,
            face_gauss_fdata=gauss_fdata)
end

# build stockage
function make_local_F(mat::AbstractMaterial, total_dofs::Int, compute_F::Bool)
    if !compute_F
        return Dict{Symbol, Vector{Vector{Float64}}}()
    end
    if length(mat.type) != 1 && mat.symmetry == :out_of_plane
        return Dict{Symbol, Vector{Vector{Float64}}}(
            bt => [zeros(total_dofs)]
            for bt in mat.B_types
        )
    else
        return Dict{Symbol, Vector{Vector{Float64}}}(
            :default => [zeros(total_dofs)]
        )
    end
end

# permutation function for couple system 
function get_node_wise_permutation(material::AbstractMaterial, n_nodes::Int)

    dofs_u = dofs_u_per_node(material)
    dofs_s = dofs_s_per_node(material)
    
    if dofs_u == 0 && dofs_s == 0
        return collect(1:n_nodes)
    end

    dofs_per_field_per_node = []
    field_types = Symbol[]

    if dofs_u > 0
        push!(dofs_per_field_per_node, dofs_u)
        push!(field_types, :mechanical)
    end
    if dofs_s > 0
        push!(dofs_per_field_per_node, dofs_s)
        push!(field_types, :scalar)
    end

    n_fields = length(dofs_per_field_per_node)
    dofs_per_node = sum(dofs_per_field_per_node)
    total_dofs = n_nodes * dofs_per_node
    p = zeros(Int, total_dofs)
    
    node_wise_field_offsets = cumsum(dofs_per_field_per_node)
    field_wise_offsets = cumsum([0; dofs_per_field_per_node]) .* n_nodes
    
    for i in 1:n_nodes
        node_offset_node_wise = (i - 1) * dofs_per_node
        
        for j in 1:n_fields
            dof_field = dofs_per_field_per_node[j]
            
            node_offset_field_wise = (i - 1) * dof_field
            field_offset_field_wise = field_wise_offsets[j]
            
            node_offset_in_node_wise = (j > 1) ? node_wise_field_offsets[j-1] : 0
            
            for k in 1:dof_field
                # Calculate the index in the node-wise ordered vector
                node_wise_idx = node_offset_node_wise + node_offset_in_node_wise + k
                
                # Calculate the corresponding index in the field-wise ordered vector
                field_wise_idx = field_offset_field_wise + node_offset_field_wise + k
                
                p[node_wise_idx] = field_wise_idx
            end
        end
    end
    
    return p
end

# Shape function and face normal 
function build_N_matrix(N_face::AbstractVector{T}, dim::Int) where {T<:Real}
    n_nodes = length(N_face)
    N = zeros(T, dim, dim * n_nodes)
    @inbounds for i in 1:n_nodes, d in 1:dim
        N[d, dim*(i-1)+d] = N_face[i]
    end
    return N
end

function face_normal!(n::AbstractVector{T}, J) where {T}
    D = length(n)
    if D == 2
       
        if size(J,1) == 2 && size(J,2) == 1
            @inbounds (tx = J[1,1]; ty = J[2,1])
        elseif size(J,1) == 1 && size(J,2) == 2
            @inbounds (tx = J[1,1]; ty = J[1,2])
        else
            @inbounds (tx = J[1,1]; ty = J[2,1])
        end
        n[1] =  ty
        n[2] = -tx
    else
        @inbounds begin
            t1x,t1y,t1z = J[1,1], J[2,1], J[3,1]
            t2x,t2y,t2z = J[1,2], J[2,2], J[3,2]
            n[1] = t1y*t2z - t1z*t2y
            n[2] = t1z*t2x - t1x*t2z
            n[3] = t1x*t2y - t1y*t2x
        end
    end
    s = sqrt(sum(n[i]^2 for i in eachindex(n)))
    if s > 0
        @inbounds @simd for i in eachindex(n)
            n[i] /= s
        end
    end
    return n
end

function face_center!(c::AbstractVector{T}, X::AbstractMatrix{T}) where {T}
    D = size(X,2); n = size(X,1)
    @inbounds for d in 1:D
        s = zero(T)
        for i in 1:n; s += X[i,d]; end
        c[d] = s / T(n)
    end
    return c
end
# Case of force application other than boundary Type
should_apply_internal_force(mat::AbstractMaterial)::Bool =
    mat.symmetry == :out_of_plane || haskey(mat.properties, :α) 

# Out-of-plane force computation
function add_out_of_plane!(f::AbstractVector{T}, B::AbstractMatrix{T}, tensor::AbstractMatrix{T}, scaling::T) where {T}
    M = zeros(eltype(f), size(tensor, 2))
    M[end] = 1.0
    B_M = tensor * M
    mul!(f, B', B_M, -scaling, 1.0)
    return nothing
end

# Thermal Dilatation
function add_thermal_expansion!(f::AbstractVector{T}, B::AbstractMatrix{T}, C::AbstractMatrix{T}, s::T) where {T}
    nstr = size(B,1); ndofs = size(B,2)
    @inbounds for j in 1:ndofs
        accj = zero(T)
        for i in 1:nstr
            vi = zero(T)
            for k in 1:(nstr-1) # ε_th[k]=1, ε_th[end]=0
                vi += C[i,k]
            end
            accj += B[i,j] * vi
        end
        f[j] += s * accj
    end
    return f
end

# Function to compute traction/volume 
function add_volume_mech!(f::AbstractVector{T}, N::AbstractVector{T}, s::T, b::AbstractVector{T}, D::Int) where {T}
    NN = length(N)
    @inbounds for i in 1:NN
        Ni = N[i]*s; base = (i-1)*D
        @simd for d in 1:D
            f[base+d] += Ni * b[d]
        end
    end
    f
end

add_volume_scalar!(f::AbstractVector{T}, N::AbstractVector{T}, s::T, q::T) where {T} =
    (@inbounds @simd for i in eachindex(N); f[i] += N[i]*q*s; end; f)

function compute_surface_force_vector(
    N_face::AbstractMatrix{T},
    scaling::T, fₜ::AbstractTraction,
    dim::Int, face_coords::Matrix{T}, n::AbstractVector{T}) where {T<:Real}

    # position x_qp
    n_nodes = size(N_face,2) ÷ dim
    x = zeros(T, dim)
    for d in 1:dim
        s = zero(T)
        for i in 1:n_nodes
            s += N_face[d, (i-1)*dim + d] * face_coords[i,d]
        end
        x[d] = s
    end
    τ = evaluate(fₜ, x)
    if τ isa Real
        τv = similar(x)
        @inbounds @simd for d in 1:dim; τv[d] = -T(τ)*n[d]; end
        return N_face' * τv * scaling
    else
        return N_face' * τ * scaling
    end
end

# Element force computation
function compute_element_force!(
    mat::Union{ElasticMaterial, PoroelasticMaterial},
    jac_all,                                    
    shp::FESD{NN,NGP,D,T},
    jac_fdata,
    shp_face::Union{Nothing,FiniteElementShapeData},
    forces,
    face_conn::Union{Nothing,AbstractVector{Int}},
    Bdict::BElem{T},
    face_coords::Union{Nothing,AbstractMatrix{T}},
    ElementMatrices,
    e::Int; face_idx::Int = 0,
    face_ElementMatrices::Union{Nothing, ElasticFaceMatrices} = nothing) where {NN,NGP,D,T<:Float64}
    
    # - FACE -
    if face_conn !== nothing
          
        @assert shp_face !== nothing && jac_fdata !== nothing && face_ElementMatrices !== nothing
        ff = face_ElementMatrices.f; n = face_ElementMatrices.n; x = face_ElementMatrices.x
        
        fill!(ff, zero(T))

        if forces !== nothing && !isnothing(forces.fₜ)
          
            fτ = forces.fₜ isa AbstractTraction ? forces.fₜ : force_computation(forces.fₜ)
            center = mean(face_coords, dims=1)
            @inbounds for qp in eachindex(shp_face.weights)
                Nf = shp_face.N[qp]
                detJ, _, J = jac_fdata[qp, face_idx]
                face_normal!(n, J)
                gauss_pos = Nf' * face_coords
                outward_check = dot(n, gauss_pos' - center[:]) < 0
                if !outward_check
                    n = -n 
                end
                scaling = abs(detJ) * shp_face.weights[qp]

                # build block-diagonal N for vector traction
                Nmat = build_N_matrix(Nf, D)
               
                # accumulate (dim*NF) vector
                ff .+= compute_surface_force_vector(Nmat, scaling, fτ, D, face_coords, n) 
            end
        end
        return ff
    end

    # VOLUME 
    f = ElementMatrices.f; fill!(f, zero(T))
    if forces !== nothing && !isnothing(forces.fᵥ)
       
        b = forces.fᵥ
        Bv_e = Bdict.voigt_gradient_operator
        @inbounds for qp in eachindex(shp.weights)
            detJ = (jac_all[qp,e])[1]::T
            scaling = detJ * shp.weights[qp]
            Bv = @view Bv_e[:, :, qp]
            
            if mat.symmetry == :out_of_plane
                add_out_of_plane!(f, Bv, mat.C, scaling)

            elseif haskey(mat.properties, :α)
                add_thermal_expansion!(f, Bv, mat.C, scaling)
            else
                add_volume_mech!(f, shp.N[qp], scaling, b, D)
            end
        end
    end

    return f
end

function compute_element_force!(
    mat::ThermalMaterial,
    jac_all::BJacs,
    shp::FESD{NN,NGP,D,T},
    jac_fdata::Union{Nothing,AbstractMatrix{<:Tuple}},
    shp_face::Union{Nothing,FiniteElementShapeData},
    forces,
    face_conn::Union{Nothing,AbstractVector{Int}},
    Bdict::BElem{T},
    face_coords::Union{Nothing,AbstractMatrix{T}},
    ElementMatrices::ThermalElementForceMatrices{NN,D,T},
    e::Int; face_idx::Int = 0,
    face_ElementMatrices::Union{Nothing,ThermalFaceMatrices} = nothing) where {NN,NGP,D,T<:Float64}

    # FACE
    if face_conn !== nothing
        @assert shp_face !== nothing && jac_fdata !== nothing && face_ElementMatrices !== nothing
        ff = face_ElementMatrices.f; n = face_ElementMatrices.n; x = face_ElementMatrices.x
        fill!(ff, zero(T))

        if forces !== nothing && !isnothing(forces.fₜ)
            fτ = forces.fₜ isa AbstractTraction ? forces.fₜ : force_computation(forces.fₜ)

            @inbounds for qp in eachindex(shp_face.weights)
                Nf = shp_face.N[qp]
                detJ, _, J = jac_fdata[qp, face_idx]
                face_normal!(n, J)
                scaling = abs(detJ) * shp_face.weights[qp]

                # evaluate scalar heat flux at x_qp
                x_from_face!(x, Nf, face_coords)
                q = evaluate(fτ, x)
                qv = q isa Real ? T(q) : T(q[1])

                # accumulate NF vector
                @inbounds @simd for i in eachindex(Nf)
                    ff[i] += Nf[i] * qv * s
                end
            end
        end
        return ff
    end

    # VOLUME 
    f = ElementMatrices.f; fill!(f, zero(T))
    if forces !== nothing && !isnothing(forces.fᵥ)
        qv = T(forces.fᵥ)
        Gs_e = Bdict.vector_gradient_operator
        @inbounds for qp in eachindex(shp.weights)
            detJ = (jac_all[qp,e])[1]::T
            scaling = abs(detJ) * shp.weights[qp]
            Gs = @view Gs_e[:, :, qp]

            if mat.symmetry == :out_of_plane
                add_out_of_plane!(f, Gs, mat.κ, scaling)
            else
                add_volume_scalar!(f, shp.N[qp], scaling, qv)
            end
        end
    end
    return f
end

function compute_element_force!(
    mat::PiezoelectricMaterial,
    jac_all,
    shp::FESD{NN,NGP,D,T},
    jac_fdata,
    shp_face::Union{Nothing,FiniteElementShapeData},
    forces,
    face_conn::Union{Nothing,AbstractVector{Int}},
    Bdict::BElem{T},
    face_coords::Union{Nothing,AbstractMatrix{T}},
    ElementMatrices::PiezoElementForceMatrices{DNN,NN,D,T},
    e::Int; face_idx::Int = 0,
    face_ElementMatrices::Union{Nothing,PiezoFaceMatrices} = nothing) where {NN,NGP,D,T,DNN}

    # Volume force calculation block
    if face_conn === nothing
        if mat.symmetry == :out_of_plane
            Bu_e = Bdict.voigt_gradient_operator
            Bϕ_e = Bdict.vector_gradient_operator
            C, e_mat, ε = mat.C, mat.e, mat.ε
    
            F_u, F_ϕ = zeros(T, DNN), zeros(T, NN)
            G_u, G_ϕ = zeros(T, DNN), zeros(T, NN)
    
            @inbounds for qp in eachindex(shp.weights)
                detJ = (jac_all[qp,e])[1]::T; scaling = abs(detJ) * shp.weights[qp]
                Bu = @view Bu_e[:, :, qp]; Bϕ = @view Bϕ_e[:, :, qp]
                
                # F vector contributions (related to epss33)
                add_out_of_plane!(F_u, Bu, C, scaling)
                add_out_of_plane!(F_ϕ, Bϕ, e_mat, scaling)
                
                # G vector contributions (related to E3)
                add_out_of_plane!(G_u, Bu, transpose(e_mat), -scaling)
                add_out_of_plane!(G_ϕ, Bϕ, ε, scaling)

            end
            
            dofs_per_node = D + 1
            f_F_interleaved = zeros(T, DNN + NN)
            f_G_interleaved = zeros(T, DNN + NN)
            
            for i in 1:NN
                base_u = (i-1)*D; base_comb = (i-1)*dofs_per_node
                for d in 1:D
                    f_F_interleaved[base_comb+d] = F_u[base_u+d]
                    f_G_interleaved[base_comb+d] = G_u[base_u+d]
                end
                f_F_interleaved[base_comb+D+1] = F_ϕ[i]
                f_G_interleaved[base_comb+D+1] = G_ϕ[i]
            end
            return (f_F_interleaved, f_G_interleaved)

        elseif !isnothing(forces) && !isnothing(forces.fᵥ) # Standard volume forces
            fu, fϕ = ElementMatrices.fu, ElementMatrices.fϕ
            fill!(fu, zero(T)); fill!(fϕ, zero(T))
            for qp in eachindex(shp.weights)
                detJ = (jac_all[qp,e])[1]::T; scaling = abs(detJ) * shp.weights[qp]
                add_volume_mech!(fu, shp.N[qp], scaling, forces.fᵥ, D)
                if hasproperty(forces, :qᵥ) && !isnothing(forces.qᵥ)
                    add_volume_scalar!(fϕ, shp.N[qp], scaling, T(forces.qᵥ))
                end
            end
            # Combine into a single vector
            combined_f = zeros(T, DNN + NN); dofs_per_node = D + 1
            for i in 1:NN
                base_u = (i-1)*D; base_comb = (i-1)*dofs_per_node
                for d in 1:D; combined_f[base_comb+d] = fu[base_u+d]; end
                combined_f[base_comb+D+1] = fϕ[i]
            end
            return combined_f
        end
    end

    # Surface force calculation block
    if face_conn !== nothing
        @assert face_ElementMatrices !== nothing && shp_face !== nothing
        ff_u, ff_ϕ = face_ElementMatrices.fu, face_ElementMatrices.fϕ
        n, xf = face_ElementMatrices.n, face_ElementMatrices.x
        fill!(ff_u, zero(T)); fill!(ff_ϕ, zero(T))

        if !isnothing(forces) && !isnothing(forces.fₜ)
            fτ = force_computation(forces.fₜ)
            for qp in eachindex(shp_face.weights)
                Nf = shp_face.N[qp]; detJ, _, J = jac_fdata[qp, face_idx]
                face_normal!(n, J); scaling = abs(detJ) * shp_face.weights[qp]
                Nmat = build_N_matrix(Nf, D)
                ff_u .+= compute_surface_force_vector(Nmat, scaling, fτ, D, face_coords, n)
            end
        end
        if !isnothing(forces) && hasproperty(forces, :f_q) && !isnothing(forces.f_q)
             f_q_traction = force_computation(forces.f_q)
             for qp in eachindex(shp_face.weights)
                Nf = shp_face.N[qp]; detJ, _, _ = jac_fdata[qp, face_idx]
                scaling = abs(detJ) * shp_face.weights[qp]
                x_from_face!(xf, Nf, face_coords)
                q_val = evaluate(f_q_traction, xf)
                qv = q_val isa Real ? T(q_val) : T(q_val[1])
                @inbounds @simd for i in eachindex(Nf); ff_ϕ[i] += Nf[i] * qv * scaling; end
            end
        end
        
        # Combine into a single vector for the face
        DNF, NF = length(ff_u), length(ff_ϕ)
        combined_f = zeros(T, DNF + NF); dofs_per_node_face = D + 1
        for i in 1:NF
            base_u = (i-1)*D; base_comb = (i-1)*dofs_per_node_face
            for d in 1:D; combined_f[base_comb+d] = ff_u[base_u+d]; end
            combined_f[base_comb+D+1] = ff_ϕ[i]
        end
        return combined_f
    end

    return mat.symmetry == :out_of_plane ? (zeros(T, DNN + NN), zeros(T, DNN + NN)) : zeros(T, DNN + NN)
end


# Point force computation
function apply_point_forces!(F::Vector{Float64}, PF::Dict{Int,<:AbstractVector{<:Real}}, dofs_per_node::Int)
    @inbounds for (node, v) in PF
        base = (node-1)*dofs_per_node
        @assert length(v) == dofs_per_node
        for i in 1:dofs_per_node
            F[base+i] += Float64(v[i])
        end
    end
    F
end

function apply_point_forces!(F::Vector{Float64}, PF::Dict{Int,<:Real}, dofs_per_node::Int)
    @inbounds for (node, s) in PF
        base = (node-1)*dofs_per_node
        F[base+1] += Float64(s)
    end
    F
end

                            
# assemble force
function assemble_F_elem!(mat::AbstractMaterial, conn, nodes::Matrix{Float64},
                          J_data::BJacs, B_dict::BElem, shp::FESD{NN,NGP,D},
                          dofs_per_node::Int,
                          local_F::Dict{Symbol, Vector{Vector{Float64}}},
                           NodalForces::Union{Nothing,Dict}, BoundaryFace::Union{Nothing,Dict},
                          jacobian_fcache, global_dofsforces, face_gauss_fdata,
                          e::Int) where {NN,NGP,D}

    # Surface Force Computation 
    
    if NodalForces !== nothing
        
        for (face_name, fsrc) in NodalForces
            isnothing(fsrc.fᵥ) && isnothing(fsrc.fₜ) && continue

            face_conn = Int[]
            jac_fcache = nothing
            face_coords = nothing
            face_idx = 0

            # Boundary face matching 
            if !isnothing(fsrc.fₜ) && BoundaryFace !== nothing && haskey(BoundaryFace, face_name)
                
                bf_faces = global_dofsforces[face_name]
                bf_jacs  = jacobian_fcache[face_name]
                
                for (i, face) in enumerate(eachrow(bf_faces))
                    if all(dof -> dof in conn, face)
                      
                        face_conn  = collect(face)
                        jac_fcache = bf_jacs
                        node_ids   = unique(floor.(Int, (face_conn .- 1) ./ dofs_per_node) .+ 1)
                        face_coords = nodes[node_ids, :]
                        face_idx = i
                        break
                    end
                end
            end

            isempty(face_conn) && continue

            face_conn_dofs = face_conn
            NF = length(face_conn_dofs) ÷ dofs_per_node 
            ElementMatrices  = allocate_force_element_matrices(mat, Val(NN), Val(D))
            f_ElementMatrices = allocate_face_element_matrices(mat, Val(NF), Val(D))

            fᵉ = compute_element_force!(
                mat,
                J_data, shp,
                jac_fcache, face_gauss_fdata,
                fsrc,
                face_conn,
                B_dict, face_coords,
                ElementMatrices, e; face_idx=face_idx, face_ElementMatrices=f_ElementMatrices)
            
            local_F[:default][1][face_conn] .+= fᵉ
        end
    end

    # Volume Force Computation
    if NodalForces === nothing || !any(f -> !isnothing(f.fᵥ), values(NodalForces))
        jac_fcache = nothing
        gauss_fdata_local = nothing
        face_conn = nothing
        face_coords = nothing
        
        default_f = should_apply_internal_force(mat) ?
                    (fᵥ = zeros(D), fₜ = nothing) :
                    (fᵥ = nothing,  fₜ = nothing)
        ElementMatrices = allocate_force_element_matrices(mat, Val(NN), Val(D))

        fᵉ = compute_element_force!(
            mat,
            J_data, shp,
            jac_fcache, gauss_fdata_local,
            default_f,
            face_conn,
            B_dict, face_coords,
            ElementMatrices, e)
            
        if length(mat.type) != 1 && mat.symmetry == :out_of_plane
            #  multiphysic case : tuple of fields
            for (i, bt) in enumerate(mat.B_types)
                local_F[bt][1][conn] .+= fᵉ[i]
            end

        else
            local_F[:default][1][conn] .+= fᵉ
        end
    end
    return nothing
end

# Force resizing
function finalize_F(local_F, ref_mat::AbstractMaterial, total_dofs::Int, dofs_per_node::Int, PointForces)
    if length(ref_mat.type) != 1 && ref_mat.symmetry == :out_of_plane
        final_forces = Dict{Symbol, Vector{Float64}}()
        for bt in ref_mat.B_types
            final_forces[bt] = reduce(+, local_F[bt])
        end
        
        if PointForces !== nothing && haskey(final_forces, :u)
            apply_point_forces!(final_forces[:u], PointForces, dofs_per_node)
        end
        
        return final_forces
        
    else 
        F = reduce(+, local_F[:default])
        
        if PointForces !== nothing
            apply_point_forces!(F, PointForces, dofs_per_node)
        end
        
        return F
    end
end

# Bloc part
# ============================================================
# 10) Assemblage multiphysique bloc (piézo)
#       -> Kuu, Kϕϕ, Kϕu (+ Muu, Fu, Fϕ)
# ============================================================
# function assemble_K_blocks_piezo(
#     connectivity::Matrix{Int}, nodes::Matrix{Float64},
#     materials::Vector{<:PiezoelectricMaterial}, dim::Int,
#     connect_elem_phase::Vector{Int}, geom::GeometricData;
#     compute_mass::Bool=false,
#     NodalForces::Union{Nothing, Dict}=nothing,
#     BoundaryFace=nothing,
#     kwargs...
# )
#     @assert all(m->m isa PiezoelectricMaterial, materials)
#     Ne, NN  = size(connectivity)
#     nnode   = size(nodes,1)
#     Nu      = nnode * dim
#     Nϕ      = nnode
#     shp, jacs, Bset = geom.shape_function_data, geom.jacobian_transformation_data, geom.differential_operator_matrices

#     unique_phases = sort(unique(connect_elem_phase))
#     phase_to_material_idx = Dict{Int, Int}(p => i for (i,p) in enumerate(unique_phases))
#     midx = [phase_to_material_idx[p] for p in connect_elem_phase]

#     # nnz par bloc (majorants simples)
#     nnz_uu = (dim*NN)^2 * Ne
#     nnz_up = (dim*NN)*NN * Ne
#     nnz_pp = (NN)^2 * Ne
#     nnz_M  = compute_mass ? (dim*NN)^2 * Ne : 0

#     Iuu = Vector{Int}(undef, nnz_uu); Juu = similar(Iuu); Vuu = Vector{Float64}(undef, nnz_uu)
#     Iup = Vector{Int}(undef, nnz_up); Jup = similar(Iup); Vup = Vector{Float64}(undef, nnz_up)   # K_{uϕ} (utilité: signe ensuite)
#     Ipu = Vector{Int}(undef, nnz_up); Jpu = similar(Ipu); Vpu = Vector{Float64}(undef, nnz_up)   # K_{ϕu}
#     Ipp = Vector{Int}(undef, nnz_pp); Jpp = similar(Ipp); Vpp = Vector{Float64}(undef, nnz_pp)

#     IM  = Vector{Int}(undef, nnz_M);  JM = similar(IM);  VM = Vector{Float64}(undef, nnz_M)

#     # ElementMatrices
#     mat0 = materials[1]
#     ElementMatricesP = allocate_element_matrices(mat0, NN)

#     # Forces bloc
#     Fu = zeros(Float64, Nu)
#     Fϕ = zeros(Float64, Nϕ)

#     # Face data
#     boundary_data = preprocess_boundary_faces(BoundaryFace, connectivity, materials, dim, true)

#     #  locaux field-wise
#     gU  = Vector{Int}(undef, dim*NN)
#     gPg = Vector{Int}(undef, NN)       # indices globaux φ (Nu + node)
#     gPl = Vector{Int}(undef, NN)       # indices locaux φ (1..Nϕ)

#     #  forces élémentaires
#     fu_e = zeros(Float64, dim*NN)
#     fp_e = zeros(Float64, NN)

#     pKuu = 1; pKup = 1; pKpu = 1; pKpp = 1; pM = 1

#     @views @inbounds for e in 1:Ne
#         mat = materials[midx[e]]
#         conn = connectivity[e,:]

#         # K local
#         compute_element_matrices_piezo!(ElementMatricesP, mat, Bset[e], shp, jacs, e; compute_mass)

#         # DOFs field-wise
#         local_mech_dofs_fieldwise!(gU, conn, dim)
#         local_scalar_dofs_fieldwise!(gPg, conn, Nu)
#         for i in 1:NN; gPl[i] = gPg[i] - Nu; end

#         nu, ns = length(gU), length(gPg)

#         # Kuu
#         for j in 1:nu
#             gUj = gU[j]; base = pKuu + (j-1)*nu - 1
#             @simd for i in 1:nu
#                 Iuu[base+i] = gU[i];  Juu[base+i] = gUj;  Vuu[base+i] = ElementMatricesP.Kuu[i,j]
#             end
#         end
#         pKuu += nu*nu

#         # K_{uϕ} (positif ici = B' e^T G) ; on renverra Kϕu et laisserons l’utilisateur faire Kuϕ = -Kϕu'
#         for j in 1:ns
#             gPj = gPg[j]; base = pKup + (j-1)*nu - 1
#             @simd for i in 1:nu
#                 Iup[base+i] = gU[i];  Jup[base+i] = gPj;  Vup[base+i] = ElementMatricesP.Kus[i,j]
#             end
#         end
#         pKup += nu*ns

#         # K_{ϕu} = (K_{uϕ})'
#         for j in 1:nu
#             gUj = gU[j]; base = pKpu + (j-1)*ns - 1
#             @simd for i in 1:ns
#                 Ipu[base+i] = gPg[i];  Jpu[base+i] = gUj;  Vpu[base+i] = ElementMatricesP.Kus[j,i]
#             end
#         end
#         pKpu += ns*nu

#         # Kϕϕ
#         for j in 1:ns
#             gPj = gPg[j]; base = pKpp + (j-1)*ns - 1
#             @simd for i in 1:ns
#                 Ipp[base+i] = gPg[i];  Jpp[base+i] = gPj;  Vpp[base+i] = ElementMatricesP.Kss[i,j]
#             end
#         end
#         pKpp += ns*ns

#         # Muu
#         if compute_mass
#             for j in 1:nu
#                 gUj = gU[j]; base = pM + (j-1)*nu - 1
#                 @simd for i in 1:nu
#                     IM[base+i] = gU[i]; JM[base+i] = gUj; VM[base+i] = (ElementMatricesP.Muu::Matrix{Float64})[i,j]
#                 end
#             end
#             pM += nu*nu
#         end

#         #---FORCES bloc-------
#         # volume + internes out_of_plane
#         fill!(fu_e, 0.0); fill!(fp_e, 0.0)
#         compute_element_force!(mat, jacs, shp, nothing, nothing,
#                                (fᵥ = nothing, f_qv = nothing, fₜ = nothing, f_q=nothing),
#                                nothing, conn, Bset[e], nothing, fu_e, fp_e, e)

#         # scatter vers Fu/Fϕ (bloc)
#         @simd for i in 1:nu; Fu[gU[i]] += fu_e[i]; end
#         @simd for i in 1:ns; Fϕ[gPl[i]] += fp_e[i]; end

#         # forces de surface (sur bords)
#         if NodalForces !== nothing
#             for (face_name, fsrc) in NodalForces
#                 (isnothing(fsrc.fᵥ) && isnothing(fsrc.fₜ) && !(hasproperty(fsrc,:f_q) && !isnothing(fsrc.f_q))) && continue
#                 bd = boundary_data
#                 if bd.face_element_map !== nothing &&
#                    haskey(bd.face_element_map, face_name) &&
#                    haskey(bd.face_element_map[face_name], e)

#                     face_idx = bd.face_element_map[face_name][e]
#                     jac_fcache = bd.jacobian_fcache[face_name]
#                     face_nodes = @view BoundaryFace[face_name].element_border[face_idx, :]
#                     Xface = @view nodes[face_nodes, :]
#                     shp_face = bd.face_gauss_fdata

#                     fill!(fu_e, 0.0); fill!(fp_e, 0.0)
#                     compute_element_force!(mat, jacs, shp, jac_fcache, shp_face,
#                                            fsrc, face_nodes, conn, Bset[e],
#                                            Xface, fu_e, fp_e, e; face_idx=face_idx)

#                     # indices globaux bloc pour la face
#                     NF = length(face_nodes)
#                     gUf = Vector{Int}(undef, dim*NF)
#                     local_mech_dofs_fieldwise!(gUf, face_nodes, dim)
#                     for k in 1:NF
#                         # phi local pour la face = node-id (1..Nϕ)
#                         Fϕ[face_nodes[k]] += fp_e[k]
#                     end
#                     @simd for i in 1:length(gUf); Fu[gUf[i]] += fu_e[i]; end
#                 end
#             end
#         end
#     end

#     # Matrices
#     Kuu = sparse(Iuu, Juu, Vuu, Nu, Nu)
#     Kϕu = sparse(Ipu, Jpu, Vpu, Nϕ, Nu)    # K_{ϕu}
#     Kϕϕ = sparse(Ipp, Jpp, Vpp, Nϕ, Nϕ)
#     Muu = compute_mass ? sparse(IM, JM, VM, Nu, Nu) : spzeros(Nu, Nu)

#     return (Kuu, Kϕϕ, Kϕu, Muu, Fu, Fϕ)
# end

# function assemble_K_blocks_piezo(
#     connectivity::Matrix{Int}, nodes::Matrix{Float64},
#     material::PiezoelectricMaterial, dim::Int, geom::GeometricData; kwargs...
# )
#     Ne = size(connectivity, 1)
#     assemble_K_blocks_piezo(connectivity, nodes, [material], dim, ones(Int,Ne), geom; kwargs...)
# end
# function assemble_K_blocks_piezo(
#     connectivity::Matrix{Int}, nodes::Matrix{Float64},
#     materials::Vector{<:PiezoelectricMaterial}, dim::Int,
#     connect_elem_phase::Vector{Int}, geom::GeometricData;
#     compute_mass::Bool=false,
#     NodalForces::Union{Nothing, Dict}=nothing,
#     BoundaryFace=nothing,
#     kwargs...
# )
#     @assert all(m->m isa PiezoelectricMaterial, materials)
#     Ne, NN  = size(connectivity)
#     nnode   = size(nodes,1)
#     Nu      = nnode * dim
#     Nϕ      = nnode
#     shp, jacs, Bset = geom.shape_function_data, geom.jacobian_transformation_data, geom.differential_operator_matrices

#     unique_phases = sort(unique(connect_elem_phase))
#     phase_to_material_idx = Dict{Int, Int}(p => i for (i,p) in enumerate(unique_phases))
#     midx = [phase_to_material_idx[p] for p in connect_elem_phase]

#     # nnz majorants
#     nnz_uu = (dim*NN)^2 * Ne
#     nnz_pu = (NN*dim*NN) * Ne
#     nnz_pp = (NN)^2 * Ne
#     nnz_M  = compute_mass ? (dim*NN)^2 * Ne : 0

#     # Triplets (we only keep Kuu, Kϕu, Kϕϕ; Kuϕ is not needed)
#     Iuu = Vector{Int}(undef, nnz_uu); Juu = similar(Iuu); Vuu = Vector{Float64}(undef, nnz_uu)
#     Ipu = Vector{Int}(undef, nnz_pu); Jpu = similar(Ipu); Vpu = Vector{Float64}(undef, nnz_pu)   # K_{ϕu}
#     Ipp = Vector{Int}(undef, nnz_pp); Jpp = similar(Ipp); Vpp = Vector{Float64}(undef, nnz_pp)
#     IM  = Vector{Int}(undef, nnz_M);  JM = similar(IM);  VM = Vector{Float64}(undef, nnz_M)

#     # ElementMatrices
#     mat0 = materials[1]
#     ElementMatricesP = allocate_element_matrices(mat0, NN)

#     # Forces bloc
#     Fu = zeros(Float64, Nu)
#     Fϕ = zeros(Float64, Nϕ)

#     # Face data
#     boundary_data = preprocess_boundary_faces(BoundaryFace, connectivity, materials, dim, true)

#     # Local DOF  (field-wise)
#     gU  = Vector{Int}(undef, dim*NN)   # 1..Nu
#     gPg = Vector{Int}(undef, NN)       # Nu + node
#     gPl = Vector{Int}(undef, NN)       # 1..Nϕ  (== node)

#     # Element force 
#     fu_e = zeros(Float64, dim*NN)
#     fp_e = zeros(Float64, NN)

#     pKuu = 1; pKpu = 1; pKpp = 1; pM = 1

#     @views @inbounds for e in 1:Ne
#         mat = materials[midx[e]]
#         conn = connectivity[e,:]

#         # Local K
#         compute_element_matrices_piezo!(ElementMatricesP, mat, Bset[e], shp, jacs, e; compute_mass)

#         # DOFs field-wise
#         local_mech_dofs_fieldwise!(gU, conn, dim)   # -> 1..Nu
#         local_scalar_dofs_fieldwise!(gPg, conn, Nu) # -> Nu + node
#         @inbounds for i in 1:NN; gPl[i] = gPg[i] - Nu; end  # -> 1..Nϕ

#         nu, ns = length(gU), length(gPl)

#         #--Kuu (Nu×Nu)-----
#         for j in 1:nu
#             gUj = gU[j]
#             base = pKuu + (j-1)*nu - 1
#             @simd for i in 1:nu
#                 Iuu[base+i] = gU[i]
#                 Juu[base+i] = gUj
#                 Vuu[base+i] = ElementMatricesP.Kuu[i,j]
#             end
#         end
#         pKuu += nu*nu

#         #--Kϕu (Nϕ×Nu)  rows=φ(local), cols=u(local)-----
#         for j in 1:nu
#             gUj = gU[j]
#             base = pKpu + (j-1)*ns - 1
#             @simd for i in 1:ns
#                 Ipu[base+i] = gPl[i]    # 1..Nϕ
#                 Jpu[base+i] = gUj       # 1..Nu
#                 Vpu[base+i] = ElementMatricesP.Kus[j,i]  # (K_{uϕ})' → K_{ϕu}
#             end
#         end
#         pKpu += ns*nu

#         #--Kϕϕ (Nϕ×Nϕ) rows=φ(local), cols=φ(local)-----
#         for j in 1:ns
#             gPj = gPl[j]
#             base = pKpp + (j-1)*ns - 1
#             @simd for i in 1:ns
#                 Ipp[base+i] = gPl[i]    # 1..Nϕ
#                 Jpp[base+i] = gPj       # 1..Nϕ
#                 Vpp[base+i] = ElementMatricesP.Kss[i,j]
#             end
#         end
#         pKpp += ns*ns

#         #--Muu (optionnel)-----
#         if compute_mass
#             for j in 1:nu
#                 gUj = gU[j]
#                 base = pM + (j-1)*nu - 1
#                 @simd for i in 1:nu
#                     IM[base+i] = gU[i]
#                     JM[base+i] = gUj
#                     VM[base+i] = (ElementMatricesP.Muu::Matrix{Float64})[i,j]
#                 end
#             end
#             pM += nu*nu
#         end

#         #--FORCES bloc (vol + internes out_of_plane)-----
#         fill!(fu_e, 0.0); fill!(fp_e, 0.0)
#         compute_element_force!(mat, jacs, shp, nothing, nothing,
#                                (fᵥ = nothing, f_qv = nothing, fₜ = nothing, f_q = nothing),
#                                nothing, conn, Bset[e], nothing, fu_e, fp_e, e)

#         @simd for i in 1:nu; Fu[gU[i]]  += fu_e[i]; end
#         @simd for i in 1:ns; Fϕ[gPl[i]] += fp_e[i]; end

#         #--Surface forces (si demandées)-----
#         if NodalForces !== nothing
#             for (face_name, fsrc) in NodalForces
#                 (isnothing(fsrc.fᵥ) && isnothing(fsrc.fₜ) && !(hasproperty(fsrc,:f_q) && !isnothing(fsrc.f_q))) && continue

#                 bd = boundary_data
#                 if bd.face_element_map !== nothing &&
#                    haskey(bd.face_element_map, face_name) &&
#                    haskey(bd.face_element_map[face_name], e)

#                     face_idx = bd.face_element_map[face_name][e]
#                     jac_fcache = bd.jacobian_fcache[face_name]
#                     face_nodes = @view BoundaryFace[face_name].element_border[face_idx, :]
#                     Xface = @view nodes[face_nodes, :]
#                     shp_face = bd.face_gauss_fdata

#                     fill!(fu_e, 0.0); fill!(fp_e, 0.0)
#                     compute_element_force!(mat, jacs, shp, jac_fcache, shp_face,
#                                            fsrc, face_nodes, conn, Bset[e],
#                                            Xface, fu_e, fp_e, e; face_idx=face_idx)

#                     # u (field-wise) for the face
#                     NF = length(face_nodes)
#                     gUf = Vector{Int}(undef, dim*NF)
#                     local_mech_dofs_fieldwise!(gUf, face_nodes, dim)
#                     @simd for i in 1:length(gUf)
#                         Fu[gUf[i]] += fu_e[i]
#                     end
#                     # φ (local indices = node ids)
#                     @inbounds for k in 1:NF
#                         Fϕ[face_nodes[k]] += fp_e[k]
#                     end
#                 end
#             end
#         end
#     end

#     # === Trim triplets to used length ===
#     Kuu = sparse(@view Iuu[1:pKuu-1], @view Juu[1:pKuu-1], @view Vuu[1:pKuu-1], Nu, Nϕ + Nu - Nϕ)  # Nu×Nu
#     Kuu = sparse(@view Iuu[1:pKuu-1], @view Juu[1:pKuu-1], @view Vuu[1:pKuu-1], Nu, Nu)

#     Kϕu = sparse(@view Ipu[1:pKpu-1], @view Jpu[1:pKpu-1], @view Vpu[1:pKpu-1], Nϕ, Nu)
#     Kϕϕ = sparse(@view Ipp[1:pKpp-1], @view Jpp[1:pKpp-1], @view Vpp[1:pKpp-1], Nϕ, Nϕ)
#     Muu = compute_mass ? sparse(@view IM[1:pM-1], @view JM[1:pM-1], @view VM[1:pM-1], Nu, Nu) : spzeros(Nu, Nu)

#     return (Kuu, Kϕϕ, Kϕu, Muu, Fu, Fϕ)
# end

# =============================================
# Blocs multiphysiques (πézo) — extract/rebuild
# =============================================

# tailles des blocs
# block_dims_piezo(nnode::Int, dim::Int) = (Nu = nnode*dim, Nϕ = nnode)

# # split RHS monolithique [u; φ]
# function split_rhs_piezo(F::AbstractVector{T}, nnode::Int, dim::Int) where {T}
#     Nu, Nϕ = block_dims_piezo(nnode, dim)
#     @assert length(F) == Nu + Nϕ "split_rhs_piezo: taille RHS incohérente."
#     Fu  = @view F[1:Nu]
#     Fϕ  = @view F[Nu+1:end]
#     return Fu, Fϕ
# end

# split K monolithique [u; φ] -> (Kuu, Kϕϕ, Kϕu)
# function split_blocks_piezo(K::SparseMatrixCSC{T,Int}, nnode::Int, dim::Int) where {T}
#     Nu, Nϕ = block_dims_piezo(nnode, dim)
#     @assert size(K,1) == size(K,2) == Nu+Nϕ "split_blocks_piezo: taille K incohérente."
#     Kuu = K[1:Nu,        1:Nu       ]
#     Kϕu = K[Nu+1:end,    1:Nu       ]
#     Kϕϕ = K[Nu+1:end,    Nu+1:end   ]
#     return Kuu, Kϕϕ, Kϕu
# end

# # variantes denses (si besoin)
# function split_blocks_piezo(K::AbstractMatrix{T}, nnode::Int, dim::Int) where {T}
#     Nu, Nϕ = block_dims_piezo(nnode, dim)
#     @assert size(K,1) == size(K,2) == Nu+Nϕ
#     Kuu = K[1:Nu,        1:Nu       ]
#     Kϕu = K[Nu+1:end,    1:Nu       ]
#     Kϕϕ = K[Nu+1:end,    Nu+1:end   ]
#     return Kuu, Kϕϕ, Kϕu
# end

# helpers internes pour composer K sparse à partir des blocs
# function append_block!(I::Vector{Int}, J::Vector{Int}, V::Vector{T},
#                         p::Int, K::SparseMatrixCSC{T,Int}, roff::Int, coff::Int) where {T}
#     rows, cols, vals = findnz(K)
#     n = length(vals)
#     @inbounds @simd for k in 1:n
#         I[p] = rows[k] + roff
#         J[p] = cols[k] + coff
#         V[p] = vals[k]
#         p += 1
#     end
#     return p
# end

# function append_transpose_block!(I::Vector{Int}, J::Vector{Int}, V::Vector{T},
#                                   p::Int, K::SparseMatrixCSC{T,Int}, roff::Int, coff::Int, sgn::T) where {T}
#     rows, cols, vals = findnz(K)
#     n = length(vals)
#     @inbounds @simd for k in 1:n
#         I[p] = cols[k] + roff
#         J[p] = rows[k] + coff
#         V[p] = sgn * vals[k]
#         p += 1
#     end
#     return p
# end

# build monolithique [u; φ] à partir des blocs (sparse)
# sign_convention = :yvonnet  => Kuφ = -Kϕu'
#                  = :symmetric => Kuφ = +Kϕu'
# function assemble_monolithic_from_blocks(Kuu::SparseMatrixCSC{T,Int},
#                                         Kϕϕ::SparseMatrixCSC{T,Int},
#                                         Kϕu::SparseMatrixCSC{T,Int};
#                                         sign_convention::Symbol = :yvonnet) where {T}
#     Nu = size(Kuu,1);  @assert size(Kuu,2) == Nu
#     Nϕ = size(Kϕϕ,1);  @assert size(Kϕϕ,2) == Nϕ
#     @assert size(Kϕu,1) == Nϕ && size(Kϕu,2) == Nu

#     nnz_total = nnz(Kuu) + nnz(Kϕϕ) + 2*nnz(Kϕu)
#     I = Vector{Int}(undef, nnz_total)
#     J = Vector{Int}(undef, nnz_total)
#     V = Vector{T}(undef, nnz_total)

#     p = 1
#     p = append_block!(I,J,V,p, Kuu, 0,   0)
#     p = append_block!(I,J,V,p, Kϕu, Nu,  0)  # bloc bas-gauche
#     sgn = sign_convention === :yvonnet ? -one(T) : one(T)
#     p = append_transpose_block!(I,J,V,p, Kϕu, 0,   Nu, sgn) # bloc haut-droit
#     p = append_block!(I,J,V,p, Kϕϕ, Nu,  Nu)
#     return sparse(I,J,V, Nu+Nϕ, Nu+Nϕ)
# end

# # build RHS monolithique [u; φ] depuis (Fu, Fϕ)
# function assemble_rhs_from_blocks(Fu::AbstractVector{T}, Fϕ::AbstractVector{T}) where {T}
#     return vcat(Fu, Fϕ)
# end

# # check rapide de consistance (optionnel)
# function check_block_shapes(Kuu, Kϕϕ, Kϕu, nnode::Int, dim::Int)
#     Nu, Nϕ = block_dims_piezo(nnode, dim)
#     @assert size(Kuu,1) == size(Kuu,2) == Nu
#     @assert size(Kϕϕ,1) == size(Kϕϕ,2) == Nϕ
#     @assert size(Kϕu,1) == Nϕ && size(Kϕu,2) == Nu
#     return true
# end
