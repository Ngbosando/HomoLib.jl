# =============================================
# Precompute_data.jl 
# =============================================

# Data structures

struct FiniteElementShapeData{NN, NGP, D, T<:Float64}
    N::Vector{SVector{NN, T}}           # N_i(ξ_q) - Shape functions at Gauss points
    ∇N_ξ::Vector{SVector{NN, T}}        # ∂N/∂ξ - Derivatives wrt ξ
    ∇N_η::Vector{SVector{NN, T}}        # ∂N/∂η - Derivatives wrt η (empty in 1D)
    ∇N_ζ::Vector{SVector{NN, T}}        # ∂N/∂ζ - Derivatives wrt ζ (empty if D<3)
    weights::SVector{NGP, T}            # w_q - Integration weights
end

const SHAPE_DATA_CACHE = Dict{
    Tuple{Symbol, Int, Int, DataType, Int},
    FiniteElementShapeData{NN,NGP,D,Float64} where {NN,NGP,D}
}()

struct GeometricData{SD<:FiniteElementShapeData, JD, BD}
    shape_function_data::SD              # Shape function data
    jacobian_transformation_data::JD     # Matrix of tuples (detJ, invJ, J) - NGP × Ne
    differential_operator_matrices::BD   # Indexable differential operators
end

struct DifferentialOperatorStorage{T}
    voigt_gradient_operator        :: Union{Nothing, Array{T,4}}  # (nrowV × (D*NN) × NGP × Ne)
    vector_gradient_operator       :: Union{Nothing, Array{T,4}}  # (D × NN × NGP × Ne)
    scalar_interpolation_operator  :: Union{Nothing, Array{T,4}}  # (1 × NN × NGP × Ne)
end

"Element-indexable handle: B[e] => view on element e."
struct BIndexable{S}; storage::S; end
struct BAtElem{S};    storage::S; e::Int; end
struct BFieldAtElem{A}; A::A; e::Int; end

# Main precomputation function

function precompute_geometric_data(
    element_type::Symbol,
    int_order::Int,
    dim::Int,
    connectivity::Matrix{Int},
    nodes::Matrix{Float64},
    material::AbstractMaterial)
    
    # Compute shape function data
    shape_function_data = compute_shape_function_data(
        element_type, int_order, dim, size(connectivity, 2)
    )  

    # Compute Jacobians (Matrix{Tuple{detJ, invJ, J}} - NGP × Ne)
    jacobian_transformation_data = compute_element_jacobian_data(
        connectivity, nodes, shape_function_data, Val(dim)
    ) 

    # Compute B matrices (differential operators)
    differential_operator_matrices = compute_differential_operator_matrices(
        material, shape_function_data, jacobian_transformation_data
    ) 

    return GeometricData(shape_function_data, jacobian_transformation_data, differential_operator_matrices)
end

# Shape function computation

function compute_shape_function_data(
    element_type_sym::Symbol,
    int_order::Int,
    dim::Int,
    n_nodes::Int)

    key = (element_type_sym, int_order, dim, Float64, n_nodes)
    return get!(SHAPE_DATA_CACHE, key) do

        gauss_pts, gauss_w = integration_rule(element_type_sym, int_order, Val(dim))
        nqp = length(gauss_w)

        N_all     = Vector{SVector{n_nodes, Float64}}(undef, nqp)
        dN_dξ_all = Vector{SVector{n_nodes, Float64}}(undef, nqp)
        dN_dη_all = Vector{SVector{n_nodes, Float64}}(undef, dim >= 2 ? nqp : 0)
        dN_dζ_all = Vector{SVector{n_nodes, Float64}}(undef, dim == 3 ? nqp : 0)

        @inbounds for q in 1:nqp
            ξηζ = gauss_pts[q]
            N, dN_dξ, dN_dη, dN_dζ = shape_functions(element_type_sym, ξηζ...)  # (ξ,η[,ζ])
            N_all[q]     = N
            dN_dξ_all[q] = dN_dξ
            if dim >= 2; dN_dη_all[q] = dN_dη; end
            if dim == 3; dN_dζ_all[q] = dN_dζ; end
        end

        return FiniteElementShapeData{n_nodes, nqp, dim, Float64}(
            N_all, dN_dξ_all, dN_dη_all, dN_dζ_all,
            SVector{nqp, Float64}(gauss_w)
        )
    end
end

# Element Jacobian computation

function compute_element_jacobian_data(
    connectivity::Matrix{Int},
    nodes::Matrix{T},
    shp::FiniteElementShapeData{NN, NGP, ED, T},
    ::Val{D},) where {NN, NGP, ED, T<:Float64, D}

    Ne = size(connectivity, 1) 

    # NGP × Ne matrix of tuples (detJ, invJ, J)
    TupleType = Tuple{T, SMatrix{D,ED,T,D*ED}, SMatrix{D,ED,T,D*ED}}
    out = Matrix{TupleType}(undef, NGP, Ne) 

    coords = MMatrix{NN, D, T}(undef)

    @inbounds for e in 1:Ne
        out_e = view(out, :, e)

        compute_jacobian_at_gauss_points!(
            out_e, coords, view(connectivity, e, :), nodes, shp, Val(D), Val(ED)
        )
    end
    return out
end

function compute_jacobian_at_gauss_points!(
    out_for_elem,              
    coords::MMatrix{NN,D,T},  
    elem_conn::AbstractVector{Int},
    nodes::Matrix{T},
    shp::FiniteElementShapeData{NN, NGP, ED, T},
    ::Val{D}, ::Val{EDkw},) where {NN,NGP,D,ED,T<:Float64,EDkw}

    @inbounds for i in 1:NN, j in 1:D
        coords[i, j] = nodes[elem_conn[i], j]
    end
    coords = SMatrix(coords)  

    @inbounds for qp in 1:NGP
        J = calculate_jacobian_at_gauss_point(coords, shp, qp, Val(D), Val(EDkw))
        det_abs, invJ = if D == EDkw
            (abs(det(J)), inv(J))                        # Standard case
        else
            (sqrt(abs(det(J' * J))), pinv(J))           
        end
        out_for_elem[qp] = (det_abs, invJ, J)
    end
    
    return nothing
end

#Reference gradient computation

function reference_gradients(
    shp::FiniteElementShapeData{NN,NGP,ED,T},
    qp::Int,
    ::Val{1},) where {NN,NGP,ED,T<:Float64}

    M = MMatrix{1,NN,T}(undef)
    @inbounds M[1, :] .= shp.∇N_ξ[qp]
    return SMatrix{1,NN,T}(M)
end

function reference_gradients(
    shp::FiniteElementShapeData{NN,NGP,ED,T},
    qp::Int,
    ::Val{2},) where {NN,NGP,ED,T<:Float64}

    M = MMatrix{2,NN,T}(undef)
    @inbounds begin
        M[1, :] .= shp.∇N_ξ[qp]
        M[2, :] .= shp.∇N_η[qp]
    end
    return SMatrix{2,NN,T}(M)
end

function reference_gradients(
    shp::FiniteElementShapeData{NN,NGP,ED,T},
    qp::Int,
    ::Val{3},) where {NN,NGP,ED,T<:Float64}

    M = MMatrix{3,NN,T}(undef)
    @inbounds begin
        M[1, :] .= shp.∇N_ξ[qp]
        M[2, :] .= shp.∇N_η[qp]
        M[3, :] .= shp.∇N_ζ[qp]
    end
    return SMatrix{3,NN,T}(M)
end

# Jacobian computation: J = coords' * ∇N_ref'
function calculate_jacobian_at_gauss_point(
    coords::SMatrix{NN, D, T},
    shp::FiniteElementShapeData{NN, NGP, ED, T},
    qp::Int,
    ::Val{D},
    ::Val{EDkw},) where {NN,NGP,D,ED,T<:Float64,EDkw}

    ∇N_ref = reference_gradients(shp, qp, Val(EDkw))        # ::SMatrix{ED,NN,T}
    return (coords' * ∇N_ref')::SMatrix{D,ED,T}               # ::SMatrix{D,ED,T}
end

#B-matrices computation (ε(u), ∇ϕ, N_p, …)

function compute_differential_operator_matrices(
    material::AbstractMaterial,
    shp_data::FiniteElementShapeData{NN, NGP, D, T},
    jac_all) where {NN, NGP, D, T}

    n_elem = size(jac_all, 2)  # NGP × Ne 
    n_dofs = D * NN
    
    fields = material.B_types
    COMPUTE_STRAIN = (:strain in fields) || (:velocity_gradient in fields)
    COMPUTE_GRADIENT = (:temperature_gradient in fields) || (:electric_field in fields)
    
    rowsV = COMPUTE_STRAIN   ? rows_voigt(D,material.symmetry) : 0
    rowsG = COMPUTE_GRADIENT ? vec_rows(D,material.symmetry)   : 0
    
    Bv = COMPUTE_STRAIN   ? zeros(T, rowsV, n_dofs, NGP, n_elem) : nothing
    Bx = COMPUTE_GRADIENT ? zeros(T, rowsG, NN, NGP, n_elem) : nothing
    Bs = nothing  # scalar place holder

    stor = DifferentialOperatorStorage{T}(Bv, Bx, Bs)
    
    # Fill the matrices
    fill_B_matrices!(stor, jac_all, shp_data, material, Val(D))
    
    return BIndexable(stor)
end

function fill_B_matrices!(
    stor::DifferentialOperatorStorage{T},
    jac_all::AbstractMatrix{<:Tuple},
    shp_data::FiniteElementShapeData{NN, NGP, D, T},
    material::AbstractMaterial,
    ::Val{Dval}) where {NN, NGP, D, T, Dval}
    
    n_elem = size(jac_all, 2)
    
    COMPUTE_STRAIN = stor.voigt_gradient_operator !== nothing
    COMPUTE_GRADIENT = stor.vector_gradient_operator !== nothing
    
    scratch_phys = [MMatrix{D, NN, T}(undef) for _ in 1:Threads.nthreads()]
    scratch_ref  = [MMatrix{D, NN, T}(undef) for _ in 1:Threads.nthreads()]
    Threads.@threads for e in 1:n_elem
        tid = Threads.threadid()
        ∇N_phys = scratch_phys[tid]
        ∇N_ref  = scratch_ref[tid]
        for qp in 1:NGP
            invJ = jac_all[qp, e][2]
            compute_physical_gradients!(∇N_phys, ∇N_ref, shp_data, qp, invJ, Val(D))
            
            if COMPUTE_STRAIN
                element_strain_storage = view(stor.voigt_gradient_operator, :, :, qp, e)
                fill_strain_B!(element_strain_storage, ∇N_phys, material.symmetry, Val(D))
            end

            if COMPUTE_GRADIENT
                element_grad_storage = view(stor.vector_gradient_operator, :, :, qp, e)
                fill_gradient_B!(element_grad_storage, ∇N_phys, Val(D))
            end
        end
    end
    
    return nothing
end

#Physical gradient computation

function compute_physical_gradients!(
    ∇N_phys::MMatrix{D, NN, T},
    ∇N_ref::MMatrix{D, NN, T},
    shp_data::FiniteElementShapeData{NN, NGP, D, T},
    qp::Int,
    invJ::SMatrix{D, D, T},
    ::Val{D}) where {NN, NGP, D, T}
    
    if D == 1
        @inbounds ∇N_ref[1, :] .= shp_data.∇N_ξ[qp]
        @inbounds for i in 1:NN
            ∇N_ref[1, i] *= invJ[1, 1]
        end
    elseif D == 2
        @inbounds begin
            ∇N_ref[1, :] .= shp_data.∇N_ξ[qp]
            ∇N_ref[2, :] .= shp_data.∇N_η[qp]
        end
        # ∇N_phys = invJ' * ∇N_ref
        mul_AtB!(∇N_phys, invJ, ∇N_ref)
    else # D == 3

        @inbounds begin
            ∇N_ref[1, :] .= shp_data.∇N_ξ[qp]
            ∇N_ref[2, :] .= shp_data.∇N_η[qp]
            ∇N_ref[3, :] .= shp_data.∇N_ζ[qp]
        end

        mul_AtB!(∇N_phys, invJ, ∇N_ref)
    end
    
    return nothing
end

#Gradient B-matrix filling

function fill_gradient_B!(
    B_storage::AbstractMatrix{T},
    ∇N_phys::MMatrix{D, NN, T},
    ::Val{D}) where {T, D, NN}
    rowsG = size(B_storage, 1)
    if rowsG == 3
        @inbounds for i in 1:NN
            B_storage[1, i] = ∇N_phys[1, i]  # ∂/∂x
            B_storage[2, i] = ∇N_phys[2, i]  # ∂/∂y
        end
    else
        copyto!(B_storage, ∇N_phys)
    end
    return nothing
end

#Strain B-matrix filling

# 1D strain B-matrix
function fill_strain_B!(
    B_storage::AbstractMatrix{T},
    ∇N_phys::MMatrix{1, NN, T},
    symmetry::Symbol,
    ::Val{1}) where {T, NN}
    
    @inbounds B_storage[1, :] .= ∇N_phys[1, :]
    return nothing
end

# 2D strain B-matrix
function fill_strain_B!(
    B_storage::AbstractMatrix{T},
    ∇N_phys::MMatrix{2, NN, T},
    symmetry::Symbol,
    ::Val{2}) where {T, NN}
    
    rowsV = size(B_storage, 1)
    @inbounds for i in 1:NN
        ix, iy = 2i-1, 2i
      
        B_storage[1, ix] = ∇N_phys[1, i]  # εxx
        B_storage[1, iy] = 0.0
        B_storage[2, ix] = 0.0
        B_storage[2, iy] = ∇N_phys[2, i]  # εyy
        B_storage[3, ix] = ∇N_phys[2, i]  # γxy
        B_storage[3, iy] = ∇N_phys[1, i]
        
    end
    return nothing
end

# 3D strain B-matrix
function fill_strain_B!(
    B_storage::AbstractMatrix{T},
    ∇N_phys::MMatrix{3, NN, T},
    symmetry::Symbol,
    ::Val{3}) where {T, NN}
    
    @inbounds for i in 1:NN
        ix, iy, iz = 3i-2, 3i-1, 3i
        B_storage[1, ix] = ∇N_phys[1, i]  # εxx
        B_storage[1, iy] = 0.0
        B_storage[1, iz] = 0.0
        
        B_storage[2, ix] = 0.0
        B_storage[2, iy] = ∇N_phys[2, i]  # εyy
        B_storage[2, iz] = 0.0
        
        B_storage[3, ix] = 0.0
        B_storage[3, iy] = 0.0
        B_storage[3, iz] = ∇N_phys[3, i]  # εzz
        
        # Shear strains
        B_storage[4, ix] = 0.0
        B_storage[4, iy] = ∇N_phys[3, i]  # γyz
        B_storage[4, iz] = ∇N_phys[2, i]
        
        B_storage[5, ix] = ∇N_phys[3, i]  # γxz
        B_storage[5, iy] = 0.0
        B_storage[5, iz] = ∇N_phys[1, i]
        
        B_storage[6, ix] = ∇N_phys[2, i]  # γxy
        B_storage[6, iy] = ∇N_phys[1, i]
        B_storage[6, iz] = 0.0
    end
    return nothing
end

#BLAS routines for small matrices

function mul_AB!(Y, A::SMatrix{D,D,T}, X::AbstractMatrix{T}) where {D,T}
    @inbounds for i in 1:D, j in axes(X,2)
        s = zero(T)
        @simd for k in 1:D
            s += A[i,k] * X[k,j]
        end
        Y[i,j] = s
    end
    return Y
end

function mul_AtB!(Y, A::SMatrix{D,D,T}, X::AbstractMatrix{T}) where {D,T}
    @inbounds for i in 1:D, j in axes(X,2)
        s = zero(T)
        @simd for k in 1:D
            s += A[k,i] * X[k,j]   # A' * X
        end
        Y[i,j] = s
    end
    return Y
end

#Access to B[e][qp]

Base.getindex(B_indexable::BIndexable, e::Int) = BAtElem(B_indexable.storage, e)
Base.size(B_indexable::BIndexable) = (size(B_indexable.storage.voigt_gradient_operator, 4),)
Base.length(B_indexable::BIndexable) = size(B_indexable, 1)
Base.iterate(B_indexable::BIndexable, state=1) = state > length(B_indexable) ? nothing : (B_indexable[state], state+1)

# Functions return the actual views
voigt_gradient_operator(B::BAtElem) = B.storage.voigt_gradient_operator !== nothing ? view(B.storage.voigt_gradient_operator, :, :, :, B.e) : nothing
vector_gradient_operator(B::BAtElem) = B.storage.vector_gradient_operator !== nothing ? view(B.storage.vector_gradient_operator, :, :, :, B.e) : nothing
scalar_interpolation_operator(B::BAtElem) = B.storage.scalar_interpolation_operator !== nothing ? view(B.storage.scalar_interpolation_operator, :, :, :, B.e) : nothing

# Make BAtElem with property access syntax
Base.propertynames(::BAtElem) = (:voigt_gradient_operator, :vector_gradient_operator, :scalar_interpolation_operator)

function Base.getproperty(B::BAtElem, name::Symbol)
    if name === :voigt_gradient_operator
        return voigt_gradient_operator(B)
    elseif name === :vector_gradient_operator
        return vector_gradient_operator(B)
    elseif name === :scalar_interpolation_operator
        return scalar_interpolation_operator(B)
    else
        return getfield(B, name)
    end
end