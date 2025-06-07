# =============================================
# Material Definitions Core
# =============================================

function Base.show(io::IO, mat::Material)
    println(io, "Material physics: ", join(mat.type, " + "))
    println(io, "Dimension: ", mat.dim, "D")
    println(io, "Symmetry: ", mat.symmetry)
    println(io, "Properties:")
    for (k, v) in mat.properties
        println(io, "  $k = $v")
    end
    println(io, "Tensors:")
    for (i, T) in enumerate(mat.tensors)
        println(io, "  Tensor[$i]: size = ", size(T))
        show(IOContext(io, :limit => true), "text/plain", T)
        println(io)
    end
    println(io, "Conjugate fields: ", mat.B_types)
end

const ALLOWED_SYMMETRIES = (
    :isotropic, :orthotropic, :anisotropic, :transversely_isotropic, :out_of_plane
)

const PHYSICS_TENSORS = Dict(
    :elastic => [:C],
    :thermal => [:κ],
    :electric => [:ϵ],
    :thermal_expansion => [:α],
)

const CONJUGATE_FIELDS = Dict(
    :elastic => [:strain],
    :thermal => [:temperature_gradient],
    :electric => [:electric_field],
    :thermal_expansion => [:temperature],
)


function material_def(physics::Vector{Symbol}, dim::Int, symmetry::Symbol; kwargs...)
    props = Dict{Symbol, Any}(kwargs...)
    dim ∈ (1, 2, 3) || error("Invalid dimension: $dim")
    symmetry ∈ ALLOWED_SYMMETRIES || error("Invalid symmetry: $symmetry")

    physics = unique(copy(physics))
    tensors = AbstractMatrix[]
    B_types = Symbol[]

    # Elastic base
    if :elastic in physics
        push!(tensors, elastic_tensor(dim, symmetry, props))
        push!(B_types, :strain)
        println(" → [Elastic] tensor C added")
    end

    # Electric
    if :electric in physics
        push!(tensors, permittivity_tensor(dim, symmetry, props))
        push!(B_types, :electric_field)
        println(" → [Electric] permittivity tensor ϵ added")
    end

    # Thermal
    if :thermal in physics
        push!(tensors, conductivity_tensor(dim, symmetry, props))
        push!(B_types, :temperature_gradient)
        println(" → [Thermal] conductivity tensor κ added")
    end

    # Thermal expansion coupling
    if haskey(props, :α)
        :elastic ∉ physics && error("Thermal expansion requires elasticity")
        α = thermal_expansion_tensor(dim, symmetry, props)
        β = -tensors[1] * α
        push!(tensors, reshape(β, :, 1))
        # push!(B_types, :temperature)
        println(" → [Coupling] thermal expansion α → -Cα added")
    end

    # Viscous
    if :viscous in physics
        μ = get(props, :μ, nothing)
        λ = get(props, :λ, 0.0)
        if isnothing(μ)
            println(" ⚠ Missing μ for viscous physics")
        else
            D = viscous_tensor(dim, μ, λ)
            push!(tensors, D)
            push!(B_types, :strain_rate)
            println(" → [Viscous] viscosity tensor D added")
        end
    end

    # Viscoelasticity: C_inf = C + Cr
    if haskey(props, :Cr)
        Cr = props[:Cr] isa Adjoint ? Matrix(props[:Cr]) : props[:Cr]
        C_inf = tensors[1] + Cr
        push!(tensors, C_inf, Cr)
        println(" → [Coupling] viscoelasticity: C_inf and Cr added")
    end

    # Piezoelectricity
    if (:elastic in physics) && (:electric in physics) && (haskey(props, :e) || haskey(props, :d))
        e_tensor = piezo_tensor(dim, symmetry, props)
        push!(tensors, e_tensor)
        println(" → [Coupling] piezoelectricity tensor e added")
    end

    # Poroelasticity
    if haskey(props, :α_p) && haskey(props, :M)
        α_poro = reshape(poroelastic_tensor(dim, symmetry, props), :, 1)
        inv_M = reshape([1 / props[:M]], :, 1)
        push!(tensors, α_poro, inv_M)
        push!(B_types, :pressure)
        println(" → [Coupling] poroelastic α and 1/M added")

        if haskey(props, :K) && haskey(props, :μ_f)
            K, μ_f = props[:K], props[:μ_f]
            K isa AbstractMatrix || error(":K must be matrix")
            push!(tensors, K / μ_f)
            push!(B_types, :fluid_velocity)
            println(" → [Coupling] fluid flow tensor K/μ_f added")
        end
    end

    # Plasticity
    if (:elastic in physics) && haskey(props, :σ_y)
        σ_y = props[:σ_y]
        H = get(props, :H, 0.0)
        println(" → Coupling: Plasticity (σ_y = $σ_y, H = $H)")
    end

    # Viscoplasticity
    if (:elastic in physics) && (:viscous in physics) && haskey(props, :σ_y)
        println(" → [Viscoplastic] coupling detected")
        props[:viscoplastic] = true
    end

    return Material(physics, dim, symmetry, props, tensors, B_types)
end


# Helper function for property validation
function validate_property!(props, key, default=nothing)
    if !haskey(props, key)
        default === nothing && error("Missing required property: $key")
        props[key] = default
    end
    return props[key]
end

# =============================================
# Core Tensor Generators with Out-of-Plane Fixes
# =============================================
function elastic_tensor(dim::Int, symmetry::Symbol, props::Dict)
    if symmetry == :isotropic
        E, ν = props[:E], props[:ν]
        E > 0 || throw(ArgumentError("E must be positive"))
        (-1 < ν < 0.5) || throw(ArgumentError("ν out of range [0, 0.5)"))
        
        G = E / (2 * (1 + ν))
        if dim == 1
            return [E;;]
        elseif dim == 2
            plane_stress = get(props, :plane_stress, false)
            if plane_stress
                factor = E / (1 - ν^2)
                return factor * [1.0 ν 0.0;
                                ν 1.0 0.0;
                                0.0 0.0 (1 - ν)/2]
            else
                factor = E / ((1 + ν) * (1 - 2ν))
                return factor * [1-ν ν 0.0;
                                ν 1-ν 0.0;
                                0.0 0.0 (1 - 2ν)/2]
            end
        else  # 3D
            λ = E * ν / ((1 + ν) * (1 - 2ν))
            return [λ+2G λ     λ     0  0  0;
                    λ     λ+2G λ     0  0  0;
                    λ     λ     λ+2G 0  0  0;
                    0     0     0     G  0  0;
                    0     0     0     0  G  0;
                    0     0     0     0  0  G]
        end
        
    elseif symmetry == :orthotropic
        required = (:E1, :E2, :E3, :ν12, :ν13, :ν23, :G12, :G13, :G23)
        missing_params = [p for p in required if !haskey(props, p)]
        !isempty(missing_params) && throw(ArgumentError("Missing parameters: $(join(missing_params, ", "))"))
        
        E1, E2, E3 = props[:E1], props[:E2], props[:E3]
        ν12, ν13, ν23 = props[:ν12], props[:ν13], props[:ν23]
        G12, G13, G23 = props[:G12], props[:G13], props[:G23]
        
        # Create compliance matrix
        S = zeros(6, 6)
        S[1,1] = 1/E1; S[2,2] = 1/E2; S[3,3] = 1/E3
        S[1,2] = S[2,1] = -ν12/E1
        S[1,3] = S[3,1] = -ν13/E1
        S[2,3] = S[3,2] = -ν23/E2
        S[4,4] = 1/G23; S[5,5] = 1/G13; S[6,6] = 1/G12
        
        return inv(Symmetric(S))
        
    elseif symmetry == :anisotropic
        haskey(props, :C_matrix) || throw(ArgumentError("C_matrix required for anisotropic material"))
        C = props[:C_matrix]
        size(C) == (6,6) || throw(ArgumentError("C_matrix must be 6×6"))
        issymmetric(C) || throw(ArgumentError("C_matrix must be symmetric"))
        return C
        
    elseif symmetry == :transversely_isotropic
        required = (:E1, :E3, :ν12, :ν23, :G12)
        missing_params = [p for p in required if !haskey(props, p)]
        !isempty(missing_params) && throw(ArgumentError("Missing parameters: $(join(missing_params, ", "))"))
        
        E1, E3, ν12, ν23, G12 = props[:E1], props[:E3], props[:ν12], props[:ν23], props[:G12]
        D = (1 + ν12) * (1 - ν12 - 2*ν23^2*(E3/E1))
        D > 0 || throw(ArgumentError("Unstable parameter combination (D = $D ≤ 0)"))
        
        C11 = (E1*(1 - ν23^2*(E3/E1)))/D
        C12 = (E1*(ν12 + ν23^2*(E3/E1)))/D
        C13 = (E1*ν23*(1 + ν12))/D
        C33 = (E3*(1 - ν12^2))/D
        C44 = E3/(2*(1 + ν23))
        
        return [C11 C12 C13 0 0 0;
                C12 C11 C13 0 0 0;
                C13 C13 C33 0 0 0;
                0 0 0 C44 0 0;
                0 0 0 0 C44 0;
                0 0 0 0 0 G12]
                
    elseif symmetry == :out_of_plane
        E, ν = props[:E], props[:ν]
        G = E /  ((1 + ν) * (1 - 2ν))
        plane_stress = get(props, :plane_stress, false)
        if plane_stress
            return E / (1 - ν^2) * [1 ν 0 ν;
                                    ν 1 0 ν; 
                                    0 0 (1 - ν) / 2 0;
                                    ν ν 0 1]
        else
            return E / ((1 + ν) * (1 - 2ν)) * [1 - ν ν 0 ν;
                                                ν 1 - ν 0 ν;
                                                0 0 (1 - 2ν) / 2 0;
                                                ν ν 0 (1 - ν)]
        end
    else
        throw(ArgumentError("Unsupported symmetry: $symmetry"))
    end
end

function piezo_tensor(dim::Int, symmetry::Symbol, props::Dict)
    haskey(props, :e) ? e = props[:e] :
    haskey(props, :d) ? e = props[:d] * elastic_tensor(dim, symmetry, props) :
    error("Piezoelectric material requires :e or :d parameter")
    
    if dim == 2
        if symmetry == :out_of_plane
            size(e) == (4, 3) || error("Piezo matrix must be 4×3 for out-of-plane in 2D")
        else
            size(e) == (3, 2) || error("Piezo matrix must be 3×2 in 2D")
        end
    else
        size(e) == (6, 3) || error("Piezo matrix must be 6×3 in 3D")
    end
    return e
end

function permittivity_tensor(dim::Int, symmetry::Symbol, props::Dict)
    haskey(props, :ϵ) || error("Missing permittivity :ϵ")
    ϵ = props[:ϵ]
    
    # Unified handling for all symmetries (same as conductivity)
    if ϵ isa Number
        # Isotropic case
        return ϵ * I(dim)
    elseif ϵ isa AbstractVector
        if symmetry == :orthotropic
            length(ϵ) == dim || throw(ArgumentError("ϵ must be [ϵ_u, ϵ_v, ϵ_w] for orthotropic"))
            return dim == 2 ? Diagonal([ϵ[1], ϵ[2]]) : symmetry == :out_of_plane ?  Diagonal([ϵ[1], ϵ[2], 0]) : Diagonal([ϵ[1], ϵ[2], ϵ[3]])
        elseif symmetry == :transversely_isotropic
            length(ϵ) == 2 || throw(ArgumentError("ϵ must be [ϵ_p, ϵ_z] for transversely_isotropic"))
            return dim == 2 ? Diagonal([ϵ[1], ϵ[1]]) : Diagonal([ϵ[1], ϵ[1], ϵ[2]])
        else
            throw(ArgumentError("Vector input not supported for symmetry $symmetry"))
        end
    elseif ϵ isa AbstractMatrix
        size(ϵ) == (dim, dim) || throw(ArgumentError("ϵ matrix must be $dim×$dim"))
        issymmetric(ϵ) || throw(ArgumentError("ϵ matrix must be symmetric"))
        return ϵ
    else
        throw(ArgumentError("Invalid ϵ type: $(typeof(ϵ))"))
    end
end

function thermal_expansion_tensor(dim::Int, symmetry::Symbol, props::Dict)
    haskey(props, :α) || error("Missing thermal expansion :α")
    α = props[:α]
    
    if symmetry == :isotropic
        if dim == 1
            return [α]
        elseif dim == 2
            return symmetry == :out_of_plane ? [α, α, 0, 0] : [α, α, 0]
        else
            return [α, α, α, 0, 0, 0]
        end
        
    elseif symmetry == :orthotropic
        if α isa Number
            if dim == 2
                return symmetry == :out_of_plane ? [α, α, 0, 0] : [α, α, 0]
            else
                return [α, α, α, 0, 0, 0]
            end
        else
            if dim == 2
                length(α) == 2 || error("α vector must have 2 components for 2D orthotropic")
                return symmetry == :out_of_plane ? [α[1], α[2], 0, 0] : [α[1], α[2], 0]
            else
                length(α) == 3 || error("α vector must have 3 components for 3D orthotropic")
                return [α[1], α[2], α[3], 0, 0, 0]
            end
        end
        
    elseif symmetry == :transversely_isotropic
        if α isa Number
            if dim == 2
                return symmetry == :out_of_plane ? [α, α, 0, 0] : [α, α, 0]
            else
                return [α, α, α, 0, 0, 0]
            end
        else
            length(α) == 2 || error("Transverse isotropic α requires [α_p, α_z]")
            if dim == 2
                return symmetry == :out_of_plane ? [α[1], α[1], 0, 0] : [α[1], α[1], 0]
            else
                return [α[1], α[1], α[2], 0, 0, 0]
            end
        end
        
    elseif symmetry == :out_of_plane
        if α isa Number
            return dim == 2 ? [α, α, 0, 0] : error("Out-of-plane only for 2D")
        else
            dim == 2 && length(α) == 4 || error("Out-of-plane α must be scalar or length-4 vector in 2D")
            return α
        end
        
    elseif symmetry == :anisotropic
        expected_len = dim == 2 ? (symmetry == :out_of_plane ? 4 : 3) : 6
        α isa AbstractVector || error("Anisotropic α must be a vector")
        length(α) == expected_len || error("α must be length $expected_len for anisotropic $dim D")
        return α
        
    else
        error("Unsupported symmetry: $symmetry")
    end
end

function conductivity_tensor(dim::Int, symmetry::Symbol, props::Dict)
    haskey(props, :κ) || error("Missing conductivity :κ")
    κ_input = props[:κ]
    
    if κ_input isa Number
        # Isotropic case
        return κ_input * I(dim)
    elseif κ_input isa AbstractVector
        if symmetry == :orthotropic
            length(κ_input) == dim || throw(ArgumentError("κ must be [κ_u, κ_v, κ_w] for orthotropic"))
            return dim == 2 ? Diagonal([κ_input[1], κ_input[2]]) : Diagonal([κ_input[1], κ_input[2], κ_input[3]])
        elseif symmetry == :transversely_isotropic
            length(κ_input) == 2 || throw(ArgumentError("κ must be [κ_p, κ_z] for transversely_isotropic"))
            return dim == 2 ? Diagonal([κ_input[1], κ_input[1]]) : Diagonal([κ_input[1], κ_input[1], κ_input[2]])
        elseif symmetry == :out_of_plane
            length(κ_input) == 3 || error("κ must be length 3y for out-of-plane")
            return Diagonal(κ_input)
        else
            throw(ArgumentError("Vector input not supported for symmetry $symmetry"))
        end
    elseif κ_input isa AbstractMatrix
        size(κ_input) == (dim, dim) || throw(ArgumentError("κ matrix must be $dim×$dim"))
        issymmetric(κ_input) || throw(ArgumentError("κ matrix must be symmetric"))
        return κ_input
    else
        throw(ArgumentError("Invalid κ type: $(typeof(κ_input))"))
    end
end

function poroelastic_tensor(dim::Int, symmetry::Symbol, props::Dict)
    haskey(props, :α_p) || error("Missing Biot coefficient :α")
    α = props[:α_p]
    
    if symmetry == :isotropic
        if dim == 1
            return [α]
        elseif dim == 2
            return symmetry == :out_of_plane ? [α, α, 0,0] : [α, α, 0]
        else
            return [α, α, α, 0, 0, 0]
        end
    else
        return thermal_expansion_tensor(dim, symmetry, props)
    end
end

function viscous_tensor(dim::Int, μ::Real, λ::Real)
    if dim == 2
        return μ * [2 0 0; 0 2 0; 0 0 1] + λ * [1 0 0; 0 1 0; 0 0 0]
    elseif dim == 3
        return μ * [2 0 0 0 0 0;
                    0 2 0 0 0 0;
                    0 0 2 0 0 0;
                    0 0 0 1 0 0;
                    0 0 0 0 1 0;
                    0 0 0 0 0 1] +
               λ * [1 0 0 0 0 0;
                    0 1 0 0 0 0;
                    0 0 1 0 0 0;
                    0 0 0 0 0 0;
                    0 0 0 0 0 0;
                    0 0 0 0 0 0]
    end
    error("Unsupported dim for viscous tensor")
end