# =============================================
# Material Definitions Core
# =============================================

# ----------------------
# Material Display
# ----------------------
function Base.show(io::IO, mat::Material)
    println(io, "Material physics: ", join(mat.type, " + "))
    println(io, "Dimension: ", mat.dim, "D")
    println(io, "Symmetry: ", mat.symmetry)

    println(io, "Properties:")
    for (k, v) in mat.properties
        println(io, "  $k = $v")
    end

    if !isnothing(mat.mass_properties)
        println(io, "Mass Properties:")
        for (k, v) in mat.mass_properties
            println(io, "  $k = $v")
        end
    end

    println(io, "Tensors:")
    for (i, T) in enumerate(mat.tensors)
        println(io, "  Tensor[$i]: size = ", size(T))
        show(IOContext(io, :limit => true), "text/plain", T)
        println(io)
    end

    println(io, "Conjugate fields: ", mat.B_types)
end

# ----------------------
# Constants
# ----------------------
const ALLOWED_SYMMETRIES = (
    :isotropic, :orthotropic, :anisotropic, :transversely_isotropic, :out_of_plane
)

const PHYSICS_TENSORS = Dict(
    :elastic           => [:C],       # stiffness
    :thermal           => [:κ],       # conductivity
    :electric          => [:ϵ],       # permittivity
    :fluid             => [:μ_f],     # fluid viscosity
    :pressure          => [:α_p]      # Biot coefficient for poroelastic (linear force)
)

const CONJUGATE_FIELDS = Dict(
    :elastic           => [:strain],
    :thermal           => [:temperature_gradient],
    :electric          => [:electric_field],
    :fluid             => [:velocity_gradient],
    :pressure          => [:pressure]
)

# ----------------------
# Main Material Definition
# ----------------------
function material_def(physics::Vector{Symbol}, dim::Int, symmetry::Symbol; 
                     mass_properties=nothing, kwargs...)
    dim ∈ (1, 2, 3) || error("Invalid dimension: $dim")
    symmetry ∈ ALLOWED_SYMMETRIES || error("Invalid symmetry: $symmetry")

    props = Dict{Symbol, Any}(kwargs...)
    physics = unique(copy(physics))
    tensors = AbstractMatrix[]
    B_types = Symbol[]

    # ---- Base Physics ----
    if :elastic in physics
        push!(tensors, elastic_tensor(dim, symmetry, props))
        push!(B_types, :strain)
    end
    if :electric in physics
        push!(tensors, permittivity_tensor(dim, symmetry, props))
        push!(B_types, :electric_field)
    end
    if :thermal in physics
        push!(tensors, conductivity_tensor(dim, symmetry, props))
        push!(B_types, :temperature_gradient)
    end
    if :fluid in physics
        push!(tensors, fluid_viscosity_tensor(dim, props))
        push!(B_types, :velocity_gradient)
    end

    # ---- Coupling Effects ----
    # Piezoelectric (elastic + electric)
    if (:elastic in physics) && (:electric in physics) && haskey(props, :e)
        push!(tensors, piezo_tensor(dim, symmetry, props))
    end

    # Thermal expansion requires elasticity
    if haskey(props, :α)
        :elastic ∉ physics && error("Thermal expansion requires elasticity")
        α_vec = thermal_expansion_tensor(dim, symmetry, props)
        β = -tensors[findfirst(isequal(elastic_tensor(dim, symmetry, props)), tensors)] * α_vec
        push!(tensors, reshape(β, :, 1))
    end
    # Stokes flow
    if :fluid in physics && haskey(props, :α_p) && :elastic ∉ physics
        α_p_vec  = poroelastic_tensor(dim, symmetry, props)
        push!(tensors, reshape(α_p_vec, :, 1))
        push!(B_types, :pressure)
    end
    
    # ---- Simple Poroelastic (pressure as force) ----
    if :elastic in physics && haskey(props, :α_p) && :fluid ∉ physics
        α_p_vec  = poroelastic_tensor(dim, symmetry, props)
        push!(tensors, reshape(α_p_vec, :, 1))
        # push!(B_types, :pressure)
    end

    # Biot poroelastic coupling (elastic + fluid)
    if (:elastic in physics) && (:fluid in physics) && haskey(props, :M) && haskey(props, :α_p)
        α_p_vec = poroelastic_tensor(dim, symmetry, props)
        invM_vec = fill(1/props[:M], length(α_p_vec))
        push!(tensors, reshape(α_p_vec, :, 1), reshape(invM_vec, :, 1))
    end

    return Material(physics, dim, symmetry, props, tensors, B_types, mass_properties)
end

# ----------------------
# Core Tensor Generators
# ----------------------

# Elastic Tensor
function elastic_tensor(dim::Int, symmetry::Symbol, props::Dict)
    # Check if full matrix is provided
    if haskey(props, :C) && props[:C] isa AbstractMatrix
        C = props[:C]
        # Validate size based on dimension and symmetry
        if dim == 1
            size(C) == (1,1) || error("C must be 1×1 in 1D")
        elseif dim == 2
            if symmetry == :out_of_plane
                size(C) == (4,4) || error("C must be 4×4 for 2D out-of-plane")
            else
                size(C) == (3,3) || error("C must be 3×3 in 2D")
            end
        else  # 3D
            size(C) == (6,6) || error("C must be 6×6 in 3D")
        end
        return C
    end
    
    # Otherwise generate automatically based on symmetry
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
            return G * [1 - ν ν 0 ν;
                        ν 1 - ν 0 ν;
                        0 0 (1 - 2ν) / 2 0;
                        ν ν 0 1 - ν]
        end
    else
        throw(ArgumentError("Unsupported symmetry: $symmetry"))
    end
end

# Permittivity Tensor
function permittivity_tensor(dim::Int, symmetry::Symbol, props::Dict)
    # Now just validates the provided matrix
    haskey(props, :ϵ) || error("Missing permittivity :ϵ (full matrix)")
    ϵ = props[:ϵ]
    if symmetry == :out_of_plane
        size(ϵ) == (3, 3) || error("ϵ must be 3×3 for out-of-plane in 2D")
    else
        size(ϵ) == (dim, dim) || error("ϵ must be $dim×$dim matrix")
    end

    issymmetric(ϵ) || @warn "Permittivity matrix is not symmetric"
    return ϵ
end


# Conductivity Tensor
function conductivity_tensor(dim::Int, symmetry::Symbol, props::Dict)
    # Check if full matrix is provided
    if haskey(props, :κ_matrix)
        κ = props[:κ_matrix]
        size(κ) == (dim, dim) || error("κ_matrix must be $dim×$dim")
        return κ
    end
    
    # Otherwise generate automatically
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
            length(κ_input) == 3 || error("κ must be length 3 for out-of-plane")
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

# Fluid Viscosity Tensor (Stokes creeping flow)
function fluid_viscosity_tensor(dim::Int, props::Dict)
    haskey(props, :μ_f) ||  error("Missing fluid viscosity :μ_f")
    μ_f = props[:μ_f]
    # simple isotropic viscous tensor for Newtonian fluid
    if dim == 2
        return μ_f * [2 0 0; 0 2 0; 0 0 1]
    elseif dim == 3
        return μ_f * [2 0 0 0 0 0;
                      0 2 0 0 0 0;
                      0 0 2 0 0 0;
                      0 0 0 1 0 0;
                      0 0 0 0 1 0;
                      0 0 0 0 0 1]
    end
    error("Unsupported dim for fluid viscosity")
end

# Piezoelectric Tensor
function piezo_tensor(dim::Int, symmetry::Symbol, props::Dict)
    # Now just validates the provided matrix
    haskey(props, :e) || error("Piezoelectric material requires :e parameter (full matrix)")
    e = props[:e]
    
    if dim == 2
        if symmetry == :out_of_plane
            size(e) == (3, 4) || error("Piezo matrix must be 3×4 for out-of-plane in 2D")
        else
            size(e) == (2, 3) || error("Piezo matrix must be 2×3 in 2D")
        end
    else
        size(e) == (3, 6) || error("Piezo matrix must be 3×6 in 3D")
    end
    return e
end

# Thermal Expansion Tensor
# function thermal_expansion_tensor(dim::Int, symmetry::Symbol, props::Dict)
#     haskey(props, :α) || error("Missing thermal expansion :α")
#     α = props[:α]
#     return α isa Number ? fill(α, (symmetry==:out_of_plane && dim==2) ? 4 : (dim==3 ? 6 : dim)) : α
# end
function thermal_expansion_tensor(dim::Int, symmetry::Symbol, props::Dict)
    # Check if full vector is provided
    if haskey(props, :α_vector)
        α = props[:α_vector]
        expected_len = dim == 2 ? (symmetry == :out_of_plane ? 4 : 3) : 6
        length(α) == expected_len || error("α_vector must be length $expected_len for $dim D $symmetry")
        return α
    end
    
    # Otherwise generate automatically
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
# Poroelastic Tensor (Biot coefficient)
function poroelastic_tensor(dim::Int, symmetry::Symbol, props::Dict)
    haskey(props, :α_p) || error("Missing Biot coefficient :α")
    α = props[:α_p]
    
    if symmetry == :isotropic
        if dim == 1
            return [α]
        elseif dim == 2
            return symmetry == :out_of_plane ? [α, α, 0,0] : [α, α,0]
        else
            return [α, α, α,0,0,0]
        end
    else
        return thermal_expansion_tensor(dim, symmetry, props)
    end
end