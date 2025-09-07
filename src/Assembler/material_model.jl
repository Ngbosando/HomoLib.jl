# =================
# Material system 
# =================

const TYPE_ELASTIC   = [:elastic]
const TYPE_THERMAL   = [:thermal]
const TYPE_PIEZO     = [:elastic, :electric]
const TYPE_PORO      = [:elastic, :pressure]
const TYPE_STOKES    = [:fluid, :pressure]

const BT_ELASTIC     = [:strain]
const BT_THERMAL     = [:temperature_gradient]
const BT_ELECTRIC    = [:electric_field]
const BT_PIEZO       = [:strain, :electric_field]
const BT_STOKES      = [:velocity_gradient, :pressure]

const EMPTY_PROPS = Dict{Symbol,Any}()  

"""
    ElasticMaterial(dim; symmetry=:isotropic, C=nothing, E=nothing, ν=nothing, ρ=0.0, stress=false, out_of_plane=false)

Elastic material.
- If "C" is provided (Voigt), it is used as is.
- Otherwise "E, ν" are used to build "C" with:
  - "stress=true" → plane stress
  - "stress=false" → plane strain - default
  - "out_of_plane=true" → 4×4 matrix with ε33
"""
struct ElasticMaterial{TC<:AbstractMatrix{<:Real}} <: AbstractMaterial
    dim::Int
    symmetry::Symbol
    C::TC
    ρ::Float64
    props::Dict{Symbol,Any}  
end

"""
    ThermalMaterial(dim; symmetry=:isotropic, κ=nothing, k=nothing, ρ=0.0, c=0.0)

Thermal material.
- If "κ" is provided (tensor), it is used as is.
- Otherwise "k" (isotropic scalar) is used to build "κ = k*I".
"""
struct ThermalMaterial{Tκ<:AbstractMatrix{<:Real}} <: AbstractMaterial
    dim::Int
    symmetry::Symbol
    κ::Tκ
    ρ::Float64
    c::Float64
    props::Dict{Symbol,Any}
end

"""
    PiezoelectricMaterial(dim; symmetry=:isotropic, C, ε, e, ρ=0.0, stress=false, out_of_plane=false)

Piezoelectric material.
- "C": elasticity tensor (Voigt)
- "ε": permittivity (dim×dim)
- "e": piezoelectric tensor (dim×nstr)
- "stress", "out_of_plane": same logic as ElasticMaterial if C built from E,ν
"""
struct PiezoelectricMaterial{TC<:AbstractMatrix{<:Real},Tε<:AbstractMatrix{<:Real},Te<:AbstractMatrix{<:Real}} <: AbstractMaterial
    dim::Int
    symmetry::Symbol
    C::TC
    ε::Tε
    e::Te
    ρ::Float64
    props::Dict{Symbol,Any}
end

"""
    PoroelasticMaterial(dim; symmetry=:isotropic, C=nothing, E=nothing, ν=nothing, α_p, M=nothing, ρs=0.0, ρf=0.0, stress=false, out_of_plane=false)

Poroelastic material (Biot).
- "C" or ("E","ν") for the matrix.
- "α_p": Biot coefficient (required)
- "M" (optional): storage modulus
- "stress", "out_of_plane": same logic as ElasticMaterial
"""
struct PoroelasticMaterial{TC<:AbstractMatrix{<:Real}} <: AbstractMaterial
    dim::Int
    symmetry::Symbol
    C::TC
    α_p::Union{Nothing,Float64}
    M::Union{Nothing,Float64}
    ρs::Float64
    ρf::Float64
    props::Dict{Symbol,Any}
end

"""
    StokesMaterial(dim; μ, ρ=1.0)

Fluid (Stokes): viscosity "μ", density "ρ".
"""
struct StokesMaterial <: AbstractMaterial
    dim::Int
    symmetry::Symbol
    μ::Float64
    ρ::Float64
    props::Dict{Symbol,Any}
end

voigt_nstr(dim::Int, out_of_plane::Bool=false) = begin
    if dim == 1
        return 1
    elseif dim == 2
        return out_of_plane ? 4 : 3
    elseif dim == 3
        return 6
    else
        error("Unsupported dim=$dim")
    end
end

function C_iso(dim::Int, symmetry::Symbol, E::Real, ν::Real, stress::Bool=false, out_of_plane::Bool=false)
    if dim == 1
        return reshape([E], 1, 1)
        
    elseif dim == 2
        if symmetry == :out_of_plane
            if stress
                E_mod = E / (1 - ν^2)
                C = zeros(4, 4)
                C[1,1] = 1; C[2,2] = 1
                C[1,2] = ν; C[2,1] = ν
                C[3,3] = (1 - ν) / 2
            
                C[1,4] = ν ; C[4,1] = ν 
                C[2,4] = ν ; C[4,2] = ν 
                C[4,4] = 1
                return E_mod * C
            else
                G = G = E /  ((1 + ν) * (1 - 2ν))
                C = [1 - ν ν 0 ν;
                            ν 1 - ν 0 ν;
                            0 0 (1 - 2ν) / 2 0;
                            ν ν 0 1 - ν]
                return G * C
            end
        else
            if stress
                E_mod = E / (1 - ν^2)
                return E_mod * [1.0 ν 0.0;
                                ν 1.0 0.0;
                                0.0 0.0 (1 - ν)/2]
            else
                G = E / ((1 + ν) * (1 - 2ν))
                return G * [1-ν ν 0.0;
                                ν 1-ν 0.0;
                                0.0 0.0 (1 - 2ν)/2]
            end
        end
        
    elseif dim == 3
        λ = E * ν / ((1 + ν) * (1 - 2ν))
        G = E / (2 * (1 + ν))
        return [λ+2G λ     λ     0  0  0;
                λ     λ+2G λ     0  0  0;
                λ     λ     λ+2G 0  0  0;
                0     0     0     G  0  0;
                0     0     0     0  G  0;
                0     0     0     0  0  G]
    else
        error("Unsupported dim=$dim")
    end
end

function k_iso(dim::Int, k::Real)
    K = zeros(dim,dim)
    @inbounds for i in 1:dim
        K[i,i] = k
    end
    return K
end

function ElasticMaterial(dim::Int; symmetry::Symbol=:isotropic, C=nothing, E=nothing, ν=nothing, ρ::Real=0.0, stress::Bool=false, out_of_plane::Bool=false)
    C_built, props = build_C_with_props(dim, symmetry, C, E, ν, stress, out_of_plane)
    ElasticMaterial(dim, symmetry, C_built, Float64(ρ), props)
end

function build_C_with_props(dim, symmetry, C, E, ν, stress, out_of_plane)
    props = Dict{Symbol,Any}()
    
    if dim == 2
        props[:stress] = stress
        props[:out_of_plane] = out_of_plane
    end

    if C !== nothing
        C_built = Matrix{Float64}(C)
        if E !== nothing; props[:E] = E; end
        if ν !== nothing; props[:ν] = ν; end
    elseif (E !== nothing) && (ν !== nothing)
        props[:E] = E
        props[:ν] = ν
        C_built = C_iso(dim, symmetry, Float64(E), Float64(ν), stress, out_of_plane)
    else
        error("Provide either C or (E, ν) for ElasticMaterial")
    end
    
    return C_built, props
end

ThermalMaterial(dim::Int; symmetry::Symbol=:isotropic, κ=nothing, k=nothing, ρ::Real=0.0, c::Real=0.0) =
    ThermalMaterial(dim, symmetry, build_κ(dim, κ, k), Float64(ρ), Float64(c), Dict{Symbol,Any}())

function build_κ(dim, κ, k)
    if κ !== nothing
        return Matrix{Float64}(κ)
    elseif k !== nothing
        return k_iso(dim, Float64(k))
    else
        error("Provide κ or k for ThermalMaterial")
    end
end

function PiezoelectricMaterial(dim::Int; symmetry::Symbol=:isotropic, C, ε, e, ρ::Real=0.0, stress::Bool=false, out_of_plane::Bool=false)
    props = Dict{Symbol,Any}()
    if dim == 2
        props[:stress] = stress
        props[:out_of_plane] = out_of_plane
    end
    PiezoelectricMaterial(dim, symmetry, Matrix{Float64}(C), Matrix{Float64}(ε), Matrix{Float64}(e), Float64(ρ), props)
end

function PoroelasticMaterial(dim::Int; symmetry::Symbol=:isotropic, C=nothing, E=nothing, ν=nothing,
                    α_p::Union{Nothing,Real}=nothing, M::Union{Nothing,Real}=nothing, ρs::Real=0.0, ρf::Real=0.0, stress::Bool=false, out_of_plane::Bool=false)
    C_built, props = build_C_with_props(dim, symmetry, C, E, ν, stress, out_of_plane)
    α_p !== nothing ? props[:α_p] = Float64(α_p) : nothing
    PoroelasticMaterial(dim, symmetry, C_built, α_p === nothing ? 0.0 : Float64(α_p),
                        M === nothing ? nothing : Float64(M), Float64(ρs), Float64(ρf), props)
end

StokesMaterial(dim::Int; symmetry::Symbol=:isotropic, μ::Real, ρ::Real=1.0) =
    StokesMaterial(dim, symmetry, Float64(μ), Float64(ρ), Dict{Symbol,Any}())

tensors_tuple(m::ElasticMaterial)      = (m.C,)
tensors_tuple(m::ThermalMaterial)      = (m.κ,)
tensors_tuple(m::PiezoelectricMaterial)= (m.C, m.ε, m.e)
tensors_tuple(m::PoroelasticMaterial)  = (m.C,)
tensors_tuple(m::StokesMaterial)       = ()

type_vec(::ElasticMaterial)       = TYPE_ELASTIC
type_vec(::ThermalMaterial)       = TYPE_THERMAL
type_vec(::PiezoelectricMaterial) = TYPE_PIEZO
type_vec(::PoroelasticMaterial)   = TYPE_PORO
type_vec(::StokesMaterial)        = TYPE_STOKES

B_types_vec(::ElasticMaterial)       = BT_ELASTIC
B_types_vec(::ThermalMaterial)       = BT_THERMAL
B_types_vec(::PiezoelectricMaterial) = BT_PIEZO
B_types_vec(::PoroelasticMaterial)   = BT_ELASTIC
B_types_vec(::StokesMaterial)        = BT_STOKES

massprops(m::ElasticMaterial)       = Dict{Symbol,Float64}(:ρ=>m.ρ)
massprops(m::ThermalMaterial)       = Dict{Symbol,Float64}(:ρ=>m.ρ, :c=>m.c)
massprops(m::PiezoelectricMaterial) = Dict{Symbol,Float64}(:ρ=>m.ρ)
massprops(m::PoroelasticMaterial)   = Dict{Symbol,Float64}(:ρs=>m.ρs, :ρf=>m.ρf, :α_p=>m.α_p)
massprops(m::StokesMaterial)        = Dict{Symbol,Float64}(:ρ=>m.ρ)

function Base.getproperty(m::AbstractMaterial, s::Symbol)
    if      s === :tensors       ; return tensors_tuple(m)
    elseif  s === :type         ; return type_vec(m)
    elseif  s === :B_types      ; return B_types_vec(m)
    elseif  s === :properties   ; return getfield(m, :props)
    elseif  s === :mass_properties
        props = getfield(m, :props)
        if !haskey(props, :massprops)
            props[:massprops] = massprops(m)
        end
        return props[:massprops]
    else
        props = getfield(m, :props)
        if haskey(props, s)
            return props[s]
        end
        
        return getfield(m, s)
    end
end

"""
    get_dofs_per_node(mat)

Number of DOFs per node based on physics.
"""
get_dofs_per_node(m::ElasticMaterial)        = m.dim
get_dofs_per_node(m::ThermalMaterial)        = 1
get_dofs_per_node(m::PiezoelectricMaterial)  = m.dim + 1
get_dofs_per_node(m::PoroelasticMaterial)    = m.dim
get_dofs_per_node(m::StokesMaterial)         = m.dim + 1

"""
    get_nstr(mat)
Number of strain components according to matrix C.
"""
get_nstr(m::AbstractMaterial) = hasfield(typeof(m), :C) ? size(m.C, 1) : (m.dim == 3 ? 6 : (m.dim == 2 ? 3 : 1))

function Base.show(io::IO, m::AbstractMaterial)
    phys = join(getproperty(m, :type), " + ")
    println(io, "Material physics: ", phys)
    println(io, "Dimension: ", m.dim, "D")
    println(io, "Symmetry: ", m.symmetry)

    if m.dim == 2
        props = getproperty(m, :properties)
        if get(props, :stress, false)
            if get(props, :out_of_plane, false)
                println(io, "Formulation: out-of-plane plane stress")
            else
                println(io, "Formulation: plane stress (σ₃₃=0)")
            end
        else
            if get(props, :out_of_plane, false)
                println(io, "Formulation: out-of-plane plane strain")
            else
                println(io, "Formulation: plane strain (ε₃₃=0)")
            end
        end
    end

    props = getproperty(m, :properties)
    if !isempty(props)
        println(io, "Properties:")
        for (k,v) in props
            k === :massprops && continue
            println(io, "  $k = $v")
        end
    end

    T = getproperty(m, :tensors)
    if !isempty(T)
        println(io, "Tensors:")
        for (i,A) in enumerate(T)
            println(io, "  T[$i] :: ", size(A))
        end
    end

    mp = getproperty(m, :mass_properties)
    if !isempty(mp)
        println(io, "Mass props:")
        for (k,v) in mp
            println(io, "  $k = $v")
        end
    end
end