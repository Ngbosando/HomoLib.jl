# ==== Spectral Decomposition & Energies (Consistent Voigt Convention) ====

"""
    get_strain_matrix(ε)

    Convert input to strain matrix form following the convention:
    - **2D**: [σ] = C [ε₁₁, ε₂₂, 2ε₁₂]ᵀ → Returns 2×2 matrix
    - **3D**: [σ] = C [ε₁₁, ε₂₂, ε₃₃, 2ε₂₃, 2ε₁₃, 2ε₁₂]ᵀ → Returns 3×3 matrix

    Handles both:
    1. Voigt vectors (2D: 3-element, 3D: 6-element)
    2. Full strain matrices (2×2 or 3×3)
"""
function get_strain_matrix(ε::AbstractVector{<:Real})
    if length(ε) == 3  # 2D case
        return @SMatrix [ε[1]    ε[3]/2;  # Factor of 2 correction for ε₁₂
                        ε[3]/2   ε[2]]
    elseif length(ε) == 6  # 3D case
        return @SMatrix [ε[1]    ε[6]/2  ε[5]/2;  # ε₁₂, ε₁₃ → /2
                        ε[6]/2   ε[2]    ε[4]/2;  # ε₂₃ → /2
                        ε[5]/2   ε[4]/2   ε[3]]
    else
        throw(ArgumentError("Voigt vector must be length 3 (2D) or 6 (3D)"))
    end
end

function get_strain_matrix(ε::AbstractMatrix{<:Real})
    if size(ε) == (2,2)  # 2D
        return ε
    elseif size(ε) == (3,3)  # 3D
        return ε
    else
        throw(ArgumentError("Matrix must be 2×2 (2D) or 3×3 (3D)"))
    end
end

"""
    get_voigt_vector(ε)

    Convert to Voigt notation following the convention:
    - **2D**: Returns [ε₁₁, ε₂₂, 2ε₁₂]ᵀ
    - **3D**: Returns [ε₁₁, ε₂₂, ε₃₃, 2ε₂₃, 2ε₁₃, 2ε₁₂]ᵀ
"""
function get_voigt_vector(ε::AbstractMatrix{<:Real})
    if size(ε) == (2,2)  # 2D
        return [ε[1,1], ε[2,2], 2ε[1,2]]  # Engineering shear strain
    elseif size(ε) == (3,3)  # 3D
        return [ε[1,1], ε[2,2], ε[3,3], 
                2ε[2,3], 2ε[1,3], 2ε[1,2]]  # All shear terms ×2
    else
        throw(ArgumentError("Matrix must be 2×2 (2D) or 3×3 (3D)"))
    end
end

get_voigt_vector(ε::AbstractVector{<:Real}) = ε  # Assume already in correct Voigt form

# Macaulay brackets (optimized, type-stable)
macaulay_positive(x::Real) = ifelse(x > 0, x, zero(x))
macaulay_negative(x::Real) = ifelse(x < 0, x, zero(x))

"""
    voigt_inner_product(a, b)

    Compute the double contraction **a:b** following continuum mechanics conventions:
    - **2D**: a₁₁b₁₁ + a₂₂b₂₂ + a₁₂b₁₂ (where Voigt stores a₁₂ = 2ε₁₂)
    - **3D**: a₁₁b₁₁ + a₂₂b₂₂ + a₃₃b₃₃ + a₂₃b₂₃ + a₁₃b₁₃ + a₁₂b₁₂
"""
function voigt_inner_product(a, b)
    εa = get_voigt_vector(a)
    εb = get_voigt_vector(b)
    
    if length(εa) == 3  # 2D
        return εa[1]*εb[1] + εa[2]*εb[2] + (εa[3]*εb[3])/2  # 2's cancel for shear
    else  # 3D
        return εa[1]*εb[1] + εa[2]*εb[2] + εa[3]*εb[3] +
               (εa[4]*εb[4] + εa[5]*εb[5] + εa[6]*εb[6])/2  # All shear terms /2
    end
end

"""
    positive_negative_energy(ε, λ, μ)

    Compute (ψ⁺, ψ⁻) for a strain state in **2D or 3D**, following the Voigt convention.
    - **ε**: Strain 
    - **λ, μ**: Lamé constants
"""
function positive_negative_energy(ε, λ::Real, μ::Real)
    E = get_strain_matrix(ε)
    ev = eigen(Symmetric(E))  # Efficient symmetric eigensolver
    
    # Split eigenvalues into positive/negative parts
    λp = max.(ev.values, zero(eltype(E)))
    λm = min.(ev.values, zero(eltype(E)))
    
    # Reconstruct E⁺ and E⁻
    Epos = sum(λp[i] * (ev.vectors[:,i] * ev.vectors[:,i]') for i in eachindex(λp))
    Eneg = sum(λm[i] * (ev.vectors[:,i] * ev.vectors[:,i]') for i in eachindex(λm))
    
    # Convert back to Voigt notation
    εpos = get_voigt_vector(Epos)
    εneg = get_voigt_vector(Eneg)
    
    # Compute trace (2D: ε₁₁+ε₂₂, 3D: ε₁₁+ε₂₂+ε₃₃)
    trε = tr(E)
    
    # Compute energies
    ψp = 0.5λ * macaulay_positive(trε)^2 + μ * voigt_inner_product(εpos, εpos)
    ψm = 0.5λ * macaulay_negative(trε)^2 + μ * voigt_inner_product(εneg, εneg)
    
    return ψp, ψm
end

"""
    compare_energy(ε, λ, μ; tol=1e-12)

    Returns `true` if ψ⁺ > ψ⁻ + total.
"""
function compare_energy(ε, λ::Real, μ::Real; tol=1e-12)
    ψp, ψm = positive_negative_energy(ε, λ, μ)
    return ψp > ψm + tol
end