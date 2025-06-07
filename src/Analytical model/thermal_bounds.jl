"""
    voigt_conductivity(volumes, conductivities)

Compute the Voigt upper bound for effective conductivity (arithmetic mean).
"""
# =============================================
# Voigt upper bound analytical model
# =============================================

function voigt_conductivity(volumes::Vector{<:Real}, conductivities::Vector{<:Real})
    return sum(v * σ for (v, σ) in zip(volumes, conductivities))
end

"""
    reuss_conductivity(volumes, conductivities)

Compute the Reuss lower bound for effective conductivity (harmonic mean).
"""
# =============================================
# Reuss lower bound analytical model
# =============================================

function reuss_conductivity(volumes::Vector{<:Real}, conductivities::Vector{<:Real})
    return 1.0 / sum(v / σ for (v, σ) in zip(volumes, conductivities))
end

"""
    hashin_shtrikman_bounds(σ1, v1, σ2, v2)

Compute the Hashin-Shtrikman bounds for a two-phase composite (σ_lower, σ_upper).
"""
# =============================================
# Hashin-Shtrikman bounds analytical model
# =============================================

function hashin_shtrikman_bounds(σ1::Real, v1::Real, σ2::Real, v2::Real)
    # Ensure σ1 <= σ2 to simplify calculations
    if σ1 > σ2
        σ1, σ2 = σ2, σ1
        v1, v2 = v2, v1
    end
    
    lower = σ1 + v2 / (1/(σ2 - σ1) + v1/(3σ1))
    upper = σ2 + v1 / (1/(σ1 - σ2) + v2/(3σ2))
    
    return (lower=lower, upper=upper)
end

"""
    maxwell_bounds(σ1, v1, σ2, v2)

Compute the Maxwell-Eucken effective conductivity bounds for a two-phase composite (lower, upper).
These coincide with Hashin-Shtrikman bounds for isotropic two-phase materials.
"""
# =============================================
# Maxwell-Eucken bounds analytical model
# =============================================

function maxwell_bounds(σ1::Real, v1::Real, σ2::Real, v2::Real)
    # When σ1 <= σ2, lower bound is matrix σ1 with inclusions σ2
    σ_low = σ1 * (σ2 + 2σ1 + 2v2*(σ2 - σ1)) / (σ2 + 2σ1 - v2*(σ2 - σ1))
    # Upper bound is matrix σ2 with inclusions σ1
    σ_high = σ2 * (σ1 + 2σ2 + 2v1*(σ1 - σ2)) / (σ1 + 2σ2 - v1*(σ1 - σ2))
    
    return (lower=σ_low, upper=σ_high)
end