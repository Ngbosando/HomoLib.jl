"""
    theorical_bound(field_type, bound_type, args...)

Compute theoretical bounds for effective material properties of composite materials.

# Arguments
- {field_type::Symbol}: Type of physical field (`:elastic` or `:thermal`)
- {bound_type::Symbol}: Type of bound to compute:
  - For elastic: {:voigt}, {:reuss}, {:hashin_shtrikman}
  - For thermal: {:voigt}, {:reuss}, {:hashin_shtrikman}, {:maxwell}
- {args...}: Material properties and volume fraction:
  - Elastic: {(κ1, μ1, κ2, μ2, f)}
  - Thermal: {(k1, k2, f)}

# Returns
- For elastic field:
  - {:voigt}/{:reuss}: Named tuple (κ=..., μ=...)
  - {:hashin_shtrikman}: Named tuple (κ_lower=..., κ_upper=..., μ_lower=..., μ_upper=...)
- For thermal field:
  - {:voigt}/{:reuss}: Named tuple (k=...)
  - {:hashin_shtrikman}/{:maxwell}: Named tuple (k_lower=..., k_upper=...)

# Theory
    Computes various bounds for effective properties of two-phase composites:
    - "Voigt bound": Arithmetic mean (upper bound for stiffness)
    - "Reuss bound": Harmonic mean (lower bound for stiffness)
    - "Hashin-Shtrikman": Tighter bounds for isotropic composites
    - "Maxwell's approximation": Effective conductivity for dilute suspensions

    # Examples
    ```julia
    # Elastic bounds
    elastic_voigt = theorical_bound(:elastic, :voigt, 3.0, 1.5, 10.0, 5.0, 0.3)
    hs_bounds = theorical_bound(:elastic, :hashin_shtrikman, 3.0, 1.5, 10.0, 5.0, 0.3)

    # Thermal bounds
    thermal_reuss = theorical_bound(:thermal, :reuss, 0.5, 5.0, 0.2)
    maxwell_bounds = theorical_bound(:thermal, :maxwell, 0.5, 5.0, 0.3)

"""

function theorical_bound(field_type::Symbol, bound_type::Symbol, args...)
    !(field_type in [:elastic, :thermal]) && 
        error("Unsupported field type: $field_type. Use :elastic or :thermal.")
    
    if field_type == :elastic
        length(args) == 5 || error("Elastic requires 5 arguments: κ1, μ1, κ2, μ2, f")
        κ1, μ1, κ2, μ2, f = args
        (0.0 <= f <= 1.0) || error("Volume fraction f must be ∈ [0,1]")
        f1, f2 = 1.0 - f, f
        
        # Hashin-Shtrikman bounds
        if bound_type == :hashin_shtrikman
            # The matrix material is phase 1, inclusion is phase 2.
            # Hashin-Shtrikman lower bound is the stiffness of the matrix with spherical
            # inclusions of phase 2.
            κ_lower = κ1 + f2 / (1/(κ2 - κ1) + 3*f1/(3κ1 + 4μ1))
            μ_lower = μ1 + f2 / (1/(μ2 - μ1) + 6*f1*(κ1 + 2μ1)/(5μ1*(3κ1 + 4μ1)))

            # The upper bound is the stiffness of the inclusion material with spherical
            # inclusions of the matrix phase.
            κ_upper = κ2 + f1 / (1/(κ1 - κ2) + 3*f2/(3κ2 + 4μ2))
            μ_upper = μ2 + f1 / (1/(μ1 - μ2) + 6*f2*(κ2 + 2μ2)/(5μ2*(3κ2 + 4μ2)))
            
            # Since voids (E2=0) are always weaker, the order is fixed.
            # For the Voigt and Reuss bounds, the order doesn't matter, so they are unchanged.
            return (κ_lower = κ_lower, κ_upper = κ_upper,
                    μ_lower = μ_lower, μ_upper = μ_upper)
        elseif bound_type == :voigt
             return (κ = f1*κ1 + f2*κ2, μ = f1*μ1 + f2*μ2)
        elseif bound_type == :reuss
            return (
                κ = 1.0 / (f1/κ1 + f2/κ2),
                μ = 1.0 / (f1/μ1 + f2/μ2)
            )
        else
            error("Unsupported bound type for elastic: $bound_type")
        end
    else
        # Thermal case is also fixed
        length(args) == 3 || error("Thermal requires 3 arguments: k1, k2, f")
        k1, k2, f = args
        (0 <= f <= 1) || error("Volume fraction f must be ∈ [0,1]")
        f1, f2 = 1.0 - f, f
        
        if bound_type == :hashin_shtrikman || bound_type == :maxwell
            k_lower = k1 + f2 / (1/(k2 - k1) + f1/(3k1))
            k_upper = k2 + f1 / (1/(k1 - k2) + f2/(3k2))
            return (k_lower = k_lower, k_upper = k_upper)
        elseif bound_type == :voigt
            return (k = f1*k1 + f2*k2)
        elseif bound_type == :reuss
            return (k = 1.0 / (f1/k1 + f2/k2))
        else
            error("Unsupported bound type for thermal: $bound_type")
        end
    end
end