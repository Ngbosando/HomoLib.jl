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
    # ========== Input Validation ========== #
    !(field_type in [:elastic, :thermal]) && 
        error("Unsupported field type: $field_type. Use :elastic or :thermal.")
    
    if field_type == :elastic
        # ========== Elastic Case ========== #
        length(args) == 5 || error("Elastic requires 5 arguments: κ1, μ1, κ2, μ2, f")
        κ1, μ1, κ2, μ2, f = args
        (0.0 <= f <= 1.0) || error("Volume fraction f must be ∈ [0,1]")
        f1, f2 = 1.0 - f, f  # f1 = matrix fraction, f2 = inclusion fraction

        # ===== Elastic Bound Types ===== #
        if bound_type == :voigt
            return (κ = f1*κ1 + f2*κ2, μ = f1*μ1 + f2*μ2)

        elseif bound_type == :reuss
            return (
                κ = 1.0 / (f1/κ1 + f2/κ2),
                μ = 1.0 / (f1/μ1 + f2/μ2)
            )

        elseif bound_type == :hashin_shtrikman
            # Order materials so κ1 ≤ κ2
            if κ1 > κ2
                κ1, κ2 = κ2, κ1
                μ1, μ2 = μ2, μ1
            end

            # Bulk modulus bounds
            κ_lower = κ1 + f2 / (1/(κ2 - κ1) + 3*f1/(3κ1 + 4μ1))
            κ_upper = κ2 + f1 / (1/(κ1 - κ2) + 3*f2/(3κ2 + 4μ2))
           
            # Shear modulus bounds
            μ_lower = μ1 + f2 / (1/(μ2 - μ1) + 6*f1*(κ1 + 2μ1)/(5μ1*(3κ1 + 4μ1)))
            μ_upper = μ2 + f1 / (1/(μ1 - μ2) + 6*f2*(κ2 + 2μ2)/(5μ2*(3κ2 + 4μ2)))

            return (
                κ_lower = κ_lower, κ_upper = κ_upper,
                μ_lower = μ_lower, μ_upper = μ_upper
            )

        else
            error("Unsupported bound type for elastic: $bound_type")
        end

    else  # thermal case
        # ========== Thermal Case ========== #
        length(args) == 3 || error("Thermal requires 3 arguments: k1, k2, f")
        k1, k2, f = args
        (0 <= f <= 1) || error("Volume fraction f must be ∈ [0,1]")
        f1, f2 = 1.0 - f, f  # f1 = matrix fraction, f2 = inclusion fraction

        # ===== Thermal Bound Types ===== #
        if bound_type == :voigt
            return (k = f1*k1 + f2*k2)

        elseif bound_type == :reuss
            return (k = 1.0 / (f1/k1 + f2/k2))

        elseif bound_type == :hashin_shtrikman
            # Order materials so k1 ≤ k2
            if k1 > k2
                k1, k2 = k2, k1
            end

            k_lower = k1 + f2 / (1/(k2 - k1) + f1/(3k1))
            k_upper = k2 + f1 / (1/(k1 - k2) + f2/(3k2))

            return (k_lower = k_lower, k_upper = k_upper)

        elseif bound_type == :maxwell
            # Compute both possible configurations
            term_low = (k1 + 2k2 + 2f2*(k1 - k2)) / (k1 + 2k2 - f2*(k1 - k2))
            σ_low = k2 * term_low
            
            term_high = (k2 + 2k1 + 2f1*(k2 - k1)) / (k2 + 2k1 - f1*(k2 - k1))
            σ_high = k1 * term_high
            
            # Ensure correct ordering
            k_lower = min(σ_low, σ_high)
            k_upper = max(σ_low, σ_high)
            return (k_lower = k_lower, k_upper = k_upper)

        else
            error("Unsupported bound type for thermal: $bound_type")
        end
    end
end