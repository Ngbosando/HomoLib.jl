"""
    mori_tanaka_elastic(K_m, G_m, K_i, G_i, v_i)

Compute effective bulk (K) and shear (G) moduli for a two-phase composite using the Mori-Tanaka method.
- `K_m`, `G_m`: Bulk/shear moduli of the matrix.
- `K_i`, `G_i`: Bulk/shear moduli of the inclusion.
- `v_i`: Volume fraction of the inclusion.
"""
# =============================================
# Mori tanaka analytical model
# =============================================

function mori_tanaka_elastic(K_m::Real, G_m::Real, K_i::Real, G_i::Real, v_i::Real)
    v_m = 1.0 - v_i
    
    # Effective bulk modulus
    denominator_K = 1.0 + (v_m * (K_i - K_m)) / (K_m + 4G_m/3)
    K_eff = K_m + v_i*(K_i - K_m) / denominator_K
    
    # Effective shear modulus
    denominator_G = 1.0 + (v_m * (G_i - G_m)) / (G_m + (G_m*(9K_m + 8G_m)) / (6*(K_m + 2G_m)))
    G_eff = G_m + v_i*(G_i - G_m) / denominator_G
    
    return (K=K_eff, G=G_eff)
end



"""
    self_consistent_elastic(K_m, G_m, K_i, G_i, v_i; maxiter=100, tol=1e-6)

Compute effective bulk (K) and shear (G) moduli for a two-phase composite using the Self-Consistent method.
- `K_m`, `G_m`: Bulk/shear moduli of phase 1 (matrix).
- `K_i`, `G_i`: Bulk/shear moduli of phase 2 (inclusion).
- `v_i`: Volume fraction of the inclusion.
- `maxiter`, `tol`: Parameters for numerical convergence.
"""
# =============================================
# self consistent ealstic analytical model
# =============================================

function self_consistent_elastic(K_m::Real, G_m::Real, K_i::Real, G_i::Real, v_i::Real; maxiter=100, tol=1e-6)
    v_m = 1.0 - v_i
    
    # Initial guess (Voigt average)
    K0 = v_m*K_m + v_i*K_i
    G0 = v_m*G_m + v_i*G_i
    
    function equations!(F, x)
        K_sc, G_sc = x
        
        # Bulk modulus equation
        denom_bulk_m = 1.0 + (3*K_sc*(K_m - K_sc)) / (3*K_sc + 4*G_sc)
        term_bulk_m = v_m*(K_m - K_sc) / denom_bulk_m
        
        denom_bulk_i = 1.0 + (3*K_sc*(K_i - K_sc)) / (3*K_sc + 4*G_sc)
        term_bulk_i = v_i*(K_i - K_sc) / denom_bulk_i
        
        F[1] = term_bulk_m + term_bulk_i
        
        # Shear modulus equation
        common_factor = (6*K_sc + 12*G_sc) / (5*(3*K_sc + 4*G_sc))
        
        denom_shear_m = 1.0 + common_factor*(G_m - G_sc)
        term_shear_m = v_m*(G_m - G_sc) / denom_shear_m
        
        denom_shear_i = 1.0 + common_factor*(G_i - G_sc)
        term_shear_i = v_i*(G_i - G_sc) / denom_shear_i
        
        F[2] = term_shear_m + term_shear_i
        
        return F
    end
    
    # Solve numerically
    result = nlsolve(equations!, [K0, G0], autodiff=:forward, iterations=maxiter, ftol=tol)
    
    if converged(result)
        return (K=result.zero[1], G=result.zero[2])
    else
        error("Self-consistent method did not converge")
    end
end

"""
    voigt_elastic(K1, G1, v1, K2, G2, v2)

Compute the Voigt upper bounds for effective bulk (K) and shear (G) moduli (arithmetic mean).
"""
# =============================================
# Voight bound analytical model
# =============================================

function voigt_elastic(K1::Real, G1::Real, v1::Real, K2::Real, G2::Real, v2::Real)
    K_eff = v1*K1 + (1-v1)*K2
    G_eff = v1*G1 + (1-v1)*G2
    return (K=K_eff, G=G_eff)
end

"""
    reuss_elastic(K1, G1, v1, K2, G2, v2)

Compute the Reuss lower bounds for effective bulk (K) and shear (G) moduli (harmonic mean).
"""
# =============================================
# Reuss Bound analytical model
# =============================================

function reuss_elastic(K1::Real, G1::Real, v1::Real, K2::Real, G2::Real, v2::Real)
    K_eff = 1.0 / (v1/K1 + (1-v1)/K2)
    G_eff = 1.0 / (v1/G1 + (1-v1)/G2)
    return (K=K_eff, G=G_eff)
end

"""
    hashin_shtrikman_bounds_elastic(K1, G1, v1, K2, G2, v2)

Compute the Hashin-Shtrikman bounds for bulk (K) and shear (G) moduli of a two-phase composite.
Returns `(K_lower, K_upper, G_lower, G_upper)`.
"""
# =============================================
# Hashin-Shtrikman bounds analytical model
# =============================================

function hashin_shtrikman_bounds_elastic(K1::Real, G1::Real, v1::Real, K2::Real, G2::Real, v2::Real)
    # Ensure K1 <= K2 (swap phases if needed)
    if K1 > K2
        K1, K2 = K2, K1
        G1, G2 = G2, G1
        v1, v2 = v2, v1
    end
    
    # Bulk modulus bounds
    K_lower = K1 + v1 / (1/(K2 - K1) + v2/(K1 + 4*G1/3))
    K_upper = K2 + v2 / (1/(K1 - K2) + v1/(K2 + 4*G2/3))
    
    # Shear modulus bounds
    G_lower = G1 + v1 / (1/(G2 - G1) + (6*v2*(K1 + 2*G1))/(5*G1*(3*K1 + 4*G1)))
    G_upper = G2 + v2 / (1/(G1 - G2) + (6*v2*(K2 + 2*G2))/(5*G2*(3*K2 + 4*G2)))
    
    return (K_lower=K_lower, K_upper=K_upper, G_lower=G_lower, G_upper=G_upper)
end