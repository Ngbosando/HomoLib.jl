# =============================================
# pyramid integration rule
# =============================================

function int_rule(element::PyramidElement, order::Int)
    # Get Gauss-Legendre points and weights for [-1,1] in all 3 dimensions
    r, wr = gausslegendre(order)
    s, ws = gausslegendre(order)
    t_interval, wt = gausslegendre(order)
    
    n = order^3
    points  = Vector{NTuple{3,Float64}}(undef, n)
    weights = Vector{Float64}(undef, n)
    idx = 1
    for i in 1:order, j in 1:order, k in 1:order
        # Cube coordinates
        ξ = r[i]
        η = s[j]
        τ = t_interval[k]
        
        # Map τ from [-1,1] to t ∈ [0,1]
        t = (τ + 1)/2
        
        # Transform to pyramid coordinates
        x = (1 - t) * ξ
        y = (1 - t) * η
        z = t
        
        # Jacobian includes both transformation and dt/dτ
        J = (1 - t)^2 / 2
        
        points[idx] = (x, y, z)
        weights[idx] = wr[i] * ws[j] * wt[k] * J
        idx += 1
    end
    
    return points, weights
end