
# =============================================
# prism integration rule
# =============================================

function integration_rule(element::PrismaticElement, order::Int)
    tri_points, tri_weights = integration_rule(Tri3(),order)  # Pour la base triangulaire
    z, wz = gausslegendre(order)   
    
   n = length(tri_points) * length(z)
    points  = Vector{NTuple{3,Float64}}(undef, n)
    weights = Vector{Float64}(undef, n)
    idx = 1
    for (p_idx, (xi, eta)) in enumerate(tri_points)
        w_tri = tri_weights[p_idx]
        for (zeta, w_z) in zip(z, wz)
            points[idx] = (xi, eta, zeta)
            weights[idx] = w_tri * w_z
            idx += 1
        end
    end
    
    
    return points, weights
end



