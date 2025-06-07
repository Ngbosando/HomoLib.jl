
# =============================================
# prism integration rule
# =============================================

function integration_rule(element::PrismaticElement, order::Int)
    tri_points, tri_weights = integration_rule(Tri3(),order)  # Pour la base triangulaire
    z, wz = gausslegendre(order)   
    
    # Transformation de [-1,1] à [0,1] 
    # z = 0.5 .* (z .+ 1)
    # wz = 0.5 .* wz

    # Produit tensoriel
    points = [((xi, eta, zeta)) for (xi, eta) in tri_points, zeta in z]
    weights = [w_tri * w_z for w_tri in tri_weights, w_z in wz]
    
    return points, weights
end



