 
# =============================================
# Line integration rule
# =============================================
function integration_rule(element::LinearElement, order::Int) 
    ξ, w = gausslegendre(order)
    return ξ, w  # 1D points stored as column vector
end

# =============================================
# quadrilateral integration rule
# =============================================

function integration_rule(element::QuadrilateralElement, order::Int) 
    # Obtenir les points et les poids pour une intégration sur l'intervalle [-1, 1]
    x, wx = gausslegendre(order)
    y, wy = gausslegendre(order)
    
    # Construire le produit tensoriel pour obtenir les points et poids 2D
    points = [(xi, yi) for xi in x, yi in y]
    weights = [wxi * wyi for wxi in wx, wyi in wy]
    
    return points, weights
end

# =============================================
# Hexaherael integration rule
# =============================================

function integration_rule(element::HexahedralElement, order::Int)
    # Points et poids 1D pour l'intervalle [-1,1]
    x, wx = gausslegendre(order)
    y, wy = gausslegendre(order)
    z, wz = gausslegendre(order)
    
    # Produit tensoriel pour obtenir les points et poids en 3D
    points = [(xi, yi, zi) for xi in x, yi in y, zi in z]
    weights = [wxi * wyi * wzi for wxi in wx, wyi in wy, wzi in wz]
    
    return points, weights
end
