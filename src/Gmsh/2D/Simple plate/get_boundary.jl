"""
    getBoundaryNodes(nodeCoordinates, L)

    Identify boundary nodes of a square domain.

    # Arguments
    - "nodeCoordinates": N×2 matrix of node coordinates
    - "L": Length of the square domain

    # Returns
    Tuple containing four arrays:
    - "leftBoundaryNodes": Nodes on x=0 boundary
    - "rightBoundaryNodes": Nodes on x=L boundary
    - "topBoundaryNodes": Nodes on y=L boundary
    - "bottomBoundaryNodes": Nodes on y=0 boundary

    # Notes
    - Uses a tolerance of 1e-12 to identify boundary nodes
    - Works for both structured and unstructured meshes
"""

function getBoundaryNodes(nodeCoordinates, L)
    # Initialisation des vecteurs pour stocker les indices des nœuds sur les frontières
    leftBoundaryNodes = Int[]  # Nœuds sur la frontière gauche (x ≈ 0)
    rightBoundaryNodes = Int[] # Nœuds sur la frontière droite (x ≈ L)
    topBoundaryNodes = Int[]   # Nœuds sur la frontière supérieure (y ≈ L)
    bottomBoundaryNodes = Int[] # Nœuds sur la frontière inférieure (y ≈ 0)

    N = size(nodeCoordinates, 1);
    # Boucle sur tous les nœuds
    for i in 1: N
        if nodeCoordinates[i, 1] < 1e-12   # Frontière gauche (x ≈ 0)
            push!(leftBoundaryNodes, i)    # Ajouter l'indice du nœud à la frontière gauche
        elseif nodeCoordinates[i, 1] > L - 1e-12  # Frontière droite (x ≈ L)
            push!(rightBoundaryNodes, i)   # Ajouter l'indice du nœud à la frontière droite
        end
    end

    for i in 1: N
        if nodeCoordinates[i, 2] < 1e-12   # Frontière inférieure (y ≈ 0)
            push!(bottomBoundaryNodes, i)  # Ajouter l'indice du nœud à la frontière inférieure
        elseif nodeCoordinates[i, 2] > L - 1e-12  # Frontière supérieure (y ≈ L)
            push!(topBoundaryNodes, i)     # Ajouter l'indice du nœud à la frontière supérieure
        end
    end
   

    return leftBoundaryNodes, rightBoundaryNodes, topBoundaryNodes, bottomBoundaryNodes
end