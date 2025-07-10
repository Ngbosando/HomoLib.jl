
"""
    Transform_boundary(length_Dof, ind_G)

Transforme les indices des degrés de liberté en indices correspondant aux coordonnées `(x, y)` pour chaque nœud de la frontière.

# Arguments
- `length_Dof::Int`: Nombre total de degrés de liberté à transformer.
- `ind_G::Matrix{Int}`: Indices des nœuds appartenant à la frontière (matrice avec une colonne représentant les indices des nœuds).

# Retour
- `Boundary_Dof::Matrix{Int}`: Matrice contenant les indices des degrés de liberté transformés. Chaque ligne contient deux colonnes représentant les indices des coordonnées `(x, y)` pour chaque nœud.

# Description
Cette fonction calcule les indices des degrés de liberté associés aux coordonnées `x` et `y` pour chaque nœud donné dans `ind_G`.  
- Le degré de liberté en `x` pour un nœud est calculé comme `node_index * 2 - 1`.
- Le degré de liberté en `y` pour un nœud est calculé comme `node_index * 2`.

"""

function Transform_boundary(ind)
    length_Dof = length(ind)
    Boundary_Dof = zeros(Int, length_Dof, 2)
    for e in 1:length_Dof
        Boundary_Dof[e, 1] = ind[e, 1] * 2 - 1
        Boundary_Dof[e, 2] = ind[e, 1] * 2
    end

    return Boundary_Dof
end
