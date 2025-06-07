
"""
    initialize_gmsh()

Initialise l'environnement Gmsh et prépare un nouveau modèle pour la génération de maillages.

# Description
Cette fonction configure Gmsh en :
1. Initialisant la bibliothèque Gmsh.
2. Activant l'affichage des messages dans le terminal.
3. Créant un nouveau modèle nommé `"TransfinitePlateWithInclusions"`.

"""

using Gmsh:gmsh
using Random

function initialize_gmsh()
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("TransfinitePlateWithInclusions")
end