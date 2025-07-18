using Gmsh:gmsh
using Random

"""
    initialize_gmsh()

    Initialize GMSH and create a new model.

    # Notes
    - Must be called before any other GMSH operations
    - Sets terminal output level to 1 (basic messages)
    - Creates a default model named "TransfinitePlateWithInclusions"
"""
function initialize_gmsh()
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("TransfinitePlateWithInclusions")
end