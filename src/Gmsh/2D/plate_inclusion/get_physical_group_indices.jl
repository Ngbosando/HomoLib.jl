"""
    get_physical_group_indices(B_G, B_D, B_C) -> (ind_G, ind_D, ind_C)

    Retrieve node indices for physical groups from a GMSH model.

    # Arguments
    - "B_G": Tag for the first physical group (typically left/top boundary)
    - "B_D": Tag for the second physical group (typically right/bottom boundary) 
    - "B_C": Tag for the third physical group (typically inclusion boundaries)

    # Returns
    Tuple containing three arrays of node indices:
    1. "ind_G": Node indices for physical group B_G
    2. "ind_D": Node indices for physical group B_D
    3. "ind_C": Node indices for physical group B_C
"""
function get_physical_group_indices(B_G, B_D, B_C)
    ind_G = convert.(Int64, gmsh.model.mesh.getNodesForPhysicalGroup(1, B_G)[1])
    ind_D = convert.(Int64, gmsh.model.mesh.getNodesForPhysicalGroup(1, B_D)[1])
    ind_C = convert.(Int64, gmsh.model.mesh.getNodesForPhysicalGroup(1, B_C)[1])
    return ind_G, ind_D, ind_C
end
