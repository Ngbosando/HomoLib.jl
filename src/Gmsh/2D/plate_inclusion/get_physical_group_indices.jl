function get_physical_group_indices(B_G, B_D, B_C)
    ind_G = convert.(Int64, gmsh.model.mesh.getNodesForPhysicalGroup(1, B_G)[1])
    ind_D = convert.(Int64, gmsh.model.mesh.getNodesForPhysicalGroup(1, B_D)[1])
    ind_C = convert.(Int64, gmsh.model.mesh.getNodesForPhysicalGroup(1, B_C)[1])
    return ind_G, ind_D, ind_C
end