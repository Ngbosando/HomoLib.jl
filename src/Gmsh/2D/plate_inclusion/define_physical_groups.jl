"""
    define_physical_groups(
        plate_lines, plate_surface, inclusion_surfaces, inclusion_borders, num_inclusions
    )

    Define physical groups in GMSH model.

    # Arguments
    - "plate_lines": List of plate boundary line tags
    - "plate_surface": Plate surface tag
    - "inclusion_surfaces": List of inclusion surface tags
    - "inclusion_borders": List of inclusion boundary curves
    - "num_inclusions": Number of inclusions

    # Returns
    Tuple containing:
    - Physical group tags for plate boundaries (B_B, B_D, B_H, B_G)
    - List of inclusion frontier tags

    # Notes
    - Creates physical groups for boundaries and domains
    - Assigns meaningful names to groups
    - Returns tags for later node extraction
"""

function define_physical_groups(
    plate_lines, plate_surface, inclusion_surfaces, inclusion_borders, num_inclusions)
    
    B_B = gmsh.model.addPhysicalGroup(1, [plate_lines[1]], 1, "Bottom1")
    B_D = gmsh.model.addPhysicalGroup(1, [plate_lines[2]], 2, "Top1")
    B_H = gmsh.model.addPhysicalGroup(1, [plate_lines[3]], 3, "Right1")
    B_G = gmsh.model.addPhysicalGroup(1, [plate_lines[4]], 4, "Left1")

    gmsh.model.addPhysicalGroup(2, [plate_surface])
    gmsh.model.setPhysicalName(2, 1, "Plate")
    gmsh.model.addPhysicalGroup(2, [inclusion_surfaces...])
    gmsh.model.setPhysicalName(2, 2, "Inclusions")

    inclusion_frontier = Int64[]
    for i in 1:num_inclusions
        BC = gmsh.model.addPhysicalGroup(1, inclusion_borders[i])
        push!(inclusion_frontier, BC)
    end

    return B_B, B_D, B_H, B_G, inclusion_frontier
end
