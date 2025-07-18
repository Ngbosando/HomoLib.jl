"""
    assign_physical_groups()

    Assign physical groups to elements in GMSH model.

    # Returns
    Tuple containing:
    - "ConnectElemPhase": Physical tag for each element
    - "PhysicalGroups": List of all physical groups

    # Notes
    - Queries GMSH for physical group assignments
    - Returns sorted element-phase pairs
"""

function assign_physical_groups()
    PhysicalGroups = gmsh.model.getPhysicalGroups(2)
    element_phys_pairs = Tuple{Int64, Int64}[]
    for (dim, phys_tag) in PhysicalGroups
        entity_tags = gmsh.model.getEntitiesForPhysicalGroup(dim, phys_tag)
        for entity_tag in entity_tags
            elem_types, elem_tag_arrays = gmsh.model.mesh.getElements(dim, entity_tag)
            for (elem_type, elem_tags) in zip(elem_types, elem_tag_arrays)
                elem_tags = convert.(Int64, elem_tags)
                for elem in elem_tags
                    push!(element_phys_pairs, (elem, phys_tag))
                end
            end
        end
    end
    # Sort the element-physical tag pairs by element ID
    sort!(element_phys_pairs, by = x -> x[1])
    # Extract physical tags in the sorted order
    ConnectElemPhase = [pair[2] for pair in element_phys_pairs]
    return ConnectElemPhase, PhysicalGroups
end


