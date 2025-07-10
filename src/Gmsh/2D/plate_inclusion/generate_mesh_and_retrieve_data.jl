function generate_mesh_and_retrieve_data(output_file,
                                         B_B,
                                         B_D,
                                         B_H,
                                         B_G,
                                         inclusion_frontier,
                                         inclusion_surfaces,
                                         element_order,
                                         element_type,
                                         plate_lines,
                                         inclusion_borders,
                                         node_div_inc,
                                         node_div_mat,
                                         voids;
                                         bl_size=0.05,
                                         bl_size_far=0.20,
                                         bl_ratio=1.3,
                                         bl_thickness=0.4)

    etype_str = string(element_type)
    if startswith(etype_str, "Quad") || startswith(etype_str, "Hex")
        gmsh.option.setNumber("Mesh.RecombineAll", 1)
        gmsh.option.setNumber("Mesh.Algorithm", 11)
        
        # gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 3)  # Blossom full-quad, test further before using this
        # gmsh.option.setNumber("Mesh.Smoothing",            10)  
        # gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 1)  # refine to all-quads, not to do too many nodes
   
        # Boundary layer field around inclusions
    fld = gmsh.model.mesh.field.add("BoundaryLayer")
    curves = reduce(vcat, inclusion_borders)
    gmsh.model.mesh.field.setNumbers(fld, "CurvesList", curves)
    gmsh.model.mesh.field.setNumber( fld, "Size",      bl_size )
    gmsh.model.mesh.field.setNumber( fld, "SizeFar",   bl_size_far )
    gmsh.model.mesh.field.setNumber( fld, "Ratio",     bl_ratio )
    gmsh.model.mesh.field.setNumber( fld, "Thickness", bl_thickness )
    gmsh.model.mesh.field.setNumber( fld, "Quads",     1 )
    gmsh.model.mesh.field.setAsBackgroundMesh(fld)
    else
        gmsh.option.setNumber("Mesh.RecombineAll", 0)
    end
    transfinite = true
    transfiniteAuto = false

    if transfinite
        NNL = node_div_mat
        NNA = node_div_inc
        for line in plate_lines
            gmsh.model.mesh.setTransfiniteCurve(line, NNL)
        end
        for arc in inclusion_borders
            for ac in arc
                gmsh.model.mesh.setTransfiniteCurve(ac, NNA)
            end
        end
    elseif transfiniteAuto
        gmsh.option.setNumber("Mesh.MeshSizeMin", 0.05)
        gmsh.option.setNumber("Mesh.MeshSizeMax", 0.05)
        gmsh.model.mesh.setTransfiniteAutomatic()
    else
        gmsh.option.setNumber("Mesh.MeshSizeMin", 0.05)
        gmsh.option.setNumber("Mesh.MeshSizeMax", 0.05)
    end
 

    # Geometry periodicity setup
    function get_entity_tag_from_physical_name(target_name, target_dim)
        # Fetch all physical groups
        physical_groups = gmsh.model.getPhysicalGroups()
        for (dim, physical_tag) in physical_groups
            # Get the name of the physical group
            name = gmsh.model.getPhysicalName(dim, physical_tag)
            if name == target_name && dim == target_dim
                # Get the entity tags associated with this physical group
                entities = gmsh.model.getEntitiesForPhysicalGroup(dim, physical_tag)
                return entities[1]  # Assumes one entity per physical group
            end
        end
        error("Physical group '$target_name' (dim=$target_dim) not found.")
    end
    left_tag = get_entity_tag_from_physical_name("Left1", 1)
    right_tag = get_entity_tag_from_physical_name("Right1", 1)
    bottom_tag = get_entity_tag_from_physical_name("Bottom1", 1)
    top_tag = get_entity_tag_from_physical_name("Top1", 1)
    translation = [1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1]
  
 
    # Apply the periodic constraint
    # gmsh.model.mesh.setPeriodic(1, [left_tag], [right_tag], translation)
    
    affine = [1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0,
              0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0]
    gmsh.model.mesh.setPeriodic(1, [left_tag], [right_tag], affine)

    gmsh.model.mesh.setPeriodic(1, [bottom_tag], [top_tag], translation)
    # gmsh.model.mesh.setPeriodic(1, [plate_lines[2]], [plate_lines[4]], translation2)  

    # Generate and save mesh
    gmsh.model.mesh.generate(2)
    gmsh.write(output_file)

    # Apply element order and retrieve mesh data
    gmsh.model.mesh.setOrder(element_order)
    tag, raw_nodes = gmsh.model.mesh.getNodes()
    nodes_x = raw_nodes[1:3:end]; nodes_y = raw_nodes[2:3:end]
    nodes = hcat(nodes_x, nodes_y)

    gmsh_elem_code = Dict(
        :Lin2   => 1, :Lin3   => 8,  :Lin4  => 26, :Lin5 => 27, :Lin6 => 28,
        :Tri3   => 2, :Tri6   => 9,  :Tri10 => 21, :Tri15 => 23, :Tri21 => 25,
        :Quad4  => 3, :Quad8  => 16, :Quad9 => 10, :Quad16 => 36, :Quad25 => 37, :Quad36 => 38,
        :Tet4   => 4, :Tet10  => 11, :Tet20 => 29, :Tet35 => 30, :Tet56 => 31,
        :Hex8   => 5, :Hex20  => 17, :Hex27 => 12, :Hex64 => 92, :Hex125 => 93, :Hex216 => 94,
        :Pri6   => 6, :Pri18  => 13, :Pri40 => 90, :Pri75 => 91, :Pri126 => 106,
        :Pyr5   => 7, :Pyr14  => 14, :Pyr30 => 118, :Pyr55 => 119, :Pyr91 => 120
    )
    connectivity = gmsh.model.mesh.getElements()
    types = connectivity[1]
  
  
    elem_code = gmsh_elem_code[element_type]  

    idx = findfirst(==(elem_code), types)
  
    if isnothing(idx)
        error("Element $elem_sym (code $elem_code) not found in mesh!")
    end
    elements = convert.(Int64, connectivity[3][idx])
    
    elements = reshape_elements(elements,element_type,element_order) 
    nodes = filter_and_sort_nodes(tag, nodes[:,1], nodes[:,2], elements, size(elements,2))
    nodes_x, nodes_y = nodes[:,1], nodes[:,2]
    elemPhase, physGroups = assign_physical_groups()
    # elements = reshape_elements(elements,element_type,element_order)



 

    ind_D, ind_G = sort_periodic_nodes(B_D, B_G, 1)
    ind_H, ind_B = sort_periodic_nodes(B_H, B_B, 1; coord_axis=:x)

    ind_BC = Int[]; ind_C = Int[]
    for group in inclusion_frontier
        push!(ind_BC, convert.(Int, gmsh.model.mesh.getNodesForPhysicalGroup(1, group)[1])...)
    end
    for surface in inclusion_surfaces
        push!(ind_C, convert.(Int, gmsh.model.mesh.getNodes(2, surface)[1])...)
    end
    ind_C = vcat(ind_C, ind_BC)
    border_nodes = extract_border_nodes_from_elements(1; box_size=(1.0,1.0,0.0), inclusion_borders)

    # Extract inclusion border elements


    return ind_G, ind_D, ind_B, ind_H, ind_C, elements,
           nodes_x, nodes_y, elemPhase, physGroups, border_nodes
end