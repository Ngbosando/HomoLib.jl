function plaque(b, h, lc,lt1,lt2, filename, E_o, element_type::Symbol; show_gui=false)

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add(filename)

    # Create the points of the rectangle
    gmsh.model.geo.addPoint(0,   0, 0, lc, 1)
    gmsh.model.geo.addPoint(b,   0, 0, lc, 2)
    gmsh.model.geo.addPoint(b,    h, 0, lc, 3)
    gmsh.model.geo.addPoint(0,    h, 0, lc, 4)

    # Create the lines of the rectangle:
    # Line 1: from point 1 to 2 (bottom)
    gmsh.model.geo.addLine(1, 2, 1)
    # Line 2: from point 2 to 3 (right)
    gmsh.model.geo.addLine(2, 3, 2)
    # Line 3: from point 3 to 4 (top)
    gmsh.model.geo.addLine(3, 4, 3)
    # Line 4: from point 4 to 1 (left)
    gmsh.model.geo.addLine(4, 1, 4)

    # Create a curve loop and the plane surface
    gmsh.model.geo.addCurveLoop([1, 2, 3, 4], 1)
    gmsh.model.geo.addPlaneSurface([1], 1)

    # Set transfinite meshing options
 
        gmsh.model.geo.mesh.setTransfiniteCurve(1, lt1)
        gmsh.model.geo.mesh.setTransfiniteCurve(2, lt2)
        gmsh.model.geo.mesh.setTransfiniteCurve(3, lt1)
        gmsh.model.geo.mesh.setTransfiniteCurve(4, lt2)
        gmsh.model.geo.mesh.setTransfiniteSurface(1)
    
    if element_type == :quadrilateral
        gmsh.model.geo.mesh.setRecombine(2, 1)
    end

    # -------------------------------
    # Create physical groups for the domain and boundaries
    # -------------------------------

    # Physical group for the domain (the surface, dim = 2)
    gmsh.model.addPhysicalGroup(2, [1], 1)
    gmsh.model.setPhysicalName(2, 1, "Domain")
    
    # Physical groups for each boundary line (dim = 1)
    # Bottom edge (line 1)
    gmsh.model.addPhysicalGroup(1, [1], 2)
    gmsh.model.setPhysicalName(1, 2, "Bottom")
    # Right edge (line 2)
    gmsh.model.addPhysicalGroup(1, [2], 3)
    gmsh.model.setPhysicalName(1, 3, "Right")
    # Top edge (line 3)
    gmsh.model.addPhysicalGroup(1, [3], 4)
    gmsh.model.setPhysicalName(1, 4, "Top")
    # Left edge (line 4)
    gmsh.model.addPhysicalGroup(1, [4], 5)
    gmsh.model.setPhysicalName(1, 5, "Left")

    # Synchronize the geometry before meshing
    gmsh.model.geo.synchronize()

    # Generate the mesh
    gmsh.model.mesh.generate(2)
    gmsh.write(filename * ".msh")

    # Define the element order (set high-order if needed)
    gmsh.model.mesh.setOrder(E_o)

    # Retrieve nodes and elements
    nodes_tuple = gmsh.model.mesh.getNodes()
    nodes_coords = nodes_tuple[2]
    nodes_x = nodes_coords[1:3:end]
    nodes_y = nodes_coords[2:3:end]
    nodes = [nodes_x nodes_y]

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
    if element_type == :triangle
        if E_o == 1
            elem_sym = :Tri3
        elseif E_o == 2
            elem_sym = :Tri6
        elseif E_o == 3
            elem_sym = :Tri10
        elseif E_o == 4
            elem_sym = :Tri15
        elseif E_o == 5
            elem_sym = :Tri21
        end
    elseif element_type == :quadrilateral
        if E_o == 1
            elem_sym = :Quad4   
        elseif E_o == 2
            elem_sym = :Quad9
        elseif E_o == 3
            elem_sym = :Quad16
        elseif E_o == 4
            elem_sym = :Quad25
        elseif E_o == 5
            elem_sym = :Quad36
        end
    end

    elem_code = gmsh_elem_code[elem_sym]

    idx = findfirst(==(elem_code), types)
  
    if isnothing(idx)
        error("Element $elem_sym (code $elem_code) not found in mesh!")
    end
    elements = convert.(Int64, connectivity[3][idx])
    
    elements = reshape_elements(elements,elem_sym,E_o) 

    border_nodes = extract_border_nodes_from_elements(1; box_size=(b,h,0.0),inclusion_borders = nothing)
    
    if show_gui
        gmsh.fltk.run()
    end
    gmsh.finalize()

   
    return nodes, elements, border_nodes
end