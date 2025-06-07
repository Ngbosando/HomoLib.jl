function plaque(b, h, lc,lt, filename, E_o, element_type::Symbol)

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add(filename)

    # Create the points of the rectangle
    gmsh.model.geo.addPoint(0,   0, 0, lc, 1)
    gmsh.model.geo.addPoint(b,   0, 0, lc, 2)
    gmsh.model.geo.addPoint(b,   h, 0, lc, 3)
    gmsh.model.geo.addPoint(0,   h, 0, lc, 4)

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
    if lt != 0
        gmsh.model.geo.mesh.setTransfiniteCurve(1, lt)
        gmsh.model.geo.mesh.setTransfiniteCurve(2, lt)
        gmsh.model.geo.mesh.setTransfiniteCurve(3, lt)
        gmsh.model.geo.mesh.setTransfiniteCurve(4, lt)
        gmsh.model.geo.mesh.setTransfiniteSurface(1)
    end
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

    elements = gmsh.model.mesh.getElements()
    connectivity = convert.(Int64, elements[3][2])
    type = elements[1]
 println("type = $type")
    if element_type == :triangle
        if E_o == 1
            connectivity = [connectivity[1:3:end] connectivity[2:3:end] connectivity[3:3:end]]
        elseif E_o == 2
            connectivity = [connectivity[1:6:end] connectivity[2:6:end] connectivity[3:6:end] connectivity[4:6:end] connectivity[5:6:end] connectivity[6:6:end]]
        elseif E_o == 3
            connectivity = [connectivity[1:10:end] connectivity[2:10:end] connectivity[3:10:end] connectivity[4:10:end] connectivity[5:10:end] connectivity[6:10:end] connectivity[7:10:end] connectivity[8:10:end] connectivity[9:10:end] connectivity[10:10:end]]
        elseif E_o == 4
            connectivity = [connectivity[1:15:end] connectivity[2:15:end] connectivity[3:15:end] connectivity[4:15:end] connectivity[5:15:end] connectivity[6:15:end] connectivity[7:15:end] connectivity[8:15:end] connectivity[9:15:end] connectivity[10:15:end] connectivity[11:15:end] connectivity[12:15:end] connectivity[13:15:end] connectivity[14:15:end] connectivity[15:15:end]]
        elseif E_o == 5
            connectivity = [connectivity[1:21:end] connectivity[2:21:end] connectivity[3:21:end] connectivity[4:21:end] connectivity[5:21:end] connectivity[6:21:end] connectivity[7:21:end] connectivity[8:21:end] connectivity[9:21:end] connectivity[10:21:end] connectivity[11:21:end] connectivity[12:21:end] connectivity[13:21:end] connectivity[14:21:end] connectivity[15:21:end] connectivity[16:21:end] connectivity[17:21:end] connectivity[18:21:end] connectivity[19:21:end] connectivity[20:21:end] connectivity[21:21:end]]
        end
    elseif element_type == :quadrilateral
        if E_o == 1
            connectivity = [connectivity[1:4:end] connectivity[2:4:end] connectivity[3:4:end] connectivity[4:4:end]]
        elseif E_o == 2
            connectivity = [connectivity[1:9:end] connectivity[2:9:end] connectivity[3:9:end] connectivity[4:9:end] connectivity[5:9:end] connectivity[6:9:end] connectivity[7:9:end] connectivity[8:9:end] connectivity[9:9:end]]
        elseif E_o == 3
            connectivity = [connectivity[1:16:end] connectivity[2:16:end] connectivity[3:16:end] connectivity[4:16:end] connectivity[5:16:end] connectivity[6:16:end] connectivity[7:16:end] connectivity[8:16:end] connectivity[9:16:end] connectivity[10:16:end] connectivity[11:16:end] connectivity[12:16:end] connectivity[13:16:end] connectivity[14:16:end] connectivity[15:16:end] connectivity[16:16:end]]
        elseif E_o == 4
            connectivity = [connectivity[1:25:end] connectivity[2:25:end] connectivity[3:25:end] connectivity[4:25:end] connectivity[5:25:end] connectivity[6:25:end] connectivity[7:25:end] connectivity[8:25:end] connectivity[9:25:end] connectivity[10:25:end] connectivity[11:25:end] connectivity[12:25:end] connectivity[13:25:end] connectivity[14:25:end] connectivity[15:25:end] connectivity[16:25:end] connectivity[17:25:end] connectivity[18:25:end] connectivity[19:25:end] connectivity[20:25:end] connectivity[21:25:end] connectivity[22:25:end] connectivity[23:25:end] connectivity[24:25:end] connectivity[25:25:end]]
        elseif E_o == 5
            connectivity = [connectivity[1:36:end] connectivity[2:36:end] connectivity[3:36:end] connectivity[4:36:end] connectivity[5:36:end] connectivity[6:36:end] connectivity[7:36:end] connectivity[8:36:end] connectivity[9:36:end] connectivity[10:36:end] connectivity[11:36:end] connectivity[12:36:end] connectivity[13:36:end] connectivity[14:36:end] connectivity[15:36:end] connectivity[16:36:end] connectivity[17:36:end] connectivity[18:36:end] connectivity[19:36:end] connectivity[20:36:end] connectivity[21:36:end] connectivity[22:36:end] connectivity[23:36:end] connectivity[24:36:end] connectivity[25:36:end] connectivity[26:36:end] connectivity[27:36:end] connectivity[28:36:end] connectivity[29:36:end] connectivity[30:36:end] connectivity[31:36:end] connectivity[32:36:end] connectivity[33:36:end] connectivity[34:36:end] connectivity[35:36:end] connectivity[36:36:end]]
        end
    end

    
    # gmsh.fltk.run()
    gmsh.finalize()

    # Return nodes, connectivity, and physical groups (you can filter these if desired)
    return nodes, connectivity
end