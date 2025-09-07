"""
    generate_irregular_quad_patch_mesh(E_o, element_type; nx=4, ny=4, perturb=0.1, filename="quad_patch.msh")

    Generate an irregular quadrilateral patch mesh with random perturbations.

    # Arguments
    - "E_o": Element order (1 for linear, 2 for quadratic, etc.)
    - "element_type": Element type (":triangle" or ":quadrilateral")
    - "nx": Number of elements in x-direction (default: 4)
    - "ny": Number of elements in y-direction (default: 4)
    - "perturb": Perturbation factor for interior nodes (default: 0.1)
    - "filename": Output filename (default: "quad_patch.msh")

    # Returns
    - "nodes": Node coordinates matrix (N×2)
    - "elements": Element connectivity matrix
    - "border_nodes": Dictionary of border node indices

    # Notes
    - Creates a structured but irregular mesh by perturbing interior nodes
    - Uses transfinite interpolation for mesh generation
    - Boundary nodes remain on their original positions
"""

function generate_irregular_quad_patch_mesh(E_o, element_type; nx=4, ny=4, perturb=0.1, filename="quad_patch.msh")
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("irregular_quad_patch")

    # Step size
    dx, dy = 1.0 / nx, 1.0 / ny

    # Point storage
    pts = Matrix{Int}(undef, nx+1, ny+1)
    coords = Matrix{Tuple{Float64, Float64}}(undef, nx+1, ny+1)

    for j in 0:ny, i in 0:nx
        x = i * dx
        y = j * dy
        # Perturb interior 
        if 0 < i < nx && 0 < j < ny
            x += perturb * dx * (2rand() - 1)
            y += perturb * dy * (2rand() - 1)
        end
        coords[i+1, j+1] = (x, y)
        pts[i+1, j+1] = gmsh.model.geo.addPoint(x, y, 0.0)
    end

    # Create lines and transfinite settings
    for j in 1:ny
        for i in 1:nx
            p1 = pts[i, j]
            p2 = pts[i+1, j]
            p3 = pts[i+1, j+1]
            p4 = pts[i, j+1]

            l1 = gmsh.model.geo.addLine(p1, p2)
            l2 = gmsh.model.geo.addLine(p2, p3)
            l3 = gmsh.model.geo.addLine(p3, p4)
            l4 = gmsh.model.geo.addLine(p4, p1)

            cl = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
            sf = gmsh.model.geo.addPlaneSurface([cl])

            gmsh.model.geo.mesh.setTransfiniteCurve(l1, 2)
            gmsh.model.geo.mesh.setTransfiniteCurve(l2, 2)
            gmsh.model.geo.mesh.setTransfiniteCurve(l3, 2)
            gmsh.model.geo.mesh.setTransfiniteCurve(l4, 2)
            gmsh.model.geo.mesh.setTransfiniteSurface(sf)
            gmsh.model.geo.mesh.setRecombine(2, sf)
        end
    end

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)
    # Define the element order 
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

    border_nodes = extract_border_nodes_from_elements(1, box_size=(1,1,0.0))
    

        # gmsh.fltk.run()

    gmsh.finalize()


    return nodes, elements, border_nodes
end