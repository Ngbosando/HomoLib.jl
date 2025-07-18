include("cube.jl")

"""
    build_reference_tetrahedron(order::Int, element_type; show_gui=false)

    Generate a reference tetrahedron mesh for finite element analysis.

    # Arguments
    - "order": Polynomial order of the elements
    - "element_type": Symbol specifying element type (e.g., :Tet4, :Tet10)
    - "show_gui": Whether to display the GMSH GUI (default: false)

    # Returns
    - "nodes": Node coordinates (N×3 matrix)
    - "connectivity1": First connectivity array
    - "connectivity2": Second connectivity array

    # Notes
    - Creates a reference tetrahedron with vertices at (0,0,0), (1,0,0), (0,1,0), (0,0,1)
    - Includes additional mesh information and properties
    - Can compute basis functions and Jacobians at integration points
"""

function build_reference_tetrahedron(order::Int,element_type); show_gui=false
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("ReferenceTetrahedron")

    # Define 4 points of the reference tetrahedron
    p1 = gmsh.model.occ.addPoint(0.0, 0.0, 0.0,2)
    p2 = gmsh.model.occ.addPoint(1.0, 0.0, 0.0,2)
    p3 = gmsh.model.occ.addPoint(0.0, 1.0, 0.0,2)
    p4 = gmsh.model.occ.addPoint(0.0, 0.0, 1.0,2)

    # Define lines
    l1 = gmsh.model.occ.addLine(p1, p2)
    l2 = gmsh.model.occ.addLine(p2, p3)
    l3 = gmsh.model.occ.addLine(p3, p1)

    l4 = gmsh.model.occ.addLine(p1, p4)
    l5 = gmsh.model.occ.addLine(p2, p4)
    l6 = gmsh.model.occ.addLine(p3, p4)

    # Define triangular surfaces
    loop1 = gmsh.model.occ.addCurveLoop([l1, l2, -l3])
    s1 = gmsh.model.occ.addPlaneSurface([loop1])  # Base: p1-p2-p3

    loop2 = gmsh.model.occ.addCurveLoop([l1, l5, -l4])
    s2 = gmsh.model.occ.addPlaneSurface([loop2])  # Side: p1-p2-p4

    loop3 = gmsh.model.occ.addCurveLoop([l2, l6, -l5])
    s3 = gmsh.model.occ.addPlaneSurface([loop3])  # Side: p2-p3-p4

    loop4 = gmsh.model.occ.addCurveLoop([l3, l4, -l6])
    s4 = gmsh.model.occ.addPlaneSurface([loop4])  # Side: p3-p1-p4

    # Create the volume from surfaces
    sl = gmsh.model.occ.addSurfaceLoop([s1, s2, s3, s4])
    vol = gmsh.model.occ.addVolume([sl])
    gmsh.model.occ.synchronize()
  
    # Set element order and generate mesh
    gmsh.option.setNumber("Mesh.ElementOrder", order)
    gmsh.model.mesh.generate(3)



     gmsh_elem_code = Dict(
        :Lin2   => 1, :Lin3   => 8,  :Lin4  => 26, :Lin5 => 27, :Lin6 => 28,
        :Tri3   => 2, :Tri6   => 9,  :Tri10 => 21, :Tri15 => 23, :Tri21 => 25,
        :Quad4  => 3, :Quad8  => 16, :Quad9 => 10, :Quad16 => 36, :Quad25 => 37, :Quad36 => 38,
        :Tet4   => 4, :Tet10  => 11, :Tet20 => 29, :Tet35 => 30, :Tet56 => 31,
        :Hex8   => 5, :Hex20  => 17, :Hex27 => 12, :Hex64 => 92, :Hex125 => 93, :Hex216 => 94,
        :Pri6   => 6, :Pri18  => 13, :Pri40 => 90, :Pri75 => 91, :Pri126 => 106,
        :Pyr5   => 7, :Pyr14  => 14, :Pyr30 => 118, :Pyr55 => 119, :Pyr91 => 120
    )

    element_configs = Dict(
        :Lin => Dict(1 => 2, 2 => 3, 3 => 4, 4 => 5, 5 => 6),
        :Tri => Dict(1 => 3, 2 => 6, 3 => 10, 4 => 15, 5 => 21),
        :Quad => Dict(1 => 4, 2 => 9, 3 => 16, 4 => 25, 5 => 36),
        :Tet => Dict(1 => 4, 2 => 10, 3 => 20, 4 => 35, 5 => 56),
        :Hex => Dict(1 => 8, 2 => 27, 3 => 64, 4 => 125, 5 => 216),
        :Pri => Dict(1 => 6, 2 => 18, 3 => 40, 4 => 75, 5 => 126),
        :Pyramid => Dict(1 => 5, 2 => 13, 3 => 30, 4 => 55, 5 => 91)
    )
      nodes_tuple = gmsh.model.mesh.getNodes()
    nodes_coords = nodes_tuple[2]
    nodes_x = nodes_coords[1:3:end]
    nodes_y = nodes_coords[2:3:end]
    nodes_z = nodes_coords[3:3:end]
    nodes = [nodes_x nodes_y nodes_z]

   elements = gmsh.model.mesh.getElements()
    types = elements[1]

  
    elem_code = gmsh_elem_code[element_type]  

    idx = findfirst(==(elem_code), types)

    if isnothing(idx)
        error("Element $elem_sym (code $elem_code) not found in mesh!")
    end
    connectivity = convert.(Int64, elements[3][idx])
 
    
    function reshape_elements(connectivity::Vector{Int}, elem_sym::Symbol, order::Int)
        # Find matching family
        family = nothing
        for fam in keys(element_configs)
            if startswith(string(elem_sym), string(fam))
                family = fam
                break
            end
        end
        family === nothing && error("Unknown element symbol: $elem_sym")
        # Number of nodes per element
        nodes_per_elem = element_configs[family][order]
        n_elem = div(length(connectivity), nodes_per_elem)
        return reshape(connectivity, nodes_per_elem, n_elem)'
    end

    function generate_tet_tensor_nodes_cartesian(p::Int)
        nodes = NTuple{3,Float64}[]
        for a in 0:p
            for b in 0:(p - a)
                for c in 0:(p - a - b)
                    d = p - a - b - c
                    λ1 = a / p
                    λ2 = b / p
                    λ3 = c / p
                    λ4 = d / p
                    # Reference tetrahedron vertices: (0,0,0), (1,0,0), (0,1,0), (0,0,1)
                    x = λ2
                    y = λ3
                    z = λ4
                    push!(nodes, (x, y, z))
                end
            end
        end
        return nodes
    end


    function get_gmsh_tensor_permutation_from_session(element_family::Symbol, order::Int)
        element_type = gmsh.model.mesh.getElementType(string(element_family), order)
        _, _, _, num_nodes, node_locs, _ = gmsh.model.mesh.getElementProperties(element_type)

        gmsh_coords = [node_locs[3i-2:3i] for i in 1:num_nodes]
      

        if element_family == :Hexahedron
            n = order + 1
            ref_pts = [(x, y, z) for z in range(-1, 1, length=n),
                                y in range(-1, 1, length=n),
                                x in range(-1, 1, length=n)]
        elseif element_family == :Tetrahedron
            ref_pts = generate_tet_tensor_nodes_cartesian(order)
       
        else
            error("Unsupported element family: $element_family")
        end

        perm = fill(-1, num_nodes)
        used = falses(num_nodes)

        for (i, ref) in enumerate(ref_pts)
            best_idx = 0
            best_dist = Inf
            for j in 1:num_nodes
                if used[j]; continue; end
                dist = norm(gmsh_coords[j] .- collect(ref))
                if dist < best_dist
                    best_dist = dist
                    best_idx = j
                end
            end
            @assert best_dist < 1e-6 "Failed to match tensor point $ref to any Gmsh node"
            perm[i] = best_idx
            used[best_idx] = true
        end

        return perm
    end

    connectivity = Matrix(reshape_elements(connectivity,element_type,element_order))
    function reorder_connectivity_to_tensor(connect::Matrix{Int}, element_family::Symbol, order::Int)
        perm = get_gmsh_tensor_permutation_from_session(element_family, order)
        return connect[:, perm]
    end
    # Infer 3D element family from symbol
    element_families_3D = Dict(
        :Hex => :Hexahedron,
        :Tet => :Tetrahedron,
        :Pri => :Prism,
        :Pyr => :Pyramid
    )

    family = nothing
    for fam in keys(element_families_3D)
        if startswith(string(element_type), string(fam))
            family = fam
            break
        end
    end

    if family !== nothing
        gmsh_family = element_families_3D[family]
        connectivity = reorder_connectivity_to_tensor(connectivity, gmsh_family, element_order)
    end

   
    interpolationOrder = order-1


    function pp(label, v, mult)
        println(" * " * string(length(v) / mult) * " " * label * ": " * string(v))
    end

    # Iterate over all the element types present in the mesh:
    elementTypes = gmsh.model.mesh.getElementTypes()

    for t in elementTypes
        # Retrieve properties for the given element type
        elementName, dim, order, numNodes, localNodeCoord, numPrimNodes =
            gmsh.model.mesh.getElementProperties(t)
        println("\n** " * elementName * " **\n")

        # Retrieve integration points for that element type, enabling exact
        # integration of polynomials of order "interpolationOrder". The "Gauss"
        # integration family returns the "economical" Gauss points if available, and
        # defaults to the "CompositeGauss" (tensor product) rule if not.
        localCoords, weights =
            gmsh.model.mesh.getIntegrationPoints(t, "Gauss" * string(interpolationOrder))
        pp("integration points to integrate order " *
        string(interpolationOrder) * " polynomials", localCoords, 3)

        # Return the basis functions evaluated at the integration points. Selecting
        # "Lagrange" and "GradLagrange" returns the isoparamtric basis functions and
        # their gradient (in the reference space of the given element type). A
        # specific interpolation order can be requested using "LagrangeN" and
        # "GradLagrangeN" with N = 1, 2, ... Other supported function spaces include
        # "H1LegendreN", "GradH1LegendreN", "HcurlLegendreN", "CurlHcurlLegendreN".
        numComponents, basisFunctions, numOrientations =
            gmsh.model.mesh.getBasisFunctions(t, localCoords, "Lagrange")
        pp("basis functions at integration points", basisFunctions, 1)
        numComponents, basisFunctions, numOrientations =
            gmsh.model.mesh.getBasisFunctions(t, localCoords, "GradLagrange")
        pp("basis function gradients at integration points", basisFunctions, 3)

        # Compute the Jacobians (and their determinants) at the integration points
        # for all the elements of the given type in the mesh. Beware that the
        # Jacobians are returned "by column": see the API documentation for details.
        jacobians, determinants, coords =
            gmsh.model.mesh.getJacobians(t, localCoords)
        pp("Jacobian determinants at integration points", determinants, 1)
    end

    if show_gui
        gmsh.fltk.run()
    end
    gmsh.finalize()
    return nodes, connectivity[1], connectivity[2]
end