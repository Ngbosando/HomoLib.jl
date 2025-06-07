function cube(filename="cube_hex.msh"; ndiv,b,h,l,element_order,element_type)
    gmsh.initialize()
    gmsh.model.add("cube")
    gmsh.option.setNumber("General.Terminal", 0)
    # 1. Points
    p1 = gmsh.model.geo.addPoint(0, 0, 0)
    p2 = gmsh.model.geo.addPoint( b, 0, 0)
    p3 = gmsh.model.geo.addPoint( b,  h, 0)
    p4 = gmsh.model.geo.addPoint(0,  h, 0)
    p5 = gmsh.model.geo.addPoint(0, 0,  l)
    p6 = gmsh.model.geo.addPoint( b, 0,  l)
    p7 = gmsh.model.geo.addPoint( b,  h,  l)
    p8 = gmsh.model.geo.addPoint(0,  h,  l)

    # 2. Lines
    l1  = gmsh.model.geo.addLine(p1, p2)
    l2  = gmsh.model.geo.addLine(p2, p3)
    l3  = gmsh.model.geo.addLine(p3, p4)
    l4  = gmsh.model.geo.addLine(p4, p1)
    l5  = gmsh.model.geo.addLine(p5, p6)
    l6  = gmsh.model.geo.addLine(p6, p7)
    l7  = gmsh.model.geo.addLine(p7, p8)
    l8  = gmsh.model.geo.addLine(p8, p5)
    l9  = gmsh.model.geo.addLine(p4, p8)
    l10 = gmsh.model.geo.addLine(p5, p1)
    l11 = gmsh.model.geo.addLine(p3, p7)
    l12 = gmsh.model.geo.addLine(p6, p2)
    
    # 3. Curve loops & Plane surfaces
    left = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
    right= gmsh.model.geo.addCurveLoop([l5, l6, l7, l8])
    bottom = gmsh.model.geo.addCurveLoop([l9, l8, l10, -l4])
    top = gmsh.model.geo.addCurveLoop([l2, l11, -l6, l12])
    back = gmsh.model.geo.addCurveLoop([l1, -l12, -l5, l10])
    front= gmsh.model.geo.addCurveLoop([l3, l9, -l7, -l11])

    sf_bottom = gmsh.model.geo.addPlaneSurface([bottom])
    sf_top    = gmsh.model.geo.addPlaneSurface([top])
    sf_front  = gmsh.model.geo.addPlaneSurface([front])
    sf_back   = gmsh.model.geo.addPlaneSurface([back])
    sf_left   = gmsh.model.geo.addPlaneSurface([left])
    sf_right  = gmsh.model.geo.addPlaneSurface([right])

    gmsh.model.geo.synchronize()

    # 4. Surface loop and volume
    surf_loop = gmsh.model.geo.addSurfaceLoop([sf_bottom, sf_top, sf_front, sf_back, sf_left, sf_right])
    vol = gmsh.model.geo.addVolume([surf_loop])

    gmsh.model.geo.synchronize()

    # 5. Physical groups and transfinite
    gmsh.model.addPhysicalGroup(2, [sf_left], 101);   gmsh.model.setPhysicalName(2, 101, "Left")
    gmsh.model.addPhysicalGroup(2, [sf_right], 102);  gmsh.model.setPhysicalName(2, 102, "Right")
    gmsh.model.addPhysicalGroup(2, [sf_bottom], 103); gmsh.model.setPhysicalName(2, 103, "Bottom")
    gmsh.model.addPhysicalGroup(2, [sf_top], 104);    gmsh.model.setPhysicalName(2, 104, "Top")
    gmsh.model.addPhysicalGroup(2, [sf_front], 105);  gmsh.model.setPhysicalName(2, 105, "Front")
    gmsh.model.addPhysicalGroup(2, [sf_back], 106);   gmsh.model.setPhysicalName(2, 106, "Back")
    gmsh.model.addPhysicalGroup(3, [vol], 1);         gmsh.model.setPhysicalName(3, 1, "CubeVolume")

    for l in (l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12)
        gmsh.model.mesh.setTransfiniteCurve(l, ndiv+1)
    end
    for sf in (sf_left, sf_right, sf_bottom, sf_top, sf_front, sf_back)
        gmsh.model.mesh.setTransfiniteSurface(sf)

        etype_str = string(element_type)
        if startswith(etype_str, "Quad") || startswith(etype_str, "Hex")
            
            gmsh.model.mesh.setRecombine(2, sf)
        end
    end
    gmsh.model.mesh.setTransfiniteVolume(vol)
    etype_str = string(element_type)
    if startswith(etype_str, "Quad") || startswith(etype_str, "Hex")

        gmsh.model.mesh.setRecombine(3, vol)
    end
    # gmsh.option.setNumber("Mesh.RecombineAll", 1)
  
    gmsh.model.mesh.generate(3)

    gmsh.model.mesh.setOrder(element_order)


    # Retrieve nodes and elements
    nodes_tuple = gmsh.model.mesh.getNodes()
    nodes_coords = nodes_tuple[2]
    nodes_x = nodes_coords[1:3:end]
    nodes_y = nodes_coords[2:3:end]
    nodes_z = nodes_coords[3:3:end]
    nodes = [nodes_x nodes_y nodes_z]

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
        nodes_ref = NTuple{3,Float64}[]
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
                    push!(nodes_ref, (x, y, z))
                end
            end
        end
        return nodes_ref
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


    ind_L = convert.(Int64, gmsh.model.mesh.getNodesForPhysicalGroup(2, 101)[1])
    ind_R = convert.(Int64, gmsh.model.mesh.getNodesForPhysicalGroup(2, 102)[1])
    ind_B = convert.(Int64, gmsh.model.mesh.getNodesForPhysicalGroup(2, 103)[1])
    ind_T = convert.(Int64, gmsh.model.mesh.getNodesForPhysicalGroup(2, 104)[1])
    ind_F = convert.(Int64, gmsh.model.mesh.getNodesForPhysicalGroup(2, 105)[1])
    ind_Ba = convert.(Int64, gmsh.model.mesh.getNodesForPhysicalGroup(2, 106)[1])
    gmsh.write(filename)
  


    gmsh.fltk.run()
    gmsh.finalize()

    return nodes, connectivity, ind_B, ind_Ba, ind_F, ind_L, ind_T, ind_R
end