function extract_border_nodes_from_elements(dim::Int; tol=1e-6, box_size=(1.0,1.0,1.0), inclusion_borders)
    # Prepare dict: side name => Vector of element connectivity vectors
    sides = Dict{Symbol, Matrix{Int}}()
    # Plate boundaries: collect element connectivity whose centroids lie on borders
    elemTypes, elemTags, nodeTags = gmsh.model.mesh.getElements(dim)
    Lx, Ly, Lz = box_size
    for side in (:left, :right, :bottom, :top)
        sides[side] = Array{Int,2}(undef, 0, 0)
    end
    # Identify typical connectivity size by first element
    if !isempty(elemTypes)
        n0 = gmsh.model.mesh.getElementProperties(elemTypes[1])[4]
        # initialize with zero rows, n0 columns
        for side in keys(sides)
            sides[side] = Array{Int,2}(undef, 0, n0)
        end
    end
    # Loop through elements
    for (etype, etags, nids) in zip(elemTypes, elemTags, nodeTags)
        n_per = gmsh.model.mesh.getElementProperties(etype)[4]
        for (i, tag) in enumerate(etags)
            conn = nids[(i-1)*n_per+1 : i*n_per]
            # centroid
            coords = [gmsh.model.mesh.getNode(n)[1] for n in conn]
            mx = mean(p[1] for p in coords)
            my = mean(p[2] for p in coords)
            # test each plate side
            if abs(mx - 0.0) < tol
                sides[:left] = vcat(sides[:left], reshape(conn, 1, n_per))
            elseif abs(mx - Lx) < tol
                sides[:right] = vcat(sides[:right], reshape(conn, 1, n_per))
            end
            if abs(my - 0.0) < tol
                sides[:bottom] = vcat(sides[:bottom], reshape(conn, 1, n_per))
            elseif abs(my - Ly) < tol
                sides[:top] = vcat(sides[:top], reshape(conn, 1, n_per))
            end
        end
    end
    # Inclusion borders: each arc defines curve elements (1D) connectivity
    for (i, arcs) in enumerate(inclusion_borders)
        key = Symbol("inclusion_$(i)")
        # determine element size for curves: first arc
        if !isempty(arcs)
            # assume all curves same etype
            etypes, etags, nodeTags1 = gmsh.model.mesh.getElements(1, arcs[1])
            n_per_curve = gmsh.model.mesh.getElementProperties(etypes[1])[4]
            sides[key] = Array{Int,2}(undef, 0, n_per_curve)
        else
            sides[key] = zeros(Int, 0, 0)
        end
        for arc in arcs
            etypes, etags, nodeTags1 = gmsh.model.mesh.getElements(1, arc)
            for (etype, etags2, nids2) in zip(etypes, etags, nodeTags1)
                n_per = gmsh.model.mesh.getElementProperties(etype)[4]
                for j in 1:length(etags2)
                    conn = nids2[(j-1)*n_per+1 : j*n_per]
                    sides[key] = vcat(sides[key], reshape(conn, 1, n_per))
                end
            end
        end
    end
    return sides
end
