function extract_border_nodes_from_elements(dim::Int; tol=1e-6, box_size=(1.0, 1.0, 1.0))
    elemTypes, elemTags, nodeTags = gmsh.model.mesh.getElements(dim)

    sides = Dict{Symbol, Vector{Vector{Int}}}(
        :left => [], :right => [],
        :front => [], :back => [],
        :bottom => [], :top => []
    )

    Lx, Ly, Lz = box_size
    for (etype, etags, nids) in zip(elemTypes, elemTags, nodeTags)
        n_per_elem = gmsh.model.mesh.getElementProperties(etype)[4]
        le = length(etags)
        for i in 1:le
            nodes = nids[(i-1)*n_per_elem+1 : i*n_per_elem]
            coords = [gmsh.model.mesh.getNode(n)[1] for n in nodes]

            mx = mean(p[1] for p in coords)
            my = mean(p[2] for p in coords)
            mz = mean(p[3] for p in coords)

            if abs(mx - 0.0) < tol
                push!(sides[:left], nodes)
            elseif abs(mx - Lx) < tol
                push!(sides[:right], nodes)
            end

            if abs(my - 0.0) < tol
                push!(sides[:bottom], nodes)
            elseif abs(my - Ly) < tol
                push!(sides[:top], nodes)
            end

            if dim == 2  # only check front/back in 3D
                if abs(mz - 0.0) < tol
                    push!(sides[:front], nodes)
                elseif abs(mz - Lz) < tol
                    push!(sides[:back], nodes)
                end
            end
        end
    end

    return (; sides...)
end