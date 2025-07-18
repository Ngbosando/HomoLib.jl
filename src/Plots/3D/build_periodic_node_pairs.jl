"""
    build_periodic_node_pairs(coords::Matrix{Float64}; tol=1e-8)

    Automatically identifies periodic master-slave node pairs along each spatial
    direction by matching nodes on opposite sides of the domain.

    Returns:
    - master_nodes::Vector{Int}
    - slave_nodes::Vector{Int}
"""
function build_periodic_node_pairs(coords::Matrix{Float64}; tol=1e-8)
    dim = size(coords, 2)
    master_nodes = Int[]
    slave_nodes = Int[]

    for d in 1:dim
        min_val = minimum(coords[:, d])
        max_val = maximum(coords[:, d])

        side_min = findall(x -> abs(x[d] - min_val) < tol, eachrow(coords))
        side_max = findall(x -> abs(x[d] - max_val) < tol, eachrow(coords))

        red_axes = setdiff(1:dim, [d])
        coords_min = [coords[i, red_axes] for i in side_min]
        coords_max = [coords[j, red_axes] for j in side_max]

        for (i, c_min) in zip(side_min, coords_min)
            distances = [norm(c_min - c_max) for c_max in coords_max]
            j_idx = argmin(distances)
            j = side_max[j_idx]

            push!(master_nodes, j)
            push!(slave_nodes, i)
        end
    end

    return master_nodes, slave_nodes
end