"""
    visualize_periodic_pairs(coords, master_nodes, slave_nodes)

Plots only the matched master-slave nodes and connecting lines.
- Red: slave
- Blue: master
- Green: connection line
"""
function visualize_periodic_pairs(
    coords::Matrix{Float64},
    master_nodes::Vector{Int},
    slave_nodes::Vector{Int})

    dim = size(coords, 2)
    fig = Figure(size=(800, 600))

    if dim == 1
        ax = Axis(fig[1, 1]; title="1D Periodic Pairs")
        for (m, s) in zip(master_nodes, slave_nodes)
            scatter!(ax, [coords[s, 1]], [0.0]; color=:red, markersize=12)
            scatter!(ax, [coords[m, 1]], [0.0]; color=:blue, markersize=12)
            lines!(ax, [coords[s, 1], coords[m, 1]], [0.0, 0.0]; color=:green)
        end
    elseif dim == 2
        ax = Axis(fig[1, 1]; title="2D Periodic Pairs")
        for (m, s) in zip(master_nodes, slave_nodes)
            p1, p2 = coords[s, :], coords[m, :]
            scatter!(ax, [p1[1]], [p1[2]]; color=:red, markersize=12)
            scatter!(ax, [p2[1]], [p2[2]]; color=:blue, markersize=12)
            lines!(ax, [p1[1], p2[1]], [p1[2], p2[2]]; color=:green)
        end
    elseif dim == 3
        ax = Axis3(fig[1, 1]; title="3D Periodic Pairs")
        for (m, s) in zip(master_nodes, slave_nodes)
            p1, p2 = coords[s, :], coords[m, :]
            scatter!(ax, [p1[1]], [p1[2]], [p1[3]]; color=:red, markersize=8)
            scatter!(ax, [p2[1]], [p2[2]], [p2[3]]; color=:blue, markersize=8)
            lines!(ax, [p1[1], p2[1]], [p1[2], p2[2]], [p1[3], p2[3]]; color=:green)
        end
    else
        error("Only supports 1D, 2D, or 3D coordinates.")
    end

    display(fig)
    return fig
end