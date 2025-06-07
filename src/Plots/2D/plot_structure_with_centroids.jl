function plot_structure_with_centroids(
    nodes::Matrix{Float64},
    elements::Matrix{Int},
    elemPhase::Vector{Int})

    fig = Figure(size = (800, 600))
    ax = Axis(fig[1, 1]; aspect = DataAspect(), title = "Mesh Structure")
   
    matrix_tag = findmax(countmap(elemPhase))[2]

    for (i, element) in enumerate(eachrow(elements))
        coords = nodes[element, 1:2]
        x = coords[:, 1]
        y = coords[:, 2]
        x_closed = [x; x[1]]
        y_closed = [y; y[1]]
        lines!(ax, x_closed, y_closed, color = :black)

        # Compute centroid
        cx = mean(x)
        cy = mean(y)

        # Phase color
        color = elemPhase[i] == matrix_tag ? :blue : :red
        scatter!(ax, [cx], [cy], color = color, markersize = 6)
    end

    return fig
end