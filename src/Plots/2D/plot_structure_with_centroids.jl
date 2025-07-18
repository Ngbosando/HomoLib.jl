"""
    plot_structure_with_centroids(nodes, elements, elemPhase)

    Visualize finite element mesh with element centroids colored by phase.

    # Arguments
    - "nodes::Matrix{Float64}": Node coordinates (n_nodes × 2)
    - "elements::Matrix{Int}": Element connectivity (n_elements × nodes_per_element)
    - "elemPhase::Vector{Int}": Phase/material assignment for each element

    # Returns
    A Makie.jl Figure object displaying:
    - Mesh edges in black
    - Element centroids colored by phase:
    - Matrix phase (most common phase) in blue
    - Inclusion phases in red

    # Usage Example
    ```julia
    nodes = rand(100, 2)  # 100 nodes in 2D
    elements = rand(1:100, 50, 3)  # 50 triangular elements
    phases = rand([1,2], 50)  # Random phase assignment

    fig = plot_structure_with_centroids(nodes, elements, phases)
    display(fig)

    Example Use Cases
    Verifying mesh quality
    Checking phase distribution in composites
    Debugging material assignments

    Visualizing RVE structures

    See Also
    plot_champ_scalaire: For field visualization on meshes
    visualize_periodic_pairs: For periodic boundary visualization
"""

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