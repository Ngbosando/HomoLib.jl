"""
    filter_and_sort_nodes(node_tag, nodes_x, nodes_y, elements)

Filter and sort nodes based on their node_tags and elementsivity.

# Arguments
- `node_tag::Array{Int64}`: Array of node node_tags.
- `nodes_x::Array{Float64}`: Array of x-coordinates of nodes.
- `nodes_y::Array{Float64}`: Array of y-coordinates of nodes.
- `elements::Array{Int64, 2}`: elementsivity matrix.

# Returns
- `Array{Float64, 2}`: Sorted array of nodes with their coordinates.
"""
function filter_and_sort_nodes(node_tag, nodes_x, nodes_y, elements, nodes_per_element)
    u_node_tag = unique(vec(elements))  
    node_tag = convert.(Int64, node_tag)
    node_tag = hcat(node_tag, nodes_x, nodes_y)
    node_tag = filter(row -> row[1] in u_node_tag, eachrow(node_tag))
    node_tag = stack(node_tag)'
    sorted_indices = sortperm(node_tag[:, 1])
    nodes = node_tag[sorted_indices, :]
    return nodes[:, 2:end]
end