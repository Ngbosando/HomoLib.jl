"""
    get_num_nodes(element_type::Symbol, element_order::Int) -> Int

    Return the number of nodes for a given element type and order.

    # Arguments
    - "element_type": Element type as a symbol. Supported types are:
    - ":lin" (1D line elements)
    - ":triangle" (2D triangular elements)
    - ":quadrilateral" (2D quadrilateral elements) 
    - ":tetrahedron" (3D tetrahedral elements)
    - ":hexahedron" (3D hexahedral elements)
    - ":prism" (3D prismatic elements)
    - ":pyramid" (3D pyramidal elements)
    - "element_order": Polynomial order of the element (1 to 5)

    # Returns
    Number of nodes for the specified element type and order.
"""
function get_num_nodes(element_type::Symbol, element_order::Int)
    node_table = Dict(
        :lin => [2, 3, 4, 5, 6],
        :triangle => [3, 6, 10, 15, 21],
        :quadrilateral => [4, 9, 16, 25, 36],
        :tetrahedron => [4, 10, 20, 35, 56],
        :hexahedron => [8, 27, 64, 125, 216],
        :prism => [6, 18, 40, 75, 126],
        :pyramid => [5, 13, 30, 55, 91]
    )

    haskey(node_table, element_type) || error("Unsupported element type: $element_type")
    values = node_table[element_type]
    element_order in 1:length(values) || error("Unsupported order $element_order for $element_type")
    return values[element_order]
end
