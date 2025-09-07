"""
    reshape_elements(elements, element_type, element_order)

    Reshape element connectivity array into proper matrix format.

    # Arguments
    - "elements": Flat array of node connectivity
    - "element_type": Element type symbol (e.g., :Tri3, :Quad4)
    - "element_order": Polynomial order of elements

    # Returns
    Matrix with shape (n_elements, nodes_per_element)

    # Notes
    - Validates element family and order
    - Checks for proper array length
    - Uses element configuration dictionary
"""

function reshape_elements(elements, element_type, element_order)
    family = get_element_family(element_type)
    n_nodes = get_node_count(family, element_order)
    check_element_length(elements, n_nodes)

    elements_matrix = reshape(elements, n_nodes, :)
  
    return permutedims(elements_matrix, (2, 1))
end

function get_element_family(element_type)
    family_match = match(r"^[A-Za-z]+", string(element_type))
    isnothing(family_match) && error("Invalid element type: $element_type")
    return Symbol(family_match.match)
end

function get_node_count(family, order)
    config = get_element_config(family)
    haskey(config, order) || error("Unsupported order $order for $family elements")
    return config[order] 
end

function check_element_length(elements, n_nodes)
    if length(elements) % n_nodes != 0
        error("Element list length $(length(elements)) is not divisible by $n_nodes")
    end
end

function get_element_config(family)
    element_configs = Dict(
        :Lin => Dict(1 => 2, 2 => 3, 3 => 4, 4 => 5, 5 => 6),
        :Tri => Dict(1 => 3, 2 => 6, 3 => 10, 4 => 15, 5 => 21),
        :Quad => Dict(1 => 4, 2 => 9, 3 => 16, 4 => 25, 5 => 36),
        :Tet => Dict(1 => 4, 2 => 10, 3 => 20, 4 => 35, 5 => 56),
        :Hex => Dict(1 => 8, 2 => 27, 3 => 64, 4 => 125, 5 => 216),
        :Pri => Dict(1 => 6, 2 => 18, 3 => 40, 4 => 75, 5 => 126),
        :Pyramid => Dict(1 => 5, 2 => 13, 3 => 30, 4 => 55, 5 => 91)
    )
    haskey(element_configs, family) || error("Unsupported element family: $family")
    return element_configs[family]
end
