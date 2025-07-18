 """
    Sort nodes on periodic boundaries (slave and master) based on spatial coordinates.
    Returns `(slave_nodes_sorted, master_nodes_sorted)`.

    Args:
        slave_phys_tag: Physical group tag of the slave boundary (e.g., "Left1").
        master_phys_tag: Physical group tag of the master boundary (e.g., "Right1").
        dim: Dimension of the entities (1 for edges, 2 for surfaces).
        coord_axis: Axis to sort by (`:y` for vertical edges, `:x` for horizontal edges).
    """

function sort_periodic_nodes(slave_phys_tag::Int32, master_phys_tag::Int32, dim::Int64; coord_axis=:y)
   
    # Get all nodes on slave and master boundaries
    slave_nodes = gmsh.model.mesh.getNodesForPhysicalGroup(dim, slave_phys_tag)[1]
    master_nodes = gmsh.model.mesh.getNodesForPhysicalGroup(dim, master_phys_tag)[1]

    # Check if boundaries have the same number of nodes
    @assert length(slave_nodes) == length(master_nodes) "Periodic boundaries must have the same number of nodes."

    # Extract coordinates for slave nodes
    slave_coords = [gmsh.model.mesh.getNode(tag)[1] for tag in slave_nodes]
    slave_keys = coord_axis == :y ? [c[2] for c in slave_coords] : [c[1] for c in slave_coords]
    slave_sorted = slave_nodes[sortperm(slave_keys)]  # Sort by y or x

    # Extract coordinates for master nodes
    master_coords = [gmsh.model.mesh.getNode(tag)[1] for tag in master_nodes]
    master_keys = coord_axis == :y ? [c[2] for c in master_coords] : [c[1] for c in master_coords]
    master_sorted = master_nodes[sortperm(master_keys)]  # Sort by y or x

    # Optional: Reverse master order if edges are oriented oppositely
    # Example: Left edge sorted bottom-to-top vs. right edge top-to-bottom
    if coord_axis == :y
        first_slave_y = slave_coords[argmin(slave_keys)][2]
        first_master_y = master_coords[argmin(master_keys)][2]
        if !isapprox(first_slave_y, first_master_y; atol=1e-6)
            master_sorted = reverse(master_sorted)
        end
    end
    slave_sorted = convert.(Int64, slave_sorted)
    master_sorted = convert.(Int64, master_sorted)

    return slave_sorted, master_sorted
end

