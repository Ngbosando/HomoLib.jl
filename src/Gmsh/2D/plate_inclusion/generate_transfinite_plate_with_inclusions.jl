include("create_circular_inclusion.jl")
include("create_elliptical_inclusion.jl")
include("create_square_inclusion.jl")
include("create_triangular_inclusion.jl")
include("create_inclusion.jl")
include("define_physical_groups.jl")
include("create_plate.jl")
include("generate_inclusions.jl")
include("generate_mesh_and_retrieve_data.jl")
include("initialize_gmsh.jl")
include("modify_connectivity.jl")
include("assign_PhysicalGroups.jl")
include("sort_periodic_nodes.jl")
include("filter_and_sort_nodes.jl")
include("extract_border_nodes_from_elements.jl")

"""
    generate_transfinite_plate_with_inclusions(
        plate_width,
        plate_height,
        volume_fraction,
        output_file,
        shape,
        N_inclu,
        element_type,
        element_order,
        node_div_inc,
        node_div_mat;
        voids=false,
        show_gui=false,
        rdn=false

    )

    Generate a 2D plate with multiple inclusions using GMSH.

    # Arguments
    - "plate_width": Width of the plate (x-direction)
    - "plate_height": Height of the plate (y-direction)
    - "volume_fraction": Target volume fraction of inclusions
    - "output_file": Base name for output files
    - "shape": Inclusion shape (":circle", ":square", ":ellipse", or ":triangle")
    - "N_inclu": Initial number of inclusions to attempt
    - "element_type": Element type symbol (e.g., :Tri3, :Quad4)
    - "element_order": Element order (1 for linear, 2 for quadratic, etc.)
    - "node_div_inc": Number of nodes along inclusion boundaries
    - "node_div_mat": Number of nodes along plate boundaries

    # Keyword Arguments
    - "voids=false": Whether inclusions should be voids (default: false)
    - "show_gui=false": Whether to display the GMSH GUI (default: false)
    - "rdn=false": Whether inclusions should be randomly distributed (default: false)

    # Returns
    Tuple containing:
    - Boundary node indices (ind_G, ind_D, ind_B, ind_H, ind_C)
    - Element connectivity matrix
    - Node coordinates (nodes_x, nodes_y)
    - Element phase assignments
    - Physical groups information
    - Border nodes dictionary

    # Notes
    - Uses transfinite interpolation for structured meshing
    - Automatically adjusts inclusion sizes to achieve target volume fraction
    - Supports periodic boundary conditions
    - Includes boundary layer refinement around inclusions for quad elements
"""
function generate_transfinite_plate_with_inclusions(
    plate_width,
    plate_height,
    volume_fraction,
    output_file,
    shape,
    N_inclu,
    element_type,
    element_order,
    node_div_inc,
    node_div_mat;
    voids=false,
    show_gui=false,
    rdn = false
)
    # Initialization of Gmsh
    initialize_gmsh()
    gmsh.option.setNumber("General.Terminal", 0)
 
    # Creation of the main plate
    plate_loop, plate_lines = create_plate(plate_width, plate_height)

    # Creation of inclusions
    inclusion_loops, inclusion_surfaces, inclusion_borders, num_inclusions = generate_inclusions(
        volume_fraction,
        plate_width,
        plate_height,
        shape,
        N_inclu,
        voids,
        rdn
    )

    # Creation of the plate surface with inclusions
    plate_surface = gmsh.model.geo.addPlaneSurface( [plate_loop; inclusion_loops...])
    

    gmsh.model.geo.synchronize()

    # Definition of physical groups
    B_B, B_D, B_H, B_G, inclusion_frontier = define_physical_groups(plate_lines,
        plate_surface,
        inclusion_surfaces,
        inclusion_borders,
        num_inclusions
    )

    # Generation of the mesh and retrieval of data
    ind_G, ind_D, ind_B, ind_H, ind_C, elements,
    nodes_x, nodes_y, elemPhase, physGroups, 
    border_nodes = generate_mesh_and_retrieve_data(output_file,
        B_B,
        B_D,
        B_H,
        B_G,
        inclusion_frontier,
        inclusion_surfaces,
        element_order,
        element_type,
        plate_lines,
        inclusion_borders,
        node_div_inc,
        node_div_mat,
        voids
    )

    if show_gui
        gmsh.fltk.run()
    end
    gmsh.finalize()

    return ind_G, ind_D, ind_B, ind_H, ind_C, elements,
    nodes_x, nodes_y, elemPhase, physGroups, border_nodes

end
