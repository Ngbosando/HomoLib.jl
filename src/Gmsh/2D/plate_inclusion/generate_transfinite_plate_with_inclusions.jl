"""
    generate_transfinite_plate_with_inclusions(plate_width, plate_height, volume_fraction, min_size, max_size, plate_divisions, inclusion_divisions, output_file, shape, n)

Generate a transfinite plate with inclusions.

# Arguments
- `plate_width::Float64`: Width of the plate.
- `plate_height::Float64`: Height of the plate.
- `volume_fraction::Float64`: Some parameter related to the inclusions.
- `min_size::Float64`: Minimum size of the inclusions.
- `max_size::Float64`: Maximum size of the inclusions.
- `plate_divisions::Int`: Number of divisions in the plate.
- `inclusion_divisions::Int`: Number of divisions in the inclusions.
- `output_file::String`: Path to the output file.
- `shape::String`: Shape of the inclusions.
- `n::Int`: Number of inclusions.
- `element_type::String`: Types of an element.
- `element_order::Int`: Order of an element.
- `Nn::Int`: Number of nodes in an element.

# Returns
- `BoundaryNodes::Vector{Int}`: Boundary nodes of the plate.
- `ind_C::Vector{Int}`: Indices of the inclusions.
- `elements::Matrix{Int}`: elementsivity matrix.
- `Nx::Int`: Number of elements in the x direction.
- `Ny::Int`: Number of elements in the y direction.
- `area::Float64`: Area of the plate.
- `elementsElemPhase::Vector{Int}`: elementsivity of the elements and phases.
- `PhysicalGroups::Vector{Int}`: Physical groups in the model.
"""
# include("reshape_element.jl")
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
function generate_transfinite_plate_with_inclusions(
    plate_width,
    plate_height,
    volume_fraction,
    min_size,
    max_size,
    output_file,
    shape,
    N_inclu,
    element_type,
    element_order,
    second_Order_Incomplete,
    node_div_inc,
    node_div_mat
)
    # Initialization of Gmsh
    initialize_gmsh()
    gmsh.option.setNumber("General.Terminal", 0)
 
    # Creation of the main plate
    plate_loop, plate_lines = create_plate(plate_width, plate_height)

    # Creation of inclusions
    inclusion_loops, inclusion_surfaces, inclusion_borders, num_inclusions = generate_inclusions(
        volume_fraction,
        min_size,
        max_size,
        plate_width,
        plate_height,
        shape,
        N_inclu
    )

    # Creation of the plate surface with inclusions
    plate_surface = gmsh.model.geo.addPlaneSurface([plate_loop; inclusion_loops...])

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
        node_div_mat)


    gmsh.fltk.run()
    gmsh.finalize()

    return ind_G, ind_D, ind_B, ind_H, ind_C, elements,
    nodes_x, nodes_y, elemPhase, physGroups, border_nodes

end
