export 
# analytical_bound
  theorical_bound
# shape_functions && quadrature
    # shapes 
    shape_functions
    # quadrature
    integration_rule
     
# assembler
    # Jacobian
    compute_jacobian
    # B-matrix
    build_B_matrices
    fill_gradient_B!
    fill_strain_B!
    build_gradient_B
    build_strain_B
    # Precompute_data
    shape_data
    jacobian_data
    # Materials
    material_def
    resolve_material_type
    # Global Assembly
    get_dofs_per_node
    assemble_global_matrix   
    assemble_global_dofs
    # Forces
    # BoundaryFace
    force_computation
    compute_volume_force_vector
    compute_volume_force_scalar
    compute_out_of_plane_force
    compute_surface_force_vector
    compute_surface_force_scalar
    Transform_boundary
# Compute properties 
    compute_effective_property
    init_element_solution
    init_global_storage
    accumulate_results!
    compute_element_contribution
    # compute_flux
    recover_field_values
    integrate_element_quantities
    extract_element_solution_fields
    find_elem_connected_to_nodes
    finalize_field_average
    init_empty_field_sum

# Solver 
    # Solve linear problems  
    BoundaryCondition 
    HomogenizationBC 
    solve!
    # Solve nonlinear problems
    # solve_nonlinear!
    # Solve eigenvalue problems
    # solve_eigenvalue!
    # Solve time-dependent problems
    # solve_time_dependent!
    # Solve optimization problems
    # solve_optimization!

# Mesh
    # Get L-R boundary .
    getBoundaryNodes
    # GMSH
    # Plate with random inclusions
    generate_transfinite_plate_with_inclusions
    create_circular_inclusion
    create_elliptical_inclusion
    create_square_inclusion
    create_triangular_inclusion
    polygons_intersect
    get_rotated_rectangle
    get_rotated_triangle
    get_rotated_ellipse
    create_inclusion
    get_axes
    generate_transfinite_plate_with_inclusions
    define_physical_groups
    create_plate
    generate_inclusions
    generate_mesh_and_retrieve_data
    initialize_gmsh
    reshape_elements
    get_element_family
    get_node_count
    check_element_length
    get_element_config
    assign_physical_groups
    sort_periodic_nodes
    extract_border_nodes_from_elements
    # Simple plate
    plaque
    generate_irregular_quad_patch_mesh
    # cube 
    cube
    build_reference_tetrahedron
# Plots 
    # 2D
    plot_champ_scalaire
    subdivide_element
    subdivide_node
    plot_structure_with_centroids
    # plot_champ_vectoriel

    # 3D
    build_periodic_node_pairs
    visualize_periodic_pairs
    # plot_champ_scalaire_3D
    # plot_champ_vectoriel_3D
