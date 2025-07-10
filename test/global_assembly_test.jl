using Revise 
using HomoLib: compute_effective_property,force_computation
using HomoLib: material_def,assemble_global_matrix,BoundaryCondition,HomogenizationBC,solve!
using HomoLib: plot_champ_scalaire,recover_field_values
using HomoLib: generate_transfinite_plate_with_inclusions,reshape_elements,extract_border_nodes_from_elements
using SparseArrays, Statistics, CairoMakie, Test
using DelaunayTriangulation

# # =============================================
# # Elasticity Test Setup
# # =============================================

    # Define material properties

    elastic_mat = material_def([:elastic], 2, :isotropic, 
                                E=1.0, ν=0.0, plane_stress=true);
                                 e_tensor = zeros(4,3);
   # Example of 2D Piezoelectric Material Definition

# Define material properties (using PZT-5H ceramic as example)
# All values in SI units (GPa for stiffness, C/m² for piezo, pF/m for permittivity)

# Elastic stiffness matrix (3x3 in Voigt notation for 2D plane stress)
C = [126.0  79.5  0;
      79.5 126.0  0;
       0     0   23.0]  # in GPa

# Piezoelectric stress matrix (3x2 for 2D)
e = [
     0   0     -6.5;    # e11, e12, e16
      -17.0    -6.5   0]    # e21, e22, e26  (in C/m²)

# Permittivity matrix (2x2)
ϵ = [25.0  0;    # in nF/m (1505 ε₀)
     0     25.0] 

# Create the material
piezo_mat = material_def([:elastic, :electric], 2, :isotropic,
    C = C * 1e-9 ,    # Convert GPa to Pa
    e = e*1e-6,          # Already in C/m²
    ϵ = ϵ* 1e3 *8.85e-12      # Convert nF/m to F/m
)

# Display the material
println("\n2D Piezoelectric Material Properties:")
show(stdout, piezo_mat)



    # Mesh 
    using HomoLib: getBoundaryNodes,plaque

    # Define mesh parameters
    b = 0.001;
    h = 0.001;
    lc = 0.1; # mesh size
    lt = 1;
    filename = "Test_plate.msh";
    element_order = 2;
    element_type = :triangle;
    nodes, elements, border_nodes= plaque(b, h, lc,lt,lt, filename, element_order, element_type;show_gui=false);

    # =============================================
    # Assembling the Global Stiffness Matrix
    # =============================================

    K_elasticity = assemble_global_matrix(
        elements,
        nodes,
        :Tri6,    # Element type
        2,         # Integration order
        piezo_mat,
        2,          # Dimension
        nothing)
        K = K_elasticity
    Kuu     = K[1:18, 1:18]
    Kuϕ     = K[1:18, 19:27]
    Kϕu     = K[19:27, 1:18]
    Kϕϕ     = K[19:27, 19:27]
    
    println("||Kuu|| = ", maximum(abs.(Kuu)))
    println("||Kuϕ|| = ", maximum(abs.(Kuϕ)))
    println("||Kϕϕ|| = ", maximum(abs.(Kϕϕ)))
    diff_norm = norm(Kϕu - Kuϕ')  # should be close to zero
    println("‖K_ϕuᵀ + K_uϕ‖ = ", diff_norm)
    # =============================================
    # Verification Checks
    # =============================================

    # Basic matrix properties check
    @test issymmetric(K_elasticity) 
    @test all(diag(K_elasticity) .>= 0) 

    # Expected matrix size check
    n_nodes = size(nodes, 1);
    expected_size = n_nodes * 2 ; # 2 DOFs per node (ux, uy)
    @test size(K_elasticity) == (expected_size, expected_size)

    # Translation:


    Δx = 1;
    Δy = 1;
    u = repeat([Δx, Δy], inner=length(nodes[:,1]));

    # Rotation:
    θ = 0.05 ; # small angle in radians
    uθ  = Float64[];
    for (x, y) in eachrow(nodes)
        push!(uθ , -θ*y)
        push!(uθ ,  θ*x)
    end
    uθ  = collect(uθ ) ; # make it an array if needed
                    
    # Apply rigid body motion
    # Check K * rigid_body_mode ≈ 0
    @test norm(K_elasticity * u) < 1e-12
    @test norm(K_elasticity * uθ ) < 1e-12
println("cond(K) = ", cond(Matrix(K)))
   

    f_global = zeros(size(K_elasticity,1))
    f = f_global
    u = f_global
    leftBoundaryNodes, rightBoundaryNodes,top,bot = getBoundaryNodes(nodes, b)
    # Define BCs

    dof_mask = [true, true,true]    # Constrain both ux (1st DOF) and uy (2nd DOF)
    values =  [0.0, 0.0, 0.0]
    dirichlet_bc1 = BoundaryCondition(
        :dirichlet,
        [1],
        dof_mask,
        values,
        3,  # dofs_per_node
    )
    dof_mask = [false, false, true]    # Constrain both ux (1st DOF) and uy (2nd DOF)
    values =  [0.0, 0.0, 100.0]
    dirichlet_bc2 = BoundaryCondition(
        :dirichlet,
        [3],
        dof_mask,
        values,
        3  # dofs_per_node
    )
      dof_mask = [false, false, true]    # Constrain both ux (1st DOF) and uy (2nd DOF)
    values =  [0.0, 0.0, 0.0]
    dirichlet_bc3 = BoundaryCondition(
        :dirichlet,
        [2,4],
        dof_mask,
        values,
        3  # dofs_per_node
    )
    # Enforce periodic BC between:
    # Solve
    U = solve!(K_elasticity, f_global, f_global, [dirichlet_bc1,dirichlet_bc2])

U_vec = copy(U)   # U = solution vector, size (243,)
n_nodes = length(U_vec) ÷ 3

u_x = U_vec[1:3:end]
u_y = U_vec[2:3:end]
ϕ   = U_vec[3:3:end]
maximum(u_y), minimum(u_y), maximum(ϕ)

    N = size(nodes,1); 

    Ux = U[1:3:end];
    Uy = U[2:3:end];
    Nx, Ny = nodes[:, 1], nodes[:, 2];
    Nx₂ = Nx + Ux;
    Ny₂ = Ny + Uy;

    # put Poisson coef = 0 then : 
    # ux1 = [ Nx[i] / b for i in 1:N];
    uxy = -6.5*100*0.01/(8.73e10*0.01)
    errors = [Ux[i] - ux1[i] for i in 1:N]
    L2_error = sqrt(sum(errors .^ 2) / N)
    Linf_error = maximum(abs.(errors))

# # =============================================

# =============================================
# Homogenization Test Setup
# =============================================
    # Define mesh parameters
        plate_width = 1;
        plate_height = 1;
        volume_fraction = 0.6;
        N_inclu = 1;
        min_size = plate_width * sqrt(volume_fraction / (N_inclu * pi));
        max_size = plate_width * sqrt(volume_fraction / (N_inclu * pi));
        output_file = "Test_plate_with_inclusions.msh";
        shape = :circle;
        element_type = :Quad36;
        element_order = 5;
        int_order = element_order*2;
        second_Order_Incomplete = false;
        node_div_inc = 5;
        node_div_mat = 5;

        ind_G, ind_D, ind_B, ind_H, ind_C, elements, Nₓ, Nᵧ,
        type_elem, physical_group,
        border_nodes = generate_transfinite_plate_with_inclusions(
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
        );

        BoundaryNodes = vcat(ind_G, ind_D, ind_B, ind_H);
        BoundaryNodes = sort(BoundaryNodes);
        BoundaryNodes = unique(BoundaryNodes);
        nodes = hcat(Nₓ, Nᵧ);
        master_nodes = vcat(ind_G,ind_B);
        slave_nodes = vcat(ind_D,ind_H);

    # thermal case
        # Material definition
            mat_thermal_1 = material_def([:thermal], 2, :isotropic; κ=1);
            # # Incluions
            mat_thermal_2 = material_def([:thermal], 2, :isotropic; κ=5);

        # Stiffness matrix
            T_thermal = assemble_global_matrix(
                elements,
                nodes,
                element_type,    # Element type
                int_order,         # Integration order
                [mat_thermal_1, mat_thermal_2],
                2,       # Dimension
                type_elem;)


        # Boudary and solving
            dirichlet_bc = HomogenizationBC(
                :dirichlet,
                [1.0 0.0],
                nothing,
                BoundaryNodes,
                [true, true],  # Assuming 1 DOF per node
                nodes,
                1
            );
            dirichlet_bc1 = HomogenizationBC(
                :dirichlet,
                [0.0 1.0],
                nothing,
                BoundaryNodes,
                [true, true],  # Assuming 1 DOF per node
                nodes,
                1
            );
            homogenization_bc = HomogenizationBC(
                :periodic,
                [1.0 0.0],
                (master_nodes,slave_nodes),
                nothing,
                [true,true],
                nodes,
                1
            );
            homogenization_bc1 = HomogenizationBC(
                :periodic,
                [0.0 1.0],
                (master_nodes,slave_nodes),
                nothing,
                [true,true],
                nodes,
                1
            );

            U₁ = solve!(T_thermal,zeros(size(T_thermal,1)),zeros(size(T_thermal,1)), [dirichlet_bc]);
            U₂ = solve!(T_thermal, zeros(size(T_thermal,1)), zeros(size(T_thermal,1)), [dirichlet_bc1]);


        # Plots temperature field and flux
            p = plot_champ_scalaire(nodes, elements, U₁, element_type)
            p2 = plot_champ_scalaire(nodes, elements, U₂, element_type)

            flux_q12 = recover_field_values(
                elements,
                nodes,
                [mat_thermal_1,mat_thermal_2],  # Prop now contains full conductivity matrices
                (T = U₁,) ,
                type_elem,
                element_type,
                int_order,
                2
            );
            Q1,Q2 = flux_q12.flux[:,1], flux_q12.flux[:,2];
            p3 = plot_champ_scalaire(nodes, elements, Q1, element_type)
            p4 = plot_champ_scalaire(nodes, elements, Q2, element_type)
            # Patch test k_m = k_i
                # Nn = length(Nₓ);
                # TH = [ Nₓ[i] / 1 for i in 1:Nn];
                # # N = length(T);
                # errors = [U₁[i] - Uᵣ₁[i] for i in 1:Nn];
                # # Normes
                # η_error = sqrt(sum(errors .^ 2) / Nn)
                # Linf_error = maximum(abs.(errors))
                
        # Effective properties computation
            solver_results = (
                U = (U₁, U₂),  # Matrix of size (n_nodes × 2)
            );
            
            κ_eff = compute_effective_property(
                [mat_thermal_1, mat_thermal_2],
                elements,
                nodes,
                type_elem,
                solver_results,  
                element_type,
                int_order,
            2
            )

    # Elastic case
        # Comparaisan with Matlab data
        # using MAT
        #     using HomoLib: getBoundaryNodes
        #     data     = matread("mesh_data.mat")
        #     nodes    = data["nodes"]    # Array{Float64,2} de taille N×2
        #     elements = data["elements"] # Array{Int64,2} de taille Nel×3
        #     elements = Int.(elements)
        #     type_elem = data["type_elem"] # Array{Symbol,1} de taille Nel
        #     type_elem = vec(Int.(type_elem))

        #     leftBoundaryNodes,rightBoundaryNodes,topBoundaryNodes,bottomBoundaryNodes = getBoundaryNodes(nodes,1)
        #     BoundaryNodes = vcat(leftBoundaryNodes, rightBoundaryNodes, topBoundaryNodes, bottomBoundaryNodes);
        #     BoundaryNodes = sort(BoundaryNodes);
        #     BoundaryNodes = unique(BoundaryNodes);
            # using HomoLib: build_periodic_node_pairs,visualize_periodic_pairs,plot_structure_with_centroids
            # using StatsBase
            # master_nodes,slave_nodes =build_periodic_node_pairs(nodes);
            # visualize_periodic_pairs(nodes, master_nodes, slave_nodes)
            # fig = plot_structure_with_centroids(nodes, elements, type_elem)
            # display(plot_structure_with_centroids(nodes, elements, fill(1, size(elements, 1))))

        # Material definition 
            elastic_mat1= material_def([:elastic], 2, :out_of_plane,  
                                    E=1.0, ν=0.45, plane_stress=false)
            elastic_mat2 = material_def([:elastic], 2, :out_of_plane, 
                                    E=50.0, ν=0.3, plane_stress=false)
        # Force definition in case of out_of_plane definition or Subc 
            NodalForces= (fᵥ = [0.0 0.0],fₜ = nothing); 
        # Stiffness matrix
            K_elast,F = assemble_global_matrix(
                        elements,
                        nodes,
                        element_type,    # Element type
                        int_order,         # Integration order
                        [elastic_mat1, elastic_mat2],
                        2,          # Dimension
                        type_elem;
                       
                        );

        # Boudary definition and solving    
            dirichlet_bc = HomogenizationBC(
                :dirichlet,
                [1.0 0.0 0.0; 0.0 0.0 0.0 ; 0.0 0.0 0.0],
                nothing,
                BoundaryNodes,
                [true, true],  # Assuming 1 DOF per node
                nodes,
                2
            );
            dirichlet_bc1 = HomogenizationBC(
                :dirichlet,
                [0.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 0.0],
                nothing,
                BoundaryNodes,
                [true, true],  # Assuming 1 DOF per node
                nodes,
                2
            );

            dirichlet_bc2 = HomogenizationBC(
                :dirichlet,
                [0.0 1/2 0.0; 1/2 0.0 0.0; 0.0 0.0 0.0],
                nothing,
                BoundaryNodes,
                [true, true],  # Assuming 1 DOF per node
                nodes,
                2
            );
            dirichlet_bc3 = HomogenizationBC(
                :dirichlet,
                [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 1.0],
                nothing,
                BoundaryNodes,
                [true, true],  # Assuming 1 DOF per node
                nodes,
                2
            );

            # periodic boundary
                homogenization_bc = HomogenizationBC(
                    :periodic,
                    [1.0 0.0; 0.0 0.0],
                    (master_nodes,slave_nodes),
                    nothing,
                    [true,true],
                    nodes,
                    2
                );
                homogenization_bc1 = HomogenizationBC(
                    :periodic,
                    [0.0 0.0; 0.0 1.0],
                    (master_nodes,slave_nodes),
                    nothing,
                    [true,true],
                    nodes,
                    2
                );
                homogenization_bc2 = HomogenizationBC(
                    :periodic,
                    [0.0 1/2; 1/2 0.0] ,
                    (master_nodes,slave_nodes),
                    nothing,
                    [true,true],
                    nodes,
                    2
                );
                homogenization_bc3 = HomogenizationBC(
                    :periodic,
                    [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 1.0],
                    (master_nodes,slave_nodes),
                    nothing,
                    [true,true],
                    nodes,
                    2
                );

            F_elast = zeros(size(K_elast,1));
            U_elast = zeros(size(K_elast,1));

            U_mech = solve!(K_elast, F_elast, U_elast, [dirichlet_bc]);
            Umech2 = solve!(K_elast, F_elast, U_elast, [dirichlet_bc1]);
            Umech3 = solve!(K_elast, F_elast, U_elast, [dirichlet_bc2]);
            Umech4 = solve!(K_elast, F, zeros(size(K_elast,1)), [dirichlet_bc3]);

        # Visualisation displacement and strain/stress
            magn=0.3;
            Ux = Umech3[1:2:end];
            Uy = Umech3[2:2:end];
            Nx₂ = nodes[:,1] + Ux * magn;
            Ny₂ = nodes[:,2] + Uy * magn ;
            nodes2 = hcat(Nx₂, Ny₂);
            # Displacement
                # points = Point2f.(nodes[:,1], nodes[:,2]);
                # On fait une triangulation Delaunay sur ces points
                # tri = triangulate(points);
                # fig, ax, _ = triplot(
                # tri;                   # Triangulation existante
                # show_points    = false,
                # triangle_color = :transparent,
                # strokecolor    = :blue,
                # strokewidth    = 0.8
                # );
                # n_major = 5;
                # x_min, x_max = minimum(nodes[:,1]), maximum(nodes[:,1]);
                # y_min, y_max = minimum(nodes[:,2]), maximum(nodes[:,2]);
                # xt = range(x_min, x_max; length = n_major);
                # yt = range(y_min, y_max; length = n_major);

                # # labels arrondis à 1 décimale
                # xtl = string.(round.(xt; digits = 1));
                # ytl = string.(round.(yt; digits = 1));

                # ax.xticks = (xt, xtl);
                # ax.yticks = (yt, ytl);
                # fig

            # Stress _strain
                stress_strain = recover_field_values(
                    elements,
                    nodes,
                    [elastic_mat1,elastic_mat2],  # Prop now contains full conductivity matrices
                    (U = Umech3,) ,
                    type_elem,
                    element_type,
                    int_order,
                    2
                );

                σ₁₁, σ₂₂, σ₁₂ = stress_strain.stress[:,1], stress_strain.stress[:,2], stress_strain.stress[:,3];
                ϵ₁₁, ϵ₂₂, ϵ₁₂ = stress_strain.strain[:,1], stress_strain.strain[:,2], stress_strain.strain[:,3];
                σ₁₁_plot = plot_champ_scalaire( nodes2, elements, σ₁₁, element_type)
                σ₂₂_plot = plot_champ_scalaire( nodes2, elements, σ₂₂, element_type)
                σ₁₂_plot = plot_champ_scalaire( nodes2, elements, σ₁₂, element_type)
                ϵ₁₁_plot = plot_champ_scalaire( nodes2, elements, ϵ₁₁, element_type)
                ϵ₂₂_plot = plot_champ_scalaire( nodes2, elements, ϵ₂₂, element_type)
                ϵ₁₂_plot = plot_champ_scalaire( nodes2, elements, ϵ₁₂, element_type)

        # Compute effective properties
            solver_results = (
                U = (U_mech, Umech2, Umech3,Umech4),  # Matrix of size (n_nodes × 6)
            );

            ℂ_eff = compute_effective_property(
                [elastic_mat1, elastic_mat2],
                elements,
                nodes,
                type_elem,
                solver_results,  # Dimension
                element_type,
                2,
            2
            );
            ℂ_eff = ℂ_eff[1]

    # In case of K_eff non-symettry to impose it
# =============================================
# Elasticity patch Test Setup
# =============================================
