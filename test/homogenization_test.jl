using LinearAlgebra, SparseArrays, Statistics
using CairoMakie, Test
using HomoLib: generate_transfinite_plate_with_inclusions,
                ThermalMaterial, ElasticMaterial, 
                PiezoelectricMaterial, PoroelasticMaterial,
                assemble_KMF, PeriodicHomogenizationBC, 
                DirichletHomogenizationBC, solve!,
                plot_champ_scalaire, force_computation,
                compute_effective_property,precompute_geometric_data
                
# Mesh generation

mutable struct MeshData
    nodes::Matrix{Float64}
    elements::Matrix{Int}
    type_elem::Vector{Int}
    Nₓ::Vector{Float64}
    Nᵧ::Vector{Float64}
    boundary::Vector{Int}
    boundary_element
    boundary_inclusion::Vector{Int}
    master::Vector{Int}
    slave::Vector{Int}
end

mutable struct ElemData
    type::Symbol
    int_order::Int
    order::Int
end

function setup_mesh(; width, height, volume_fraction,   
                    n_inclusions, Elem::ElemData,
                    node_divisions, shape, output_file, voids, rdn, show_gui, to_rotate)  # Renamed activate -> voids

    ind_G, ind_D, ind_B, ind_H, ind_C,
    elements, Nₓ, Nᵧ, type_elem, _, boundary_element =
        generate_transfinite_plate_with_inclusions(
            width, height, volume_fraction,
            output_file, shape, n_inclusions,
            Elem.type, Elem.order,
            node_divisions[1], node_divisions[2];
            voids=voids,  # Pass voids parameter directly
            show_gui=show_gui,
            rdn=rdn,  # Disable randomization for reproducibility
            to_rotate=to_rotate
        )
  
    boundary_element = boundary_element
    boundary = unique(sort(vcat(ind_G, ind_D, ind_B, ind_H)))
    boundary_inclusion = ind_C
    master   = vcat(ind_G, ind_B)
    slave    = vcat(ind_D, ind_H)
    nodes    = hcat(Nₓ, Nᵧ)

    MeshData(nodes, elements, type_elem, Nₓ, Nᵧ, boundary, boundary_element, boundary_inclusion, master, slave)
end

# Thermal homogenization

function run_thermal_case(mesh::MeshData, Elem::ElemData)
    # Material definition

    mat_thermal_1 = ThermalMaterial(2, k=1.0)
    mat_thermal_2 = ThermalMaterial(2, k=5.0)
 
    dim = 2
    # Precompute data 
    Geometric_Data = precompute_geometric_data(
            Elem.type, Elem.int_order, dim,
            mesh.elements, mesh.nodes, mat_thermal_1
        );
      
        
    # Stiffness matrix
    T_thermal = assemble_KMF(
        mesh.elements,
        mesh.nodes,
        [mat_thermal_1, mat_thermal_2],
        2,
        mesh.type_elem,
        Geometric_Data
    )

    # Periodic boundary conditions
    bc = PeriodicHomogenizationBC
    coor = (mesh.master, mesh.slave)
    bc_x = bc([1.0 0.0], coor, mesh.nodes, 1, [true,true]);
    bc_y = bc([0.0 1.0], coor, mesh.nodes, 1, [true,true]);

    # Solve fields
    U1 = solve!(T_thermal.K, zeros(size(T_thermal.K,1)), zeros(size(T_thermal.K,1)), [bc_x])
    U2 = solve!(T_thermal.K, zeros(size(T_thermal.K,1)), zeros(size(T_thermal.K,1)), [bc_y])

    # Plot temperature fields
    # p1 = plot_champ_scalaire(mesh.nodes, mesh.elements, U1, Elem.type)
    # p2 = plot_champ_scalaire(mesh.nodes, mesh.elements, U2, Elem.type)

    # Recover flux and plot
    # flux_q12 = recover_field_values(
    #     mesh.elements,
    #     [mat_thermal_1, mat_thermal_2],
    #     (T=U1,),
    #     mesh.type_elem,
    #     Elem.type,
    #     Geometric_Data
    # )
    # Q1, Q2 = flux_q12.flux[:,1], flux_q12.flux[:,2]
    # p3 = plot_champ_scalaire(mesh.nodes, mesh.elements, Q1, Elem.type)
    # p4 = plot_champ_scalaire(mesh.nodes, mesh.elements, Q2, Elem.type)

    # Effective thermal conductivity
    κ_eff,_ = compute_effective_property(
        [mat_thermal_1, mat_thermal_2],
        mesh.elements,
        mesh.type_elem,
        (U=(U1, U2),),
        2,
        Geometric_Data
    )
    return κ_eff.K
end


# Elastic homogenization

function run_elastic_case(mesh::MeshData, Elem::ElemData)
    # Material definition
    
    elastic_mat1 = ElasticMaterial(2, E=1.0, ν=0.45, symmetry=:out_of_plane, stress= false)
    elastic_mat2 = ElasticMaterial(2, E=50.0, ν=0.3, symmetry=:out_of_plane, stress= false)

    dim = 2
    # Precompute data 
    Geometric_Data = precompute_geometric_data(
            Elem.type, Elem.int_order, dim,
            mesh.elements, mesh.nodes, elastic_mat1
        );
    # Stiffness matrix
    K_elast = assemble_KMF(
        mesh.elements,
        mesh.nodes,
        [elastic_mat1, elastic_mat2],
        2,
        mesh.type_elem,
        Geometric_Data
    ); 

    # Periodic BCs
    bc = PeriodicHomogenizationBC
    coor = (mesh.master, mesh.slave)
    homogenization_bc  = bc([1.0 0.0; 0.0 0.0], coor,
                                                     mesh.nodes, 2, [true,true]);
    homogenization_bc1 = bc([0.0 0.0; 0.0 1.0], coor,
                                                     mesh.nodes, 2, [true,true]);
    homogenization_bc2 = bc([0.0 1/2; 1/2 0.0], coor,
                                                     mesh.nodes, 2, [true,true]);
    homogenization_bc3 = bc([0.0 0.0; 0.0 0.0], coor,
                                                    mesh.nodes, 2, [true,true]);                                                 

    # Solve mechanical cases
    F_elast = zeros(size(K_elast.K,1))
    U_init  = zeros(size(K_elast.K,1))
    U_mech  = solve!(K_elast.K, F_elast, U_init, [homogenization_bc])
    Umech2  = solve!(K_elast.K, F_elast, U_init, [homogenization_bc1])
    Umech3  = solve!(K_elast.K, F_elast, U_init, [homogenization_bc2])
    Umech4  = solve!(K_elast.K, K_elast.F,      U_init, [homogenization_bc3])

    # Visualization of stress/strain
    # magn   = 0.3
    # Ux     = Umech3[1:2:end]
    # Uy     = Umech3[2:2:end]
    # Nx2    = mesh.nodes[:,1] .+ Ux*magn
    # Ny2    = mesh.nodes[:,2] .+ Uy*magn
    # nodes2 = hcat(Nx2, Ny2)

    # stress_strain = recover_field_values(
    #     mesh.elements,
    #     [elastic_mat1, elastic_mat2],
    #     (U=Umech3,),
    #     mesh.type_elem,
    #     Elem.type,
    #     Geometric_Data
    # )
    # σ₁₁, σ₂₂, σ₁₂ = stress_strain.stress[:,1], stress_strain.stress[:,2], stress_strain.stress[:,3]
    # ϵ₁₁, ϵ₂₂, ϵ₁₂ = stress_strain.strain[:,1], stress_strain.strain[:,2], stress_strain.strain[:,3]
    # plot_champ_scalaire(nodes2, mesh.elements, σ₁₁, Elem.type)
    # plot_champ_scalaire(nodes2, mesh.elements, σ₂₂, Elem.type)
    # plot_champ_scalaire(nodes2, mesh.elements, σ₁₂, Elem.type)
    # plot_champ_scalaire(nodes2, mesh.elements, ϵ₁₁, Elem.type)
    # plot_champ_scalaire(nodes2, mesh.elements, ϵ₂₂, Elem.type)
    # plot_champ_scalaire(nodes2, mesh.elements, ϵ₁₂, Elem.type)

    # Compute effective stiffness
    solver_results = (U=(U_mech, Umech2, Umech3, Umech4),)
    ℂ_eff,_ = compute_effective_property(
        [elastic_mat1, elastic_mat2],
        mesh.elements,
        mesh.type_elem,
        solver_results,
        2,
        Geometric_Data
    )
    return ℂ_eff[1]
end


# Piezoelectric homogenization

function run_piezo_case(mesh::MeshData, Elem::ElemData)
    # Material definition
    C1 = [8.0 4.4 0 4.4; 4.4 8.0 0 4.4; 0 0 3.6 0; 4.4 4.4 0 8]
    C2 = [154.837 83.237 0 82.712; 83.237 154.837 0 82.712;
          0 0 35.8 0; 82.712 82.712 0 131.39]
    e_tensor1 = zeros(3,4)
    e_tensor2 = zeros(3,4); e_tensor2[3,1:2] .= -2.120582; e_tensor2[3,4] = 9.52183
    ϵ1 = diagm(0 => [3.72e-2, 3.72e-2, 3.72e-2])
    ϵ2 = diagm(0 => [4.065, 4.065, 2.079])

    piezo_mat1 = PiezoelectricMaterial(2, C=C1, ε=ϵ1, e=e_tensor1,symmetry=:out_of_plane)
    piezo_mat2 = PiezoelectricMaterial(2, C=C2, ε=ϵ2, e=e_tensor2,symmetry=:out_of_plane)
    dim = 2
    # Precompute data 
    Geometric_Data = precompute_geometric_data(
            Elem.type, Elem.int_order, dim,
            mesh.elements, mesh.nodes, piezo_mat1
        );
    # Assemble coupled stiffness
    K_piezo,_, F = assemble_KMF(
        mesh.elements,
        mesh.nodes,
        [piezo_mat1, piezo_mat2],
        2,
        mesh.type_elem,
        Geometric_Data
    ); 
    F_u, F_ϕ =  F[:strain], F[:electric_field]   
    isapprox(K_piezo, K_piezo', rtol=1e-16) || error("K_piezo is not symmetric")
    
    # Dirichlet BC cases
    bc = PeriodicHomogenizationBC
    coor = (mesh.master, mesh.slave)
    dirichlet_u = [
        [1.0 0.0; 0.0 0.0; 0.0 0.0],
        [0.0 0.0; 0.0 1.0; 0.0 0.0],
        [0.0 0.5; 0.5 0.0; 0.0 0.0],
        [0.0 0.0; 0.0 0.0; 0.0 0.0]
    ]
    dirichlet_ϕ = [
        [0.0 0.0; 0.0 0.0; 1.0 0.0],
        [0.0 0.0; 0.0 0.0; 0.0 1.0],
        [0.0 0.0; 0.0 0.0; 0.0 0.0]
    ]
    bc_u = [bc(val, coor, mesh.nodes, 3, [true,true,true])
            for val in dirichlet_u]
    bc_ϕ = [bc(val, coor, mesh.nodes, 3, [true,true,true])
            for val in dirichlet_ϕ]

    # Solve mechanical and electric loads
    FU = zeros(size(K_piezo,1))
    U1 = solve!(K_piezo, FU, FU, [bc_u[1]])
    U2 = solve!(K_piezo, FU, FU, [bc_u[2]])
    U3 = solve!(K_piezo, FU, FU, [bc_u[3]])
    U4 = solve!(K_piezo, F_u, FU, [bc_u[4]])
    U5 = solve!(K_piezo, FU, FU, [bc_ϕ[1]])
    U6 = solve!(K_piezo, FU, FU, [bc_ϕ[2]])
    U7 = solve!(K_piezo, F_ϕ, FU, [bc_ϕ[3]])

    # Energy verification
    function verify_energy(K, U, F)
        W_int = 0.5 * U' * K * U
        W_ext = 0.5 * U' * F
        @assert abs(W_int - W_ext) < 1e-6*abs(W_ext)
    end
    verify_energy(K_piezo, U4, F_u)
    verify_energy(K_piezo, U7, F_ϕ)

    # Visualization displacement and fields
    # magn   = 0.3
    # Ux     = U1[1:3:end]; Uy = U1[2:3:end]; Up = U1[3:3:end]
    # Nx2    = mesh.nodes[:,1] .+ Ux*magn; Ny2 = mesh.nodes[:,2] .+ Uy*magn
    # nodes2 = hcat(Nx2, Ny2)
    # stress_strain = recover_field_values(
    #     mesh.elements, nodes2,
    #     [piezo_mat1, piezo_mat2],
    #     (U=U1, ϕ=U1),
    #     mesh.type_elem, Elem.type, Elem.int_order, 2,
    #     Geometric_Data
    # )
    # σ₁₁, σ₂₂, σ₁₂ = stress_strain.stress[:,1], stress_strain.stress[:,2], stress_strain.stress[:,3]
    # ϵ₁₁, ϵ₂₂, ϵ₁₂ = stress_strain.strain[:,1], stress_strain.strain[:,2], stress_strain.strain[:,3]
    # E₁, E₂, E₃   = stress_strain.elec[:,1], stress_strain.elec[:,2], stress_strain.elec[:,3]
    # plot_champ_scalaire(nodes2, mesh.elements, E₁, Elem.type)
    # plot_champ_scalaire(nodes2, mesh.elements, E₂, Elem.type)
    # plot_champ_scalaire(nodes2, mesh.elements, E₃, Elem.type)
    # plot_champ_scalaire(nodes2, mesh.elements, ϵ₁₁, Elem.type)
    # plot_champ_scalaire(nodes2, mesh.elements, ϵ₂₂, Elem.type)
    # plot_champ_scalaire(nodes2, mesh.elements, ϵ₁₂, Elem.type)

    # Compute effective piezoelectric properties
    solver_results = (
        U_mech = [U1, U2, U3, U4, U5, U6, U7],
        V_elec = [U1, U2, U3, U4, U5, U6, U7]
    )
    props,_ = compute_effective_property(
        [piezo_mat1, piezo_mat2], mesh.elements,
        mesh.type_elem, solver_results, 2,
        Geometric_Data
    );
    return (C=props.C, e=props.e, ϵ=props.ϵ)
end


# Poroelastic homogenization
 
function run_poro_case(mesh::MeshData, Elem::ElemData)
    # Material definition
    poro_mat = PoroelasticMaterial(2, E=10, ν=1/3,stress = false)
    dim = 2
    # Precompute data no poro
    Geometric_Data = precompute_geometric_data(
            Elem.type, Elem.int_order, dim,
            mesh.elements, mesh.nodes, poro_mat
        );

    # force initialisation
    border_elem = mesh.boundary_element;
    NodalForces = Dict(
    :inclusion => (fᵥ = nothing, fₜ = force_computation((x, y)  -> 1)));
    BoundaryFace = Dict(
    :inclusion => (
        element_border = border_elem[:inclusion_1],
        element_type = :Lin6,
        dim = 2,
        nodes = mesh.nodes,
        int_order = 5
    ));

    # Assemble coupled stiffness
    K_poro, _, F_u = assemble_KMF(
        mesh.elements,
        mesh.nodes,
        poro_mat,
        2,
        Geometric_Data;
        NodalForces,
        BoundaryFace
    );

    isapprox(K_poro, K_poro', rtol=1e-16) || error("K_piezo is not symmetric")

    # Dirichlet BC cases
    dirichlet_u = [
        [1.0 0.0; 0.0 0.0],
        [0.0 0.0; 0.0 1.0],
        [0.0 0.5; 0.5 0.0],
        [0.0 0.0; 0.0 0.0]
    ];
    bc = DirichletHomogenizationBC
    coor = mesh.boundary

    bc_u = [bc(val, coor, mesh.nodes, 2, [true,true])
            for val in dirichlet_u];

    # Solve mechanical loads
    F = zeros(size(K_poro,1));
    U1 = solve!(K_poro, F, F, [bc_u[1]]);
    U2 = solve!(K_poro, F, F, [bc_u[2]]);
    U3 = solve!(K_poro, F, F, [bc_u[3]]);
    U4 = solve!(K_poro, F_u, F, [bc_u[4]]);
  
    # Energy verification
    function verify_energy(K, U, F)
        W_int = 0.5 * U' * K * U
        W_ext = 0.5 * U' * F
        @assert abs(W_int - W_ext) < 1e-10*abs(W_ext)
    end
    verify_energy(K_poro, U4, F_u)
 

    # Visualization displacement and fields
    # magn   = 0.3
    # Ux     = U4[1:2:end]; Uy = U4[2:2:end];
    # Nx2    = mesh.nodes[:,1] .+ Ux*magn; Ny2 = mesh.nodes[:,2] .+ Uy*magn;
    # nodes2 = hcat(Nx2, Ny2);
    # stress_strain = recover_field_values(
    #     mesh.elements,
    #     poro_mat,
    #     (U=U4,),
    #     mesh.type_elem, 
    #     2,
    #     Geometric_Data
    # );
    # σ₁₁, σ₂₂, σ₁₂ = stress_strain.stress[:,1], stress_strain.stress[:,2], stress_strain.stress[:,3];
    # ϵ₁₁, ϵ₂₂, ϵ₁₂ = stress_strain.strain[:,1], stress_strain.strain[:,2], stress_strain.strain[:,3];

    # plot_champ_scalaire(nodes2, mesh.elements, σ₁₁, Elem.type)
    # plot_champ_scalaire(nodes2, mesh.elements, σ₂₂, Elem.type)
    # plot_champ_scalaire(nodes2, mesh.elements, σ₁₂, Elem.type)
    # plot_champ_scalaire(nodes2, mesh.elements, ϵ₁₁, Elem.type)
    # plot_champ_scalaire(nodes2, mesh.elements, ϵ₂₂, Elem.type)
    # plot_champ_scalaire(nodes2, mesh.elements, ϵ₁₂, Elem.type)

    # Compute effective piezoelectric properties
    # Compute effective stiffness
    solver_results = (U=(U1, U2, U3, U4),);
    
    props,volume = compute_effective_property(
        poro_mat,
        mesh.elements,
        solver_results,
        2,
        Geometric_Data
    );
    
    vf = 1-volume
    vt = vf + volume
    Id= [1,1,0]
    f = vf / vt * Id
    fs = volume / vt

    Β_eff = fs * props.Β + f
  
    S = inv(poro_mat.tensors[1])
    i_n = dot(fs * props.Β, S * Id)

    E = 10
    ν = 1/3

    λ_s = E * ν / (1 + ν) / (1 - 2 * ν)
    μ_s = E / (2 * (1 + ν))
    κ_s = λ_s + 2 * μ_s / 3
 
    b = 0.2*(1 + (3*κ_s) /(4*μ_s))
    D = 2
    kr = (κ_s .* (1 .- f)) ./ (1 .+ (f .* D .* κ_s) ./ (2 * (D - 1) * μ_s))
    br = 1 .- kr ./ κ_s
    i_nd = ( 3 * f) / (4 * μ_s)
    i_nr = (br - f) / κ_s
    return (C=props.C, B=Β_eff, b_dilute=b, b_refine=br, solid_biot = i_n,s_biot_dilute =i_nd[1], s_biot_refine=i_nr[1] )
end


# stokes homogenization

#  Elem = ElemData(:Tri3, 1, 1);
#     mesh = setup_mesh(;
#         width=1.0, height=1.0,
#         volume_fraction=0.6, n_inclusions=1,
#         Elem, node_divisions=(2, 2), shape=:circle,
#         output_file="Test_plate_with_inclusions.msh",
#         voids = true,
#         rdn = false,
#         show_gui = true,
#         to_rotate = false
#     );
function run_stokesFlow_case(mesh::MeshData, Elem::ElemData)
    #   place holder
end

# Run Test
# Reference: Yvonnet, Computational Homogenization of Heterogeneous Materials 
# with Finite Elements (2019), Periodic BCs, Fine Mesh

# Thermal Conductivity Tests
@testset "Thermal Conductivity Verification" begin
    # Reference values from Table 3.1 (PER columns for periodic BCs)
    ref_values = Dict(
        1e-3 => 1.000,
        0.05 => 1.0672,
        0.1 => 1.1419,
        0.2 => 1.3055,
        0.3 => 1.4989,
        0.4 => 1.7279,
        0.5 => 2.0104,
        0.6 => 2.3704,
        0.7 => 2.8643
    )
   
    # Test multiple volume fractions
    for (vf, expected_κ) in ref_values

        Elem = ElemData(:Quad36, 10, 5)  # Fine mesh
        mesh = setup_mesh(;
        width=1.0, height=1.0,
        volume_fraction=vf, n_inclusions=1,
        Elem, node_divisions=(3, 3), shape=:circle,
        output_file="Test_plate_with_inclusions_non_poro.msh",
        voids=false,  # Filled inclusions
        rdn=false,  # Disable randomization for reproducibility
        show_gui=false,  # Disable GUI for tests
        to_rotate=false  # Disable rotation for tests
        )
        
        κ_eff = run_thermal_case(mesh, Elem)
   
        @test isapprox(κ_eff[1,1], expected_κ, rtol=0.01)  # 1% tolerance
    end
end

# Elasticity Tests
@testset "Elasticity Verification" begin
    # Reference values from KUBC/PER table for f=0.6
    ref_values = (
        C1111 = (10.961, 11.608),  # (PER, KUBC)
        C1122 = (5.606, 5.495),
        C1212 = (0.955, 2.213),
        C3333 = (34.554, 34.643),
        C1133 = (5.893, 6.047)
    )
    
    Elem = ElemData(:Quad36, 10, 5)  # Fine mesh
    mesh = setup_mesh(;
        width=1.0, height=1.0,
        volume_fraction=0.6, n_inclusions=1,
        Elem, node_divisions=(3, 3), shape=:circle,
        output_file="Test_plate_with_inclusions_non_poro.msh",
        voids=false,  # Filled inclusions
        rdn=false,  # Disable randomization for reproducibility
        show_gui=false,  # Disable GUI for tests
        to_rotate=false  # Disable rotation for tests
    )
    
    C_eff = run_elastic_case(mesh, Elem)
    
    # Test periodic BC results (first value in tuples)
    @test isapprox(C_eff[1,1], ref_values.C1111[1], rtol=0.01)
    @test isapprox(C_eff[1,2], ref_values.C1122[1], rtol=0.01)
    @test isapprox(C_eff[3,3], ref_values.C1212[1], rtol=0.02)  # Higher tolerance for shear
    @test isapprox(C_eff[4,4], ref_values.C3333[1], rtol=0.01)
    @test isapprox(C_eff[1,4], ref_values.C1133[1], rtol=0.01)
end

# Piezoelectricity Tests
@testset "Piezoelectricity Verification" begin
    # Reference values from Table 5.1 (f=0.6, fine mesh Per. b.c.)
    ref_values = (
        C = [25.36 8.71 0 11.92;
             8.71 25.36 0 11.92;
             0 0 0 0;
             11.92 11.92 0 54.62],
        e = [0 0 0 0;
             0 0 0 0;
             -0.20 -0.20 0 6.45],
        ϵ = [0.155 0 0;
             0 0.155 0;
             0 0 1.281]
    )
    
    Elem = ElemData(:Quad36, 10, 5)  # Fine mesh
    mesh = setup_mesh(;
        width=1.0, height=1.0,
        volume_fraction=0.6, n_inclusions=1,
        Elem, node_divisions=(3, 3), shape=:circle,
        output_file="Test_plate_with_inclusions_non_poro.msh",
        voids=false,  # Filled inclusionsts
        rdn=false,  # Disable randomization for reproducibility
        show_gui=false,  # Disable GUI for tests
        to_rotate=false  # Disable rotation for tests
    )

    results = run_piezo_case(mesh, Elem)

    # Test stiffness components
    @test isapprox(results.C[1,1], ref_values.C[1,1], rtol=0.01)
    @test isapprox(results.C[1,2], ref_values.C[1,2], rtol=0.01)
    @test isapprox(results.C[1,3], ref_values.C[1,3], atol=1e-5)
    @test isapprox(results.C[4,4], ref_values.C[4,4], rtol=0.01)

  
    # Test piezoelectric coefficients
    @test isapprox(results.e[3,1], ref_values.e[3,1], rtol=0.05)  # e311
    @test isapprox(results.e[3,3], ref_values.e[3,3], atol=1e-6)  # e312
    @test isapprox(results.e[3,4], ref_values.e[3,4], rtol=0.05)  # e333
    
    # Test dielectric coefficients
    @test isapprox(results.ϵ[1,1], ref_values.ϵ[1,1], rtol=0.02)
    @test isapprox(results.ϵ[3,3], ref_values.ϵ[3,3], rtol=0.01)
end

# Poroelasticity Tests
@testset "Poroelasticity Verification" begin
    # Analytical solution parameters
    E = 10
    ν = 1/3
    α_p = 0.7
    vf = 0.2  # Volume fraction
    
    # Analytical solutions
    λ_s = E * ν / ((1 + ν) * (1 - 2ν))
    μ_s = E / (2 * (1 + ν))
    κ_s = λ_s + 2μ_s/3
    b_analytical = α_p*(1 + (3κ_s)/(4μ_s))
    
    Elem = ElemData(:Quad36, 10, 5)  # Fine mesh
    mesh = setup_mesh(;
        width=1.0, height=1.0,
        volume_fraction=0.2, n_inclusions=1,
        Elem, node_divisions=(10, 5), shape=:circle,
        output_file="Test_plate_with_inclusions_poro.msh",
        voids=true,  # empty inclusions
        rdn=false,  # Disable randomization for reproducibility
        show_gui=false,  # Disable GUI for tests
        to_rotate=false  # Disable rotation for tests
    )
    
    results = run_poro_case(mesh, Elem)

    #  Test Biot’s tensor Compare with analytical solution
    @test isapprox(results.b_refine, results.B, rtol=0.05)
    results.solid_biot
    # Test  solid Biot modulus N with analytical solution
    @test isapprox(results.solid_biot, results.s_biot_refine, rtol=0.07)
  
end



        