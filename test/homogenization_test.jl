using Revise
using LinearAlgebra, SparseArrays, Statistics
using CairoMakie, Test
using HomoLib: generate_transfinite_plate_with_inclusions,
                  material_def, assemble_global_matrix,
                  HomogenizationBC, BoundaryCondition, solve!,
                  plot_champ_scalaire, recover_field_values,
                  compute_effective_property,force_computation

# =============================================
# Main workflow
# =============================================
function main()
    Elem = ElemData(:Quad36, 5*2, 5)
    mesh_non_poro = setup_mesh(;
        width=1.0, height=1.0,
        volume_fraction=0.6, n_inclusions=1,
        Elem, node_divisions=(3, 3), shape=:circle,
        output_file="Test_plate_with_inclusions_non_poro.msh",
        voids=false  # Filled inclusions
    )
  
    # Mesh for poro case (void inclusions)
    mesh_poro = setup_mesh(;
        width=1.0, height=1.0,
        volume_fraction=0.2, n_inclusions=1,
        Elem, node_divisions=(10, 5), shape=:circle,
        output_file="Test_plate_with_inclusions_poro.msh",
        voids=true  # Void inclusions
    )

    # println("thermal start")
    # κ_eff     = run_thermal_case(mesh_non_poro, Elem)
    # @info "Thermal effective conductivity:" κ_eff

    # println("elastic start")
    # C_eff     = run_elastic_case(mesh_non_poro, Elem)
    # @info "Elastic effective stiffness tensor:" C_eff

    # println("piezo start")
    # piezo_eff = run_piezo_case(mesh_non_poro, Elem)
    # @info "Piezoelectric effective properties:" piezo_eff

    println("poro start")
    Biot_tensor = run_poro_case(mesh_poro, Elem)
    @info "porileasticity effective properties:" Biot_tensor
end

# =============================================
# Mesh generation
# =============================================
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
                    node_divisions, shape, output_file, voids)  # Renamed activate -> voids

    ind_G, ind_D, ind_B, ind_H, ind_C,
    elements, Nₓ, Nᵧ, type_elem, _, boundary_element =
        generate_transfinite_plate_with_inclusions(
            width, height, volume_fraction,
            output_file, shape, n_inclusions,
            Elem.type, Elem.order,
            node_divisions[1], node_divisions[2];
            voids=voids,  # Pass voids parameter directly
            show_gui=false
        )
  
    boundary_element = boundary_element
    boundary = unique(sort(vcat(ind_G, ind_D, ind_B, ind_H)))
    boundary_inclusion = ind_C
    master   = vcat(ind_G, ind_B)
    slave    = vcat(ind_D, ind_H)
    nodes    = hcat(Nₓ, Nᵧ)

    MeshData(nodes, elements, type_elem, Nₓ, Nᵧ, boundary, boundary_element, boundary_inclusion, master, slave)
end

# =============================================
# Thermal homogenization
# =============================================
function run_thermal_case(mesh::MeshData, Elem::ElemData)
    # Material definition
    mat_thermal_1 = material_def([:thermal], 2, :isotropic; κ=1)
    mat_thermal_2 = material_def([:thermal], 2, :isotropic; κ=5)

    # Stiffness matrix
    T_thermal = assemble_global_matrix(
        mesh.elements,
        mesh.nodes,
        Elem.type,
        Elem.int_order,
        [mat_thermal_1, mat_thermal_2],
        2,
        mesh.type_elem
    )

    # Periodic boundary conditions
    bc_x = HomogenizationBC(
        :periodic,
        [1.0 0.0],
        (mesh.master, mesh.slave),
        nothing,
        nothing,
        mesh.nodes,
        1
    )
    bc_y = HomogenizationBC(
        :periodic,
        [0.0 1.0],
        (mesh.master, mesh.slave),
        nothing,
        nothing,
        mesh.nodes,
        1
    )

    # Solve fields
    U1 = solve!(T_thermal, zeros(size(T_thermal,1)), zeros(size(T_thermal,1)), [bc_x])
    U2 = solve!(T_thermal, zeros(size(T_thermal,1)), zeros(size(T_thermal,1)), [bc_y])

    # Plot temperature fields
    p1 = plot_champ_scalaire(mesh.nodes, mesh.elements, U1, Elem.type)
    p2 = plot_champ_scalaire(mesh.nodes, mesh.elements, U2, Elem.type)

    # Recover flux and plot
    flux_q12 = recover_field_values(
        mesh.elements,
        mesh.nodes,
        [mat_thermal_1, mat_thermal_2],
        (T=U1,),
        mesh.type_elem,
        Elem.type,
        Elem.int_order,
        2
    )
    Q1, Q2 = flux_q12.flux[:,1], flux_q12.flux[:,2]
    p3 = plot_champ_scalaire(mesh.nodes, mesh.elements, Q1, Elem.type)
    p4 = plot_champ_scalaire(mesh.nodes, mesh.elements, Q2, Elem.type)

    # Effective thermal conductivity
    κ_eff,_ = compute_effective_property(
        [mat_thermal_1, mat_thermal_2],
        mesh.elements,
        mesh.nodes,
        mesh.type_elem,
        (U=(U1, U2),),
        Elem.type,
        Elem.int_order,
        2
    )
    return κ_eff
end

# =============================================
# Elastic homogenization
# =============================================
function run_elastic_case(mesh::MeshData, Elem::ElemData)
    # Material definition
    elastic_mat1 = material_def([:elastic], 2, :out_of_plane, E=1.0, ν=0.45)
    elastic_mat2 = material_def([:elastic], 2, :out_of_plane, E=50.0, ν=0.3)

    # Stiffness matrix
    K_elast, F = assemble_global_matrix(
        mesh.elements,
        mesh.nodes,
        Elem.type,
        Elem.int_order,
        [elastic_mat1, elastic_mat2],
        2,
        mesh.type_elem
    )

    # Periodic BCs
    homogenization_bc  = HomogenizationBC(:periodic, [1.0 0.0; 0.0 0.0],
                                          (mesh.master, mesh.slave), nothing,
                                          [true,true], mesh.nodes, 2)
    homogenization_bc1 = HomogenizationBC(:periodic, [0.0 0.0; 0.0 1.0],
                                          (mesh.master, mesh.slave), nothing,
                                          [true,true], mesh.nodes, 2)
    homogenization_bc2 = HomogenizationBC(:periodic, [0.0 1/2; 1/2 0.0],
                                          (mesh.master, mesh.slave), nothing,
                                          [true,true], mesh.nodes, 2)
    homogenization_bc3 = HomogenizationBC(:periodic, [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 1.0],
                                          (mesh.master, mesh.slave), nothing,
                                          [true,true], mesh.nodes, 2)

    # Solve mechanical cases
    F_elast = zeros(size(K_elast,1))
    U_init  = zeros(size(K_elast,1))
    U_mech  = solve!(K_elast, F_elast, U_init, [homogenization_bc])
    Umech2  = solve!(K_elast, F_elast, U_init, [homogenization_bc1])
    Umech3  = solve!(K_elast, F_elast, U_init, [homogenization_bc2])
    Umech4  = solve!(K_elast, F,      U_init, [homogenization_bc3])

    # Visualization of stress/strain
    magn   = 0.3
    Ux     = Umech3[1:2:end]
    Uy     = Umech3[2:2:end]
    Nx2    = mesh.nodes[:,1] .+ Ux*magn
    Ny2    = mesh.nodes[:,2] .+ Uy*magn
    nodes2 = hcat(Nx2, Ny2)

    stress_strain = recover_field_values(
        mesh.elements,
        mesh.nodes,
        [elastic_mat1, elastic_mat2],
        (U=Umech3,),
        mesh.type_elem,
        Elem.type,
        Elem.int_order,
        2
    )
    σ₁₁, σ₂₂, σ₁₂ = stress_strain.stress[:,1], stress_strain.stress[:,2], stress_strain.stress[:,3]
    ϵ₁₁, ϵ₂₂, ϵ₁₂ = stress_strain.strain[:,1], stress_strain.strain[:,2], stress_strain.strain[:,3]
    plot_champ_scalaire(nodes2, mesh.elements, σ₁₁, Elem.type)
    plot_champ_scalaire(nodes2, mesh.elements, σ₂₂, Elem.type)
    plot_champ_scalaire(nodes2, mesh.elements, σ₁₂, Elem.type)
    plot_champ_scalaire(nodes2, mesh.elements, ϵ₁₁, Elem.type)
    plot_champ_scalaire(nodes2, mesh.elements, ϵ₂₂, Elem.type)
    plot_champ_scalaire(nodes2, mesh.elements, ϵ₁₂, Elem.type)

    # Compute effective stiffness
    solver_results = (U=(U_mech, Umech2, Umech3, Umech4),)
    ℂ_eff,_ = compute_effective_property(
        [elastic_mat1, elastic_mat2],
        mesh.elements,
        mesh.nodes,
        mesh.type_elem,
        solver_results,
        Elem.type,
        Elem.int_order,
        2
    )
    return ℂ_eff[1]
end

# =============================================
# Piezoelectric homogenization
# =============================================
function run_piezo_case(mesh::MeshData, Elem::ElemData)
    # Material definition
    C1 = [8.0 4.4 0 4.4; 4.4 8.0 0 4.4; 0 0 3.6 0; 4.4 4.4 0 8]
    C2 = [154.837 83.237 0 82.712; 83.237 154.837 0 82.712;
          0 0 35.8 0; 82.712 82.712 0 131.39]
    e_tensor1 = zeros(3,4)
    e_tensor2 = zeros(3,4); e_tensor2[3,1:2] .= -2.120582; e_tensor2[3,4] = 9.52183
    ϵ1 = diagm(0 => [3.72e-2, 3.72e-2, 3.72e-2])
    ϵ2 = diagm(0 => [4.065, 4.065, 2.079])

    piezo_mat1 = material_def([:elastic, :electric], 2, :out_of_plane,
                               C=C1, e=e_tensor1, ϵ=ϵ1)
    piezo_mat2 = material_def([:elastic, :electric], 2, :out_of_plane,
                               C=C2, e=e_tensor2, ϵ=ϵ2)

    # Assemble coupled stiffness
    K_piezo, F_u, F_ϕ = assemble_global_matrix(
        mesh.elements, mesh.nodes,
        Elem.type, Elem.int_order,
        [piezo_mat1, piezo_mat2],
        2, mesh.type_elem
    )
    isapprox(K_piezo, K_piezo', rtol=1e-16) || error("K_piezo is not symmetric")

    # Dirichlet BC cases
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
    bc_u = [HomogenizationBC(:periodic,val,(mesh.master,mesh.slave),nothing,nothing,mesh.nodes,3)
            for val in dirichlet_u]
    bc_ϕ = [HomogenizationBC(:periodic,val,(mesh.master,mesh.slave),nothing,nothing,mesh.nodes,3)
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
        @assert abs(W_int - W_ext) < 1e-10*abs(W_ext)
    end
    verify_energy(K_piezo, U4, F_u)
    verify_energy(K_piezo, U7, F_ϕ)

    # Visualization displacement and fields
    magn   = 0.3
    Ux     = U1[1:3:end]; Uy = U1[2:3:end]; Up = U1[3:3:end]
    Nx2    = mesh.nodes[:,1] .+ Ux*magn; Ny2 = mesh.nodes[:,2] .+ Uy*magn
    nodes2 = hcat(Nx2, Ny2)
    stress_strain = recover_field_values(
        mesh.elements, nodes2,
        [piezo_mat1, piezo_mat2],
        (U=U5, ϕ=U5),
        mesh.type_elem, Elem.type, Elem.int_order, 2
    )
    σ₁₁, σ₂₂, σ₁₂ = stress_strain.stress[:,1], stress_strain.stress[:,2], stress_strain.stress[:,3]
    ϵ₁₁, ϵ₂₂, ϵ₁₂ = stress_strain.strain[:,1], stress_strain.strain[:,2], stress_strain.strain[:,3]
    E₁, E₂, E₃   = stress_strain.elec[:,1], stress_strain.elec[:,2], stress_strain.elec[:,3]
    plot_champ_scalaire(nodes2, mesh.elements, E₁, Elem.type)
    plot_champ_scalaire(nodes2, mesh.elements, E₂, Elem.type)
    plot_champ_scalaire(nodes2, mesh.elements, E₃, Elem.type)
    plot_champ_scalaire(nodes2, mesh.elements, ϵ₁₁, Elem.type)
    plot_champ_scalaire(nodes2, mesh.elements, ϵ₂₂, Elem.type)
    plot_champ_scalaire(nodes2, mesh.elements, ϵ₁₂, Elem.type)

    # Compute effective piezoelectric properties
    solver_results = (
        U_mech = [U1, U2, U3, U4, U5, U6, U7],
        V_elec = [U1, U2, U3, U4, U5, U6, U7]
    )
    props,_ = compute_effective_property(
        [piezo_mat1, piezo_mat2], mesh.elements, mesh.nodes,
        mesh.type_elem, solver_results, Elem.type, Elem.int_order, 2
    )
    return (C=props.C, e=props.e, ϵ=props.ϵ)
end

# =============================================
# Poroelastic homogenization
# ============================================= 
function run_poro_case(mesh::MeshData, Elem::ElemData)
    # Material definition
    poro_mat = material_def([:elastic], 2, :isotropic,plane_stress = false, E=10, ν=1/3, α_p=0.7);

    # force initialisation
    border_elem = mesh.boundary_element;
    NodalForces = Dict(
    :inclusion => (fᵥ = nothing, fₜ = force_computation((x, y)  -> 1)));
    BoundaryFace = Dict(
    :inclusion => (
        element_border = border_elem[:inclusion_1],
        element_type = :Lin6,
        dim = 1,
        nodes = mesh.nodes,
        int_order = 5
    ));

    # Assemble coupled stiffness
    K_poro, F_u = assemble_global_matrix(
        mesh.elements, mesh.nodes,
        Elem.type, Elem.int_order,
        poro_mat,
        2, mesh.type_elem;
        NodalForces,BoundaryFace
    );


    isapprox(K_poro, K_poro', rtol=1e-16) || error("K_piezo is not symmetric")

    # Dirichlet BC cases
    dirichlet_u = [
        [1.0 0.0; 0.0 0.0],
        [0.0 0.0; 0.0 1.0],
        [0.0 0.5; 0.5 0.0],
        [0.0 0.0; 0.0 0.0]
    ];

    bc_u = [HomogenizationBC(:dirichlet,val,nothing,mesh.boundary,[true,true],mesh.nodes,2)
            for val in dirichlet_u];

     dirichlet_bc3 = BoundaryCondition(
        :dirichlet,
        mesh.boundary,
        [true, true],
        [0.0,0.0],
        2  # dofs_per_node
    );

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
    magn   = 0.3
    Ux     = U4[1:2:end]; Uy = U4[2:2:end];
    Nx2    = mesh.nodes[:,1] .+ Ux*magn; Ny2 = mesh.nodes[:,2] .+ Uy*magn;
    nodes2 = hcat(Nx2, Ny2);
    stress_strain = recover_field_values(
        mesh.elements, mesh.nodes,
        poro_mat,
        (U=U4,),
        mesh.type_elem, Elem.type, Elem.int_order, 2
    );
    σ₁₁, σ₂₂, σ₁₂ = stress_strain.stress[:,1], stress_strain.stress[:,2], stress_strain.stress[:,3];
    ϵ₁₁, ϵ₂₂, ϵ₁₂ = stress_strain.strain[:,1], stress_strain.strain[:,2], stress_strain.strain[:,3];

    plot_champ_scalaire(nodes2, mesh.elements, σ₁₁, Elem.type)
    plot_champ_scalaire(nodes2, mesh.elements, σ₂₂, Elem.type)
    plot_champ_scalaire(nodes2, mesh.elements, σ₁₂, Elem.type)
    plot_champ_scalaire(nodes2, mesh.elements, ϵ₁₁, Elem.type)
    plot_champ_scalaire(nodes2, mesh.elements, ϵ₂₂, Elem.type)
    plot_champ_scalaire(nodes2, mesh.elements, ϵ₁₂, Elem.type)

    # Compute effective piezoelectric properties
    # Compute effective stiffness
    solver_results = (U=(U1, U2, U3, U4),);
    
    props,volume = compute_effective_property(
        poro_mat,
        mesh.elements,
        mesh.nodes,
        mesh.type_elem,
        solver_results,
        Elem.type,
        Elem.int_order,
        2
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
    # elseif haskey(mat.properties, :α_p) && :pressure ∉ mat.type
    return (C=props.C, B=Β_eff, b_dilute=b, b_refine=br, solid_biot = i_n,s_biot_dilute =i_nd[1], s_biot_refine=i_nr[1] )
end

# =============================================
# stokes homogenization
# =============================================
#  Elem = ElemData(:Tri6, 2, 2);
#     mesh = setup_mesh(;
#         width=1.0, height=1.0,
#         volume_fraction=0.2, n_inclusions=1,
#         Elem, node_divisions=(5, 5), shape=:circle,
#         output_file="Test_plate_with_inclusions.msh",
#         voids = false
#     );
function run_stokesFlow_case(mesh::MeshData, Elem::ElemData)
    # Material definition
    C1 = [8.0 4.4 0 4.4; 4.4 8.0 0 4.4; 0 0 3.6 0; 4.4 4.4 0 8]
    C2 = [154.837 83.237 0 82.712; 83.237 154.837 0 82.712;
          0 0 35.8 0; 82.712 82.712 0 131.39]
    e_tensor1 = zeros(3,4)
    e_tensor2 = zeros(3,4); e_tensor2[3,1:2] .= -2.120582; e_tensor2[3,4] = 9.52183
    ϵ1 = diagm(0 => [3.72e-2, 3.72e-2, 3.72e-2])
    ϵ2 = diagm(0 => [4.065, 4.065, 2.079])

    piezo_mat1 = material_def([:elastic, :electric], 2, :out_of_plane,
                               C=C1, e=e_tensor1, ϵ=ϵ1)
    piezo_mat2 = material_def([:elastic, :electric], 2, :out_of_plane,
                               C=C2, e=e_tensor2, ϵ=ϵ2)

    # Assemble coupled stiffness
    K_piezo, F_u, F_ϕ = assemble_global_matrix(
        mesh.elements, mesh.nodes,
        Elem.type, Elem.int_order,
        [piezo_mat1, piezo_mat2],
        2, mesh.type_elem
    )
    isapprox(K_piezo, K_piezo', rtol=1e-16) || error("K_piezo is not symmetric")

    # Dirichlet BC cases
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
    bc_u = [HomogenizationBC(:periodic,val,(mesh.master,mesh.slave),nothing,nothing,mesh.nodes,3)
            for val in dirichlet_u]
    bc_ϕ = [HomogenizationBC(:periodic,val,(mesh.master,mesh.slave),nothing,nothing,mesh.nodes,3)
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
        @assert abs(W_int - W_ext) < 1e-10*abs(W_ext)
    end
    verify_energy(K_piezo, U4, F_u)
    verify_energy(K_piezo, U7, F_ϕ)

    # Visualization displacement and fields
    magn   = 0.3
    Ux     = U1[1:3:end]; Uy = U1[2:3:end]; Up = U1[3:3:end]
    Nx2    = mesh.nodes[:,1] .+ Ux*magn; Ny2 = mesh.nodes[:,2] .+ Uy*magn
    nodes2 = hcat(Nx2, Ny2)
    stress_strain = recover_field_values(
        mesh.elements, nodes2,
        [piezo_mat1, piezo_mat2],
        (U=U5, ϕ=U5),
        mesh.type_elem, Elem.type, Elem.int_order, 2
    )
    σ₁₁, σ₂₂, σ₁₂ = stress_strain.stress[:,1], stress_strain.stress[:,2], stress_strain.stress[:,3]
    ϵ₁₁, ϵ₂₂, ϵ₁₂ = stress_strain.strain[:,1], stress_strain.strain[:,2], stress_strain.strain[:,3]
    E₁, E₂, E₃   = stress_strain.elec[:,1], stress_strain.elec[:,2], stress_strain.elec[:,3]
    plot_champ_scalaire(nodes2, mesh.elements, E₁, Elem.type)
    plot_champ_scalaire(nodes2, mesh.elements, E₂, Elem.type)
    plot_champ_scalaire(nodes2, mesh.elements, E₃, Elem.type)
    plot_champ_scalaire(nodes2, mesh.elements, ϵ₁₁, Elem.type)
    plot_champ_scalaire(nodes2, mesh.elements, ϵ₂₂, Elem.type)
    plot_champ_scalaire(nodes2, mesh.elements, ϵ₁₂, Elem.type)

    # Compute effective piezoelectric properties
    solver_results = (
        U_mech = [U1, U2, U3, U4, U5, U6, U7],
        V_elec = [U1, U2, U3, U4, U5, U6, U7]
    )
    props = compute_effective_property(
        [piezo_mat1, piezo_mat2], mesh.elements, mesh.nodes,
        mesh.type_elem, solver_results, Elem.type, Elem.int_order, 2
    )
    return (C=props.C, e=props.e, ϵ=props.ϵ)
end
# Run script
main();


