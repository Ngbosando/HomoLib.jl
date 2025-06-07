using HomoLib: compute_element_force,NodalForces,BoundaryFace
using Test
using FastGaussQuadrature: gausslegendre

# @testset "Test traction surfacique sur un élément bord (Line2)" begin
    # --- Création du matériau avec material_def ---
    material = material_def(:elastic, 2, :isotropic, E=100.0, ν=0.3, plane_stress=true)

    # --- Géométrie de la face ---
    face_conn = [1, 2]
    face_nodes = [0.0 0.0; 1.0 0.0]

    # --- Traction appliquée ---
    fₜ = [100.0, 0.0]

    # --- Gauss pour Line2 (ordre 2) ---
    ξ = [-1/sqrt(3), 1/sqrt(3)]
    w = [1.0, 1.0]
    N = [(0.5 * (1 - ξi), 0.5 * (1 + ξi)) for ξi in ξ]
    gauss_fdata = (; N = N, weights = w)

    # --- Jacobien (longueur / 2 = 0.5) ---
    jacobian_fdata = [[(0.5, nothing), (0.5, nothing)]]

    # --- Forces (aucune volumique) ---
    forces = NodalForces([0.0, 0.0], fₜ)

    # --- Calcul force élémentaire ---
    f = compute_element_force(
        Val(:elastic), material, 1,
        nothing, nothing,
        jacobian_fdata, gauss_fdata,
        forces,
        face_conn, face_conn
    )

    # --- Vérification des résultats ---
    @test isapprox(f[1], 50.0; atol=1e-6)  # Fx node 1
    @test isapprox(f[2], 0.0; atol=1e-6)   # Fy node 1
    @test isapprox(f[3], 50.0; atol=1e-6)  # Fx node 2
    @test isapprox(f[4], 0.0; atol=1e-6)   # Fy node 2

    @test isapprox(sum(f[1:2:end]), 100.0; atol=1e-6)
    @test isapprox(sum(f[2:2:end]), 0.0; atol=1e-6)
# end
@testset "Test traction surfacique sur bord Line6" begin
    # Matériau
    mat = material_def(:elastic, 2, :isotropic, E=1.0, ν=0.3)

    # Connexion nodale fictive
    face_conn = 1:6

    # Coordonnées des 6 nœuds sur ligne [0,1]
    face_nodes = hcat(collect(LinRange(0, 1, 6)), zeros(6))'  # (x, y) pour chaque nœud

    # Traction appliquée
    fₜ = [60.0, 0.0]

    # Intégration de surface pour Line6
    ξ, w = gausslegendre(6)
    Nvals = [shape_functions(:Lin5, ξi) for ξi in ξ]  # Line5 = 6 nœuds

    gauss_fdata = (; N = Nvals, weights = w)

    # Jacobien : longueur = 1.0 → J = 0.5 (sur ref [-1,1])
    jacobian_fdata = [ [(0.5, nothing) for _ in ξ] ]

    # Forces (fᵥ = 0, fₜ ≠ 0, pas de pression)
    forces = NodalForces([0.0, 0.0], fₜ)

    # Calcul
    f = compute_element_force(
        Val(:elastic), mat, 1,
        nothing, nothing,
        jacobian_fdata, gauss_fdata,
        forces,
        collect(face_conn), collect(face_conn)
    )

    # --- Vérifications ---
    Fx = f[1:2:end]
    Fy = f[2:2:end]

    @test isapprox(sum(Fx), 60.0; atol=1e-6)
    @test isapprox(sum(Fy), 0.0; atol=1e-6)

    @info "Force nodale par nœud (Fx) = $(round.(Fx, digits=4))"
end
