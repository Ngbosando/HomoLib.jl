using Revise
using HomoLib: material_def, assemble_global_matrix
using Test

@testset "assemble_global_matrix: two Tri3 elements" begin
    nodes = [
        0.0 0.0;
        1.0 0.0;
        0.0 1.0;
        1.0 1.0
    ]
    connectivity = [
        1 2 3;
        2 4 3
    ]

    dim = 2
    order = 1
    element_type = :Tri3

    mat = material_def([:elastic], dim, :isotropic; E=1.0, ν=0.3, plane_stress=true)

    K_full = assemble_global_matrix(connectivity, nodes, element_type, order, mat, dim, nothing)

    K1 = assemble_global_matrix(connectivity[1:1, :], nodes, element_type, order, mat, dim, nothing)
    K2 = assemble_global_matrix(connectivity[2:2, :], nodes, element_type, order, mat, dim, nothing)
    K_expected = K1 + K2

    n_nodes = size(nodes, 1)
    @test size(K_full) == (n_nodes * dim, n_nodes * dim)
    @test Matrix(K_full) ≈ Matrix(K_expected)
    @test issymmetric(Matrix(K_full))

    translation = repeat([1.0, 1.0], n_nodes)
    @test norm(Matrix(K_full) * translation) < 1e-12
end

