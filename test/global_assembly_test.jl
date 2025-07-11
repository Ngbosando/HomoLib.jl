using Revise
using HomoLib: material_def, assemble_global_matrix,
            shape_data,jacobian_data,build_B_matrices
using LinearAlgebra, Statistics
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
    # Precompute data no poro
    gauss_data = shape_data(element_type, order, dim)
    jacobian_cache = jacobian_data(connectivity, nodes, gauss_data, dim)
    B_dicts = build_B_matrices(nodes, connectivity, mat, gauss_data, jacobian_cache)
    Geometric_Data = (gauss_data = gauss_data, jacobian_cache = jacobian_cache, B_dicts = B_dicts)

    K_full = assemble_global_matrix(connectivity, nodes, element_type, order, mat, dim, nothing, Geometric_Data)

    # For K1 and K2, build geometric data for each element separately
    conn1 = connectivity[1:1, :]
    conn2 = connectivity[2:2, :]

    jacobian_cache1 = jacobian_data(conn1, nodes, gauss_data, dim)
    B_dicts1 = build_B_matrices(nodes, conn1, mat, gauss_data, jacobian_cache1)
    Geometric_Data1 = (gauss_data = gauss_data, jacobian_cache = jacobian_cache1, B_dicts = B_dicts1)

    jacobian_cache2 = jacobian_data(conn2, nodes, gauss_data, dim)
    B_dicts2 = build_B_matrices(nodes, conn2, mat, gauss_data, jacobian_cache2)
    Geometric_Data2 = (gauss_data = gauss_data, jacobian_cache = jacobian_cache2, B_dicts = B_dicts2)

    K1 = assemble_global_matrix(conn1, nodes, element_type, order, mat, dim, nothing, Geometric_Data1)
    K2 = assemble_global_matrix(conn2, nodes, element_type, order, mat, dim, nothing, Geometric_Data2)
    K_expected = K1 + K2

    n_nodes = size(nodes, 1)
    @test size(K_full) == (n_nodes * dim, n_nodes * dim)
    @test Matrix(K_full) ≈ Matrix(K_expected)
    # @test issymmetric(Matrix(K_full))
    @test isapprox(K_full, K_full', rtol=1e-100)
    translation = repeat([1.0, 1.0], n_nodes)
    @test norm(Matrix(K_full) * translation) < 1e-12
end

