using Test
using HomoLib: material_def, compute_effective_property,
               shape_data, build_B_matrices, jacobian_data

@testset "Effective Property (Thermal)" begin
    nodes = [
        0.0 0.0;
        1.0 0.0;
        0.0 1.0
    ]
    connectivity = [1 2 3]
    elements = [1 2 3]
    element_type = :Tri3
    dim = 2

    mat = material_def([:thermal], 2, :isotropic; κ=1.0)
    gauss_data = shape_data(element_type, 1, dim)
    jacobian_cache = jacobian_data(connectivity, nodes, gauss_data, dim)
    B_dicts = build_B_matrices(nodes, connectivity, mat, gauss_data, jacobian_cache)

    Geometric_Data = (
        gauss_data = gauss_data,
        jacobian_cache = jacobian_cache,
        B_dicts = B_dicts
    )

    solver_results = (U = (zeros(3), zeros(3)),)

    result, _ = compute_effective_property(
        mat,
        elements,
        nodes,
        nothing,
        solver_results,
        element_type,
        1,
        dim,
        Geometric_Data
    )

    @test size(result.K) == (2, 2)
end
