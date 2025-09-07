using Test
using HomoLib: ThermalMaterial, compute_effective_property,
               precompute_geometric_data, assemble_KMF

@testset "Effective Property (Thermal)" begin
    nodes = [
        0.0 0.0;
        1.0 0.0;
        0.0 1.0
    ]
    connectivity = [1 2 3]
    element_type = :Tri3
    dim = 2

    mat = ThermalMaterial(2, k=1.0)

    solver_results = (U = (zeros(3), zeros(3)),)
    Geometric_Data = precompute_geometric_data(
        element_type, 1, dim,
        connectivity, nodes, mat
    );
    result, _ = compute_effective_property(
        mat,
        connectivity,
        solver_results,
        dim,
        Geometric_Data
    )

    @test size(result.K) == (2, 2)
    @test result.K ==[0.0 0.0;0.0 0.0]
end
