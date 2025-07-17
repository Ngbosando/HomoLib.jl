using HomoLib: jacobian_data,shape_data,material_def,build_B_matrices
using Test
using LinearAlgebra

@testset "Precompute Data Tests" begin
    # Simple square 2D domain for testing Quad4
    nodes = [
        0.0 0.0;
        1.0 0.0;
        1.0 1.0;
        0.0 1.0
    ]

    connectivity = [
        1 2 3 4
    ]

    material = material_def([:elastic], 2, :isotropic, 
                            E=1.0, ν=0.0, plane_stress=true)


    element_type = :Quad4
    int_order = 2

    # Step 1: Shape data
    gauss = shape_data(element_type, int_order, material.dim)
    @test length(gauss.weights) > 0
    @test length(gauss.shape_ξ) == length(gauss.weights)

    # Step 2: Jacobian
    jac_data = jacobian_data(connectivity, nodes, gauss, material.dim)
    @test length(jac_data) == 1
    for (detJ, invJ) in jac_data[1]
        @test isapprox(detJ, 0.25; atol=1e-10)  # Area of 1 for unit square split into 4 Gauss pts
        @test isapprox(det(invJ), 1/detJ; atol=1e-10)
    end

    # Step 3: B matrices
    B = build_B_matrices(nodes, connectivity, material, gauss, jac_data)
    @test length(B) == 1
    for Btype in material.B_types
        Bmat = B[1][Btype]
        @test length(Bmat) == length(gauss.weights)
        for mat in Bmat
            @test size(mat, 2) == material.dim * size(connectivity, 2)
        end
    end
end
