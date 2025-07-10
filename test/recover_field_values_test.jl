
using Test
using HomoLib: recover_field_values, material_def

@testset "Recover Field Values Simple Mesh" begin
    # 2D triangular single element
    coords = [0.0 0.0; 1.0 0.0; 0.0 1.0]
    conn = reshape([1,2,3], 1, 3)
    element_type = :Tri3
    order = 1
    dim = 2
    mat = material_def([:elastic], dim, :isotropic; E=200.0, ν=0.3, plane_stress=true)

    α = 1e-3
    β = 2e-3
    Uvec = [0.0,0.0, α,0.0, 0.0,β]
    Uresult = (U = Uvec,)

    recovered = recover_field_values(conn, coords, mat, Uresult, nothing,
                                     element_type, order, dim)

    for i in 1:3
        @test isapprox(recovered.strain[i, 1], α; atol=1e-12)
        @test isapprox(recovered.strain[i, 2], β; atol=1e-12)
        @test isapprox(recovered.strain[i, 3], 0.0; atol=1e-12)
    end

    C = (mat.properties[:E]/(1 - mat.properties[:ν]^2)) *
        [1.0 mat.properties[:ν] 0.0;
         mat.properties[:ν] 1.0 0.0;
         0.0 0.0 (1 - mat.properties[:ν])/2]
    σ_expected = C * [α, β, 0.0]

    for i in 1:3
        @test isapprox(recovered.stress[i, 1], σ_expected[1]; atol=1e-10)
        @test isapprox(recovered.stress[i, 2], σ_expected[2]; atol=1e-10)
        @test isapprox(recovered.stress[i, 3], σ_expected[3]; atol=1e-10)
    end
end
