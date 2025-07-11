
using Test
using HomoLib: recover_field_values, material_def,
            shape_data,jacobian_data,build_B_matrices,
            shape_data,jacobian_data,build_B_matrices,
            get_node_wise_permutation
   
@testset "Recover Field Values Thermal" begin
    coords = [0.0 0.0; 1.0 0.0; 0.0 1.0]
    conn = reshape([1,2,3], 1, 3)
    element_type = :Tri3
    order = 1
    dim = 2
    mat = material_def([:thermal], dim, :isotropic; κ=1.5)
    # Precompute data
    gauss_data = shape_data(element_type, order, dim)
    jacobian_cache = jacobian_data(conn, coords, gauss_data, dim)
    B_dicts = build_B_matrices(coords, conn, mat, gauss_data, jacobian_cache)
    Geometric_Data = (gauss_data = gauss_data,jacobian_cache = jacobian_cache,B_dicts = B_dicts )
    a = 2.0
    b = 3.0
    Tvec = [0.0, a, b]
    Uresult = (T = Tvec,)

    recovered = recover_field_values(conn, coords, mat, Uresult, nothing,
                                     element_type, order, dim,Geometric_Data)

    expected_grad = [a, b]
    expected_flux = -1.5 .* expected_grad

    for i in 1:3
        @test isapprox(recovered.grad_temp[i,1], expected_grad[1]; atol=1e-12)
        @test isapprox(recovered.grad_temp[i,2], expected_grad[2]; atol=1e-12)
        @test isapprox(recovered.flux[i,1], expected_flux[1]; atol=1e-12)
        @test isapprox(recovered.flux[i,2], expected_flux[2]; atol=1e-12)
    end
end


@testset "Recover Field Values Piezoelectric" begin
    coords = [0.0 0.0; 1.0 0.0; 0.0 1.0]
    conn = reshape([1,2,3], 1, 3)
    element_type = :Tri3
    order = 1
    dim = 2

    e_tensor = [0.2 0.1 0.0;
                -0.1 0.3 0.0]
    ϵ_tensor = Matrix{Float64}(I, 2, 2)
    mat = material_def([:elastic, :electric], dim, :isotropic;
                       E=1.0, ν=0.0, plane_stress=true,
                       e=e_tensor, ϵ=ϵ_tensor)
    # Precompute data
    gauss_data = shape_data(element_type, order, dim)
    jacobian_cache = jacobian_data(conn, coords, gauss_data, dim)
    B_dicts = build_B_matrices(coords, conn, mat, gauss_data, jacobian_cache)
    Geometric_Data = (gauss_data = gauss_data,jacobian_cache = jacobian_cache,B_dicts = B_dicts )                   

    α = 1e-3
    β = 2e-3
    c = 0.5
    d = -0.1
    Uvec = [0.0, 0.0, α, 0.0, 0.0, β]
    ϕvec = [0.0, c, d]
    Uϕ = vcat(Uvec, ϕvec)
    nn = div(length(Uϕ),3)
    perm = get_node_wise_permutation(mat, nn)

    Uresult = (U = Uϕ[perm], ϕ = Uϕ[perm])
    recovered = recover_field_values(conn, coords, mat, Uresult, nothing,
                                     element_type, order, dim, Geometric_Data)

    ε_expected = [α, β, 0.0]
    E_expected = [c, d]
    σ_expected = mat.tensors[1] * ε_expected - transpose(e_tensor) * E_expected
    D_expected = e_tensor * ε_expected + ϵ_tensor * E_expected

    for i in 1:3
        @test isapprox(recovered.strain[i,1], ε_expected[1]; atol=1e-12)
        @test isapprox(recovered.strain[i,2], ε_expected[2]; atol=1e-12)
        @test isapprox(recovered.strain[i,3], ε_expected[3]; atol=1e-12)

        @test isapprox(recovered.stress[i,1], σ_expected[1]; atol=1e-12)
        @test isapprox(recovered.stress[i,2], σ_expected[2]; atol=1e-12)
        @test isapprox(recovered.stress[i,3], σ_expected[3]; atol=1e-12)

        @test isapprox(recovered.elec[i,1], E_expected[1]; atol=1e-12)
        @test isapprox(recovered.elec[i,2], E_expected[2]; atol=1e-12)

        @test isapprox(recovered.elec_disp[i,1], D_expected[1]; atol=1e-12)
        @test isapprox(recovered.elec_disp[i,2], D_expected[2]; atol=1e-12)
    end
end