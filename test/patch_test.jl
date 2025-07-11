using Revise
using HomoLib: assemble_global_matrix, material_def,
            shape_data,jacobian_data,build_B_matrices
using Test

@testset "3D elastic Hex8 patch" begin
    nodes = [
        0.0 0.0 0.0;
        1.0 0.0 0.0;
        1.0 1.0 0.0;
        0.0 1.0 0.0;
        0.0 0.0 1.0;
        1.0 0.0 1.0;
        1.0 1.0 1.0;
        0.0 1.0 1.0
    ]
    # hex_coords(n) = [collect((x, y, z)) for z in range(-1.0, 1.0; length=n),
    #                                      y in range(-1.0, 1.0; length=n),
    #                                      x in range(-1.0, 1.0; length=n)] |> vec
    # nodes = hcat(hex_coords(2)...)'  # Matrix{Float64} of size (n^3, 3)

    # Build connectivity for Hex8 elements in a structured 2x2x2 grid
    connectivity = Matrix(reshape(1:8, :, 8))
    mat = material_def([:elastic], 3, :isotropic; E=1.0, ν=0.3)
    # Precompute data
    gauss_data = shape_data(:Hex8, 2, 3)
    jacobian_cache = jacobian_data(connectivity, nodes, gauss_data, 3)
    B_dicts = build_B_matrices(nodes, connectivity, mat, gauss_data, jacobian_cache)
    Geometric_Data = (gauss_data = gauss_data,jacobian_cache = jacobian_cache,B_dicts = B_dicts )
    K = assemble_global_matrix(Matrix(connectivity), nodes, :Hex8, 3, mat, 3, nothing, Geometric_Data)
    ε = 1.0e-3
    U = Float64[]
    for (x,y,z) in eachrow(nodes)
        push!(U, ε)
        push!(U, 0.0)
        push!(U, 0.0)
    end
    res = K * U
    @test norm(res) < 1e-6
end

@testset "2D piezoelectric Quad4 patch" begin
    nodes = [
        0.0 0.0;
        1.0 0.0;
        1.0 1.0;
        0.0 1.0
    ]
    connectivity = reshape(1:4, 1, 4)
    mat = material_def([:elastic, :electric], 2, :isotropic;
                        E=1.0, ν=0.3, plane_stress=true,
                        e=zeros(2,3), ϵ=Matrix(0.1I(2)))
                        # Precompute data
    gauss_data = shape_data(:Quad4, 1, 2)
    jacobian_cache = jacobian_data(Matrix(connectivity), nodes, gauss_data, 2)
    B_dicts = build_B_matrices(nodes, Matrix(connectivity), mat, gauss_data, jacobian_cache)
    Geometric_Data = (gauss_data = gauss_data,jacobian_cache = jacobian_cache,B_dicts = B_dicts )
    K = assemble_global_matrix(Matrix(connectivity), nodes, :Quad4, 1, mat, 2, nothing, Geometric_Data)
    α = 1.0e-3
    U_mech = Float64[]
    for (x,y) in eachrow(nodes)
        push!(U_mech, α)
        push!(U_mech, 0.0)
    end
    ϕ = zeros(4)
    U = vcat(U_mech, ϕ)
    using HomoLib:  get_node_wise_permutation
    perm = get_node_wise_permutation(mat, 4)
    U = U[perm]
    res = K * U
    @test norm(res) < 1e-6
end

