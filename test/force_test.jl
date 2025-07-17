using HomoLib: compute_element_force,material_def
using HomoLib: compute_element_force, material_def, shape_data, jacobian_data, build_B_matrices
using Test
using LinearAlgebra, SparseArrays, Statistics
 
@testset "Element volume force" begin
    # Unit square Quad4 element
    nodes = [
        0.0 0.0;
        1.0 0.0;
        1.0 1.0;
        0.0 1.0
    ]
    connectivity = reshape(1:4, 1, 4)
    mat = material_def([:elastic], 2, :isotropic; E=1.0, ν=0.3, plane_stress=true)
    shp = shape_data(:Quad4, 2, 2)
    jac = jacobian_data(connectivity, nodes, shp, 2)
    B_dicts = build_B_matrices(nodes, Matrix(connectivity), mat, shp, jac)
    forces = (fᵥ = [1.0, 2.0], fₜ = nothing)
    f = compute_element_force(Val(:elastic), mat,
        jac[1], shp, nothing, nothing,
        forces, nothing, connectivity[1, :], B_dicts[1], nothing)
    expected = repeat([0.25, 0.5], 4)
    @test isapprox(f, expected; atol=1e-8)
end