
using Test
include("../src/HomoLib.jl")  # Adjust to actual path if different

# Basic square mesh (4 node quad)
nodes = [
    0.0 0.0;
    1.0 0.0;
    1.0 1.0;
    0.0 1.0
]
connectivity = [1 2 3 4]
dim = 2
order = 2

function run_physics_test(B_types::Vector{Symbol})
    mat = material_def(B_types, dim, :isotropic, E=1.0, ν=0.3, plane_stress=true, α=1e-5, T_ref=0.0)
    K = assemble_global_matrix(connectivity, nodes, :Quad4, order, mat, dim, nothing)
    @test issymmetric(K)
    @test all(isfinite, K)
end

@testset "Assembly Tests for 2D Isotropic (Elastic, Piezo, Thermoelastic, Thermal, etc.)" begin
    run_physics_test([:elastic])
    run_physics_test([:elastic, :electric_field])
    run_physics_test([:elastic, :temperature_gradient])
    run_physics_test([:temperature_gradient])
    run_physics_test([:strain])  # fallback to elastic
end
