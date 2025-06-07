using Revise 
using HomoLib: material_def,Material
using LinearAlgebra, StaticArrays, SparseArrays, Tensors

using Test


@testset "Material Definitions Tests" begin
    # Isotropic Elasticity
    mat1 = material_def([:elastic], 2, :isotropic, E=200e9, ν=0.3)
    @test length(mat1.tensors) == 1
    @test mat1.B_types == [:strain]

    # Isotropic Thermal
    mat2 = material_def([:thermal], 2, :isotropic, κ=45.0)
    @test length(mat2.tensors) == 1
    @test mat2.B_types == [:temperature_gradient]

    # Isotropic Electric
    mat3 = material_def([:electric], 2, :isotropic, ϵ=8.85e-12)
    @test length(mat3.tensors) == 1
    @test mat3.B_types == [:electric_field]

    # Elastic + Thermal expansion
    mat4 = material_def([:elastic], 2, :isotropic, E=70e9, ν=0.33, α=2.5e-6)
    @test length(mat4.tensors) == 2
    @test :temperature in mat4.B_types

    # Viscoelastic with Cr
    Cr = Matrix(0.1I(3))
    mat5 = material_def([:elastic], 2, :isotropic, E=200e9, ν=0.3, Cr=Cr)
    @test length(mat5.tensors) == 3

    # Piezoelectricity
    e_tensor = rand(3, 2)
    mat6 = material_def([:elastic, :electric], 2, :isotropic, E=100e9, ν=0.25, e=e_tensor,ϵ=0.0)
    @test length(mat6.tensors) == 3
    @test :electric_field in mat6.B_types

    # Poroelasticity without flow
    mat7 = material_def([:elastic], 2, :isotropic, E=50e9, ν=0.25, α_p=0.7, M=1e9)
    @test length(mat7.tensors) == 3
    @test :pressure in mat7.B_types

    # Poroelasticity with flow
    K = Matrix(1.0I, 2, 2)
    mat8 = material_def([:elastic], 2, :isotropic, E=50e9, ν=0.25, α_p=0.7, M=1e9, K=K, μ_f=1e-3)
    @test length(mat8.tensors) == 4
    @test :fluid_velocity in mat8.B_types

    # Viscosity only
    mat9 = material_def([:viscous], 2, :isotropic, μ=10.0, λ=2.0)
    @test length(mat9.tensors) == 1
    @test mat9.B_types == [:strain_rate]

    # Plasticity
    mat10 = material_def([:elastic], 2, :isotropic, E=100e9, ν=0.3, σ_y=250e6, H=1e9)
    @test haskey(mat10.properties, :σ_y)

    # Viscoplasticity
    mat11 = material_def([:elastic, :viscous], 2, :isotropic, E=100e9, ν=0.3, μ=10.0, λ=2.0, σ_y=300e6)
    @test mat11.properties[:viscoplastic] == true
end

