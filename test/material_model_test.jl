using Revise 
using HomoLib: material_def,Material
using LinearAlgebra, StaticArrays, SparseArrays, Tensors

using Test


@testset "Material Definitions Tests" begin
    # Isotropic Elasticity
    mat1 = material_def([:elastic], 2, :isotropic, E=200e9, ν=0.3,mass_properties=Dict(:ρ=>2700))
    @test length(mat1.tensors) == 1
    @test mat1.B_types == [:strain]

    # Isotropic Thermal
    mat2 = material_def([:thermal], 2, :isotropic, κ=45.0)
    @test length(mat2.tensors) == 1
    @test mat2.B_types == [:temperature_gradient]

    # Isotropic Electric
    ϵ_tensor= Matrix(I(2)*8.85e-12)
    mat3 = material_def([:electric], 2, :isotropic, ϵ=ϵ_tensor)
    @test length(mat3.tensors) == 1
    @test mat3.B_types == [:electric_field]

    # Elastic + Thermal expansion
    mat4 = material_def([:elastic], 2, :isotropic, E=70e9, ν=0.33, α=2.5e-6)
    @test length(mat4.tensors) == 2
    @test length(mat4.tensors[2]) == 3
    @test :strain in mat4.B_types

    # Piezoelectricity
    e_tensor = rand(2, 3)
    ϵ_tensor = zeros(2,2)
    mat5 = material_def([:elastic, :electric], 2, :isotropic, E=100e9, ν=0.25, e=e_tensor,ϵ=ϵ_tensor)
    @test length(mat5.tensors) == 3
    @test :electric_field in mat5.B_types

    # Poroelasticity 
    mat6 = material_def([:elastic], 2, :isotropic, E=50e9, ν=0.25, α_p=[0.7,0])
    @test length(mat6.tensors) == 2
    @test :pressure in mat6.B_types

    # Biot Poroelasticity 
    mat7 = material_def([:fluid, :elastic], 2, :isotropic, E=50e9, ν=0.25, α_p=[0.7,0], M=1e9, μ_f=1e-3)
    @test length(mat7.tensors) == 4
    @test :velocity_gradient in mat7.B_types

    # stokes flow
    mat8 = material_def([:fluid],2,:isotropic,α_p=[0.7,0], μ_f=1e-3 )
    @test length(mat8.tensors) == 2
    @test :pressure in mat8.B_types

end

