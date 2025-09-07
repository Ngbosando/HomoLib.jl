using HomoLib: ElasticMaterial,ThermalMaterial,PiezoelectricMaterial,PoroelasticMaterial,StokesMaterial,
                get_dofs_per_node,get_dofs_per_node
using Test
using LinearAlgebra


@testset "Material Definitions Tests" begin
    # 1. Elastic Material Tests
    @testset "Elastic Materials" begin
        # 1D isotropic elastic
        mat_1d = ElasticMaterial(1, E=200e9, ν=0.3, ρ=2700)
        @test length(mat_1d.tensors) == 1
        @test mat_1d.B_types == [:strain]
        @test get_dofs_per_node(mat_1d) == 1
        @test mat_1d.mass_properties[:ρ] ≈ 2700

        # 2D isotropic elastic (plane strain)
        mat_2d = ElasticMaterial(2, E=70e9, ν=0.33, ρ=2700)
        @test length(mat_2d.tensors) == 1
        @test size(mat_2d.tensors[1]) == (3, 3)
        @test get_dofs_per_node(mat_2d) == 2

        # 2D out-of-plane elastic
        mat_2d_oop = ElasticMaterial(2, symmetry=:out_of_plane, E=70e9, ν=0.33, ρ=2700)
        @test size(mat_2d_oop.tensors[1]) == (4, 4)

        # 3D isotropic elastic
        mat_3d = ElasticMaterial(3, E=200e9, ν=0.3, ρ=7800)
        @test length(mat_3d.tensors) == 1
        @test size(mat_3d.tensors[1]) == (6, 6)
        @test get_dofs_per_node(mat_3d) == 3

        # Orthotropic elastic
        C_ortho = [200e9 50e9 0; 50e9 100e9 0; 0 0 30e9]
        mat_ortho = ElasticMaterial(2, symmetry=:orthotropic, C=C_ortho, ρ=2700)
        @test mat_ortho.tensors[1] ≈ C_ortho
    end

    # 2. Thermal Material Tests
    @testset "Thermal Materials" begin
        # 1D isotropic thermal
        mat_1d = ThermalMaterial(1, k=50.0, ρ=2000, c=900)
        @test length(mat_1d.tensors) == 1
        @test mat_1d.B_types == [:temperature_gradient]
        @test get_dofs_per_node(mat_1d) == 1
        @test mat_1d.mass_properties[:ρ] ≈ 2000
        @test mat_1d.mass_properties[:c] ≈ 900

        # 2D isotropic thermal
        mat_2d = ThermalMaterial(2, k=45.0, ρ=2000, c=900)
        @test length(mat_2d.tensors) == 1
        @test size(mat_2d.tensors[1]) == (2, 2)
        @test get_dofs_per_node(mat_2d) == 1

        # 2D anisotropic thermal
        κ_ortho = [60.0 0; 0 30.0]
        mat_2d_aniso = ThermalMaterial(2, symmetry=:orthotropic, κ=κ_ortho, ρ=2000, c=900)
        @test mat_2d_aniso.tensors[1] ≈ κ_ortho

        # 3D isotropic thermal
        mat_3d = ThermalMaterial(3, k=50.0, ρ=2000, c=900)
        @test length(mat_3d.tensors) == 1
        @test size(mat_3d.tensors[1]) == (3, 3)
    end

    # 3. Piezoelectric Material Tests
    @testset "Piezoelectric Materials" begin
        # 2D piezoelectric
        C_piezo = [120e9 75e9 0; 75e9 110e9 0; 0 0 20e9]
        ε_piezo = [1.3e-8 0; 0 1.1e-8]
        e_piezo = [0 0 -0.59; -0.61 0 0]
        
        mat_2d = PiezoelectricMaterial(2, C=C_piezo, ε=ε_piezo, e=e_piezo, ρ=7600)
        @test length(mat_2d.tensors) == 3
        @test mat_2d.B_types == [:strain, :electric_field]
        @test get_dofs_per_node(mat_2d) == 3  # 2 displacement + 1 electric
        @test mat_2d.tensors[1] ≈ C_piezo
        @test mat_2d.tensors[2] ≈ ε_piezo
        @test mat_2d.tensors[3] ≈ e_piezo

        # 3D piezoelectric
        C_3d = [166e9 77e9 78e9 0 0 0;
                77e9 166e9 78e9 0 0 0;
                78e9 78e9 162e9 0 0 0;
                0 0 0 43e9 0 0;
                0 0 0 0 43e9 0;
                0 0 0 0 0 44.5e9]
        ε_3d = [1.53e-8 0 0; 0 1.53e-8 0; 0 0 1.5e-8]
        e_3d = [0 0 0 0 12.0 0;
                0 0 0 12.0 0 0;
                -4.1 -4.1 14.1 0 0 0]
        
        mat_3d = PiezoelectricMaterial(3, C=C_3d, ε=ε_3d, e=e_3d, ρ=7500)
        @test get_dofs_per_node(mat_3d) == 4  # 3 displacement + 1 electric
    end

    # 4. Poroelastic Material Tests
    @testset "Poroelastic Materials" begin
        # 2D poroelastic
        mat_2d = PoroelasticMaterial(2, E=50e9, ν=0.25, α_p=0.7, ρs=2000, ρf=1000)
        @test length(mat_2d.tensors) == 1  # Only C tensor
        @test mat_2d.B_types == [:strain]  # Pressure handled separately
        @test get_dofs_per_node(mat_2d) == 2  # 2 displacement 
        @test mat_2d.α_p ≈ 0.7
        @test mat_2d.mass_properties[:ρs] ≈ 2000
        @test mat_2d.mass_properties[:ρf] ≈ 1000

        # 3D poroelastic
        mat_3d = PoroelasticMaterial(3, E=50e9, ν=0.25, α_p=0.7, ρs=2000, ρf=1000)
        @test get_dofs_per_node(mat_3d) == 3  # 3 displacement, pressure handled separately
    end

    # 5. Stokes Material Tests
    @testset "Stokes Materials" begin
        # 2D Stokes
        mat_2d = StokesMaterial(2, μ=1e-3, ρ=1000)
        @test length(mat_2d.tensors) == 0  # No tensors stored
        @test mat_2d.B_types == [:velocity_gradient, :pressure]
        @test get_dofs_per_node(mat_2d) == 3  # 2 velocity + 1 pressure
        @test mat_2d.μ ≈ 1e-3
        @test mat_2d.mass_properties[:ρ] ≈ 1000

        # 3D Stokes
        mat_3d = StokesMaterial(3, μ=1e-3, ρ=1000)
        @test get_dofs_per_node(mat_3d) == 4  # 3 velocity + 1 pressure
    end

    # 6. Custom Properties Tests
    @testset "Custom Properties" begin
        # Elastic with custom properties
        mat = ElasticMaterial(2, E=200e9, ν=0.3, ρ=2700)
        mat.properties[:yield_stress] = 250e6
        mat.properties[:thermal_expansion] = 2.5e-6
        
        @test mat.properties[:yield_stress] ≈ 250e6
        @test mat.properties[:thermal_expansion] ≈ 2.5e-6

        # Thermal with custom properties
        mat_thermal = ThermalMaterial(2, k=45.0, ρ=2000, c=900)
        mat_thermal.properties[:melting_point] = 1600.0
        @test mat_thermal.properties[:melting_point] ≈ 1600.0
    end

    # 7. Type and Physics Tests
    @testset "Type and Physics" begin
        mat_elastic = ElasticMaterial(2, E=200e9, ν=0.3, ρ=2700)
        mat_thermal = ThermalMaterial(2, k=45.0, ρ=2000, c=900)
        mat_piezo = PiezoelectricMaterial(2, C=ones(3,3), ε=ones(2,2), e=ones(2,3), ρ=2700)
        
        @test mat_elastic.type == [:elastic]
        @test mat_thermal.type == [:thermal]
        @test mat_piezo.type == [:elastic, :electric]
    end
end