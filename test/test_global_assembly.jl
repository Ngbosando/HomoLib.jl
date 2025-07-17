using HomoLib:material_def,assemble_global_matrix
using LinearAlgebra
using CairoMakie
using Test



    @testset "Global Matrix Assembly Additional Tests" begin
        # 2D Tests: Triangle element (Tri3)
        nodes = [0.0 0.0;
                 1.0 0.0;
                 0.0 1.0]
        connectivity = [1 2 3]
        dim = 2
        order = 1
        element_type = :Tri3

        # A helper function that includes extra checks for the stiffness matrix
        function run_test_extra(B_types::Vector{Symbol}, symmetry; kwargs...)
            mat = material_def(B_types, dim, symmetry; kwargs...)
            # Precompute data
            gauss_data = shape_data(element_type, order, dim)
            jacobian_cache = jacobian_data(connectivity, nodes, gauss_data, dim)
            B_dicts = build_B_matrices(nodes, connectivity, mat, gauss_data, jacobian_cache)
            Geometric_Data = (gauss_data = gauss_data,jacobian_cache = jacobian_cache,B_dicts = B_dicts )
            result = assemble_global_matrix(connectivity, nodes, element_type, order, mat, dim, nothing,Geometric_Data)
            K = isa(result, Tuple) ? result[1] : result
            @test size(K, 1) == size(K, 2)
            @test issymmetric(Matrix(K))
            @test norm(K) > 1e-12
            @test all(isfinite, K)
        end

        # Existing tests with additional internal checks
        run_test_extra([:elastic], :out_of_plane; E=1.0, ν=0.3, plane_stress=true, α=1e-5)
        run_test_extra([:elastic, :electric], :isotropic; E=1.0, ν=0.3, plane_stress=true, e=zeros(2,3), ϵ=[1.0 0;0 0])
        run_test_extra([:thermal], :isotropic; κ=1)
        Cr = Matrix(0.1I(3))
        run_test_extra([:elastic], :isotropic; E=200e9, ν=0.3, Cr=Cr)

        # 3D Test: Tetrahedral element (Tet4) assumed available in HomoLib
        nodes3D = [0.0 0.0 0.0;
                   1.0 0.0 0.0;
                   0.0 1.0 0.0;
                   0.0 0.0 1.0]
        connectivity3D = [1 2 3 4]
        dim3 = 3
        order3 = 1
        element_type3 = :Tet4  # This element should be supported by HomoLib

        function run_test_3D(B_types::Vector{Symbol}, symmetry; kwargs...)
            mat = material_def(B_types, dim3, symmetry; kwargs...)
            # Precompute data
            gauss_data = shape_data(element_type3, order3, dim3)
            jacobian_cache = jacobian_data(connectivity3D, nodes3D, gauss_data, dim3)
            B_dicts = build_B_matrices(nodes3D, connectivity3D, mat, gauss_data, jacobian_cache)
            Geometric_Data = (gauss_data = gauss_data,jacobian_cache = jacobian_cache,B_dicts = B_dicts )
            result = assemble_global_matrix(connectivity3D, nodes3D, element_type3, order3, mat, dim3, nothing,Geometric_Data)
            K = isa(result, Tuple) ? result[1] : result
            @test size(K, 1) == size(K, 2)
            @test issymmetric(Matrix(K))
            @test norm(K) > 1e-12
            @test all(isfinite, K)
        end
        run_test_3D([:elastic], :isotropic; E=1.0, ν=0.3)

        # Error handling: When a non-existent material type is passed, expect an error.
        # @test_throws ArgumentError material_def([:nonexistent], dim, :isotropic)

        # Dummy check: Verify that the stiffness matrix has (nearly) nonnegative eigenvalues.
        function dummy_matrix_check(B_types::Vector{Symbol}, symmetry; kwargs...)
            mat = material_def(B_types, dim, symmetry; kwargs...)
             # Precompute data
            gauss_data = shape_data(element_type, order, dim)
            jacobian_cache = jacobian_data(connectivity, nodes, gauss_data, dim)
            B_dicts = build_B_matrices(nodes, connectivity, mat, gauss_data, jacobian_cache)
            Geometric_Data = (gauss_data = gauss_data,jacobian_cache = jacobian_cache,B_dicts = B_dicts )
            result = assemble_global_matrix(connectivity, nodes, element_type, order, mat, dim, nothing,Geometric_Data)
            K = isa(result, Tuple) ? result[1] : result
            s = eigvals(Matrix(K))
            @test minimum(real(s)) > -1e-10
        end
        dummy_matrix_check([:elastic], :out_of_plane; E=1.0, ν=0.3, plane_stress=true, α=1e-5)
    end
  
