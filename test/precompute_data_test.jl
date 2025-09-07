using HomoLib: ElasticMaterial, PiezoelectricMaterial, compute_shape_function_data, 
                 compute_element_jacobian_data, compute_differential_operator_matrices,
                 BIndexable, BAtElem
using Test
using LinearAlgebra
using StaticArrays

function create_test_material(type::Symbol)
    if type == :elastic
        return ElasticMaterial(2, E=1.0, ν=0.0, ρ=1.0)
    elseif type == :piezoelectric

        C = ones(4, 4); e = ones(3, 4); eps = ones(3, 3)
        return PiezoelectricMaterial(2, :out_of_plane, Dict(), C, e, eps)
    else
        error("Test material type not recognized.")
    end
end


@testset "Precompute Data Tests" begin
    # Rigid Body Motion Test 
    nodes = [0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0]
    connectivity = [1 2 3 4]
    element_type = :Quad4
    int_order = 2
    dim = 2
    #  Pre-calculate common data 
    shape_data = compute_shape_function_data(element_type, int_order, dim, size(connectivity, 2))
    jac_data = compute_element_jacobian_data(connectivity, nodes, shape_data, Val(dim))
    shape_data = compute_shape_function_data(element_type, int_order, dim, size(connectivity, 2))
    jac_data = compute_element_jacobian_data(connectivity, nodes, shape_data, Val(dim))
    @testset "Rigid Body Motion Test (Elastic)" begin
        material = create_test_material(:elastic)
        B_storage = compute_differential_operator_matrices(material, shape_data, jac_data)
        B_elem = B_storage[1]

        # Pure translation in x-direction 
        u_translation = [1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0]
        
        for qp in 1:length(shape_data.weights)
            B_voigt = B_elem.voigt_gradient_operator[:,:,qp]
            strain = B_voigt * u_translation
            @test all(isapprox.(strain, 0.0, atol=1e-12))
        end

        # Pure translation in y-direction 
        u_translation_y = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
        
        for qp in 1:length(shape_data.weights)
            B_voigt = B_elem.voigt_gradient_operator[:,:,qp]
            strain = B_voigt * u_translation_y
            @test all(isapprox.(strain, 0.0, atol=1e-12))
        end
        
        # Small rigid body rotation 
        theta = 0.01
        u_rotation = zeros(8)
        for i in 1:4 
            x, y = nodes[i, 1], nodes[i, 2]
            u_rotation[2*i-1] = -theta * y
            u_rotation[2*i]   =  theta * x
        end

        for qp in 1:length(shape_data.weights)
            B_voigt = B_elem.voigt_gradient_operator[:,:,qp]
            strain = B_voigt * u_rotation
            @test all(isapprox.(strain, 0.0, atol=1e-12))
        end
    end
end
