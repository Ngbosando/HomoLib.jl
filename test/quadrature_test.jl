using HomoLib: integration_rule
using Test


function test_integration_rule_only()
    @testset "Integration Rules up to 3rd order" begin
        elements = [
            :Lin2, :Lin3, :Lin4,
            :Tri3, :Tri6, :Tri10, :Tri15,
            :Quad4, :Quad9, :Quad16,
            :Tet4, :Tet10, :Tet20,
            :Hex8, :Hex27, :Hex64,
            # :Pyr5, :Pyr14,
            :Pri6, :Pri18
        ]

        for elem in elements
            for order in 1:3
                qorder = get_quadrature_order(elem, order)
                @testset "$elem (order $order → q$order)" begin
                    qp, w = integration_rule(elem, qorder)

                    @test length(qp) == length(w)
                    @test all(isfinite, w)
                    @test all(p -> all(isfinite, p), qp)

                    vol = sum(w)
                    println("$vol")
                    @show reference_volume(elem)
                    @test isapprox(vol, reference_volume(elem); atol=1e-8)
                end
            end
        end
    end
end

# Required helpers
function get_quadrature_order(elem::Symbol, p::Int)
    tensor_elems = (:Lin2, :Lin3, :Lin4,
                    :Quad4, :Quad9, :Quad16,
                    :Hex8, :Hex27, :Hex64)
    dim = reference_dimension(elem)
    return elem in tensor_elems ? max(1, ceil(Int, p * dim / 2)) : p
end

function reference_dimension(elem::Symbol)::Int
    elem in (:Lin2, :Lin3, :Lin4) && return 1
    elem in (:Tri3, :Tri6, :Tri10, :Tri15,
             :Quad4, :Quad8, :Quad9, :Quad16) && return 2
    return 3
end

function reference_volume(shape::Symbol)::Float64
    Dict(
        :Lin2 => 2.0,
        :Lin3 => 2.0,
        :Lin4 => 2.0,
        :Tri3 => 0.5,
        :Tri6 => 0.5,
        :Tri10 => 0.5,
        :Tri15 => 0.5,
        :Quad4 => 4.0,
        :Quad9 => 4.0,
        :Quad16 => 4.0,
        :Tet4 => 1.0,
        :Tet10 => 1.0,
        :Tet20 => 1.0,
        :Hex8 => 8.0,
        :Hex27 => 8.0,
        :Hex64 => 8.0,
        :Pyr5 => 4/3,
        :Pyr14 => 4/3,
        :Pri6 => 1.0,
        :Pri18 => 1.0,
    )[shape]
end

# Run test
test_integration_rule_only()
