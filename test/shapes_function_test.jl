using HomoLib: shape_functions
using Test
using LinearAlgebra
using Printf
using Dates

# Helper function to safely sum gradients and check their properties
function validate_shape_functions(element::Symbol, coords::Vector{NTuple}, atol=1e-6)
    @info "Testing element: $element with $(length(coords)) sample points"
    t_start = now()
    for ξ in coords
        N, ∇N = shape_functions(element, ξ...)
        @test isapprox(sum(N), 1.0; atol=atol)
        ∇sum = sum(∇N; dims=1)
        @test all(isapprox.(∇sum, zero(∇sum); atol=atol))
    end
    println("  -> [$(element)] passed in $(now() - t_start)")
end

# Reference coordinates for 1D elements
lin_coords(n) = [(ξ,) for ξ in range(-1.0, 1.0; length=n)]

# Reference coordinates for 2D elements (Quad, Tri)
quad_coords(n) = [(ξ, η) for η in range(-1.0, 1.0; length=n), ξ in range(-1.0, 1.0; length=n)] |> vec
tri_coords() = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (1/3, 1/3), (0.25, 0.5)]

# Reference coordinates for 3D elements (Hex, Tet, Pyr, Pri)
hex_coords(n) = [(x, y, z) for z in range(-1.0, 1.0; length=n),
                             y in range(-1.0, 1.0; length=n),
                             x in range(-1.0, 1.0; length=n)] |> vec

tet_coords() = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0), (1/3, 1/3, 1/3)]
pyr_coords() = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0), (0.25, 0.25, 0.25)]
pri_coords() = [(0.0, 0.0, -1.0), (1.0, 0.0, -1.0), (0.0, 1.0, -1.0), (1/3, 1/3, 0.0)]

# Element mapping and their reference samples
element_tests = Dict(
    :Lin2 => lin_coords(3),
    :Lin3 => lin_coords(4),
    :Lin4 => lin_coords(5),
    :Lin5 => lin_coords(6),

    :Tri3 => tri_coords(),
    :Tri6 => tri_coords(),
    :Tri10 => tri_coords(),
    :Tri15 => tri_coords(),
    :Tri21 => tri_coords(),

    :Quad4 => quad_coords(2),
    # :Quad9 => quad_coords(3),
    :Quad16 => quad_coords(4),
    # :Quad25 => quad_coords(5),
    :Quad36 => quad_coords(6),

    :Hex8 => hex_coords(2),
    :Hex27 => hex_coords(3),
    :Hex64 => hex_coords(4),
    :Hex125 => hex_coords(5),
    :Hex216 => hex_coords(6),

    :Tet4 => tet_coords(),
    :Tet10 => tet_coords(),
    :Tet20 => tet_coords(),
    :Tet35 => tet_coords(),
    :Tet56 => tet_coords(),

    # :Pyr5 => pyr_coords(),
    # :Pyr14 => pyr_coords(),
    # :Pyr29 => pyr_coords(),

    :Pri6 => pri_coords(),
    :Pri18 => pri_coords(),
    :Pri36 => pri_coords()
    
)

@testset "Shape Function Tests" begin
    for (el, coords) in element_tests
        validate_shape_functions(el, coords)
    end
end
