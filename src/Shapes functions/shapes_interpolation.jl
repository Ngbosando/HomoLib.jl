include("Hex_shapes.jl")
include("Pyr_shapes.jl")
include("Pri_shapes.jl")
include("Quad_shapes.jl")
include("Tri_shapes.jl")
include("Lin_shapes.jl")
include("Tet_shapes.jl")

"""
    HomoLib.ShapeFunctions

    Module providing shape functions and their derivatives for various finite element types,
    implementing standard and hierarchical basis functions for accurate finite element analysis.

    # Supported Element Families
    - "Triangular Elements":
    - Linear (Tri3), Quadratic (Tri6), Cubic (Tri10) 
    - Quartic (Tri15), Quintic (Tri21)
    - "Quadrilateral Elements":
    - Bilinear (Quad4), Serendipity (Quad8/Quad9)
    - Lagrangian (Quad16, Quad25, Quad36)
    - "Hexahedral Elements":
    - Trilinear (Hex8), Quadratic (Hex27)
    - Cubic (Hex64), Quartic (Hex125), Quintic (Hex216)
    - "Tetrahedral Elements":
    - Linear (Tet4), Quadratic (Tet10)
    - Cubic (Tet20), Quartic (Tet35), Quintic (Tet56)
    - "Pyramid Elements":
    - Linear (Pyr5), Quadratic (Pyr14)
    - Cubic (Pyr29), Quartic (Pyr50), Quintic (Pyr77)
    - "Prismatic Elements":
    - Linear (Pri6), Quadratic (Pri18)
    - Cubic (Pri36), Quartic (Pri56), Quintic (Pri78)
    - "Line Elements":
    - 1-6 node versions (Lin1-Lin6)

    # Provided Quantities
    For each element type:
    - Shape functions N(ξ,η,ζ)
    - Derivatives ∂N/∂ξ, ∂N/∂η, ∂N/∂ζ
    - Consistent node ordering conventions
    - Support for isoparametric mapping

    # Usage
    ```julia
    # Get shape functions for quadratic triangle at point (0.2, 0.3)
    N, dNξ, dNη = shape_functions(:Tri6, 0.2, 0.3)

    # Get shape functions for trilinear hex at point (0.1, 0.2, 0.3)
    N, dNξ, dNη, dNζ = shape_functions(:Hex8, 0.1, 0.2, 0.3)
    Key Features
    Consistent Evaluation: Returns (N, dNξ, dNη[, dNζ]) tuple for all elements

    Optimal Basis Functions:

    Lagrangian polynomials for tensor-product elements

    Hierarchical bases for simplex elements

    Reference Coordinates:

    Triangles: Area coordinates (ξ,η)

    Quads/Hexes: Natural coordinates [-1,1]^d

    Tets: Volume coordinates (ξ,η,ζ)

    Efficient Computation: Closed-form expressions used where available

    Theory
    Shape functions satisfy:

    Partition of unity: ∑Nᵢ = 1

    Kronecker delta property: Nᵢ(ξⱼ) = δᵢⱼ

    Polynomial completeness: Can represent all polynomials up to element order

    References
    Zienkiewicz, O.C. & Taylor, R.L. (2005). "The Finite Element Method"
"""
const ELEMENT_SYMBOL_MAP = Dict{Symbol, AbstractElementType}(
    # Triangles
    :Tri3  => Tri3(),
    :Tri6  => Tri6(),
    :Tri10 => Tri10(),
    :Tri15 => Tri15(),
    :Tri21 => Tri21(),
    # Hexahedra
    :Hex8   => Hex8(),
    :Hex27  => Hex27(),
    :Hex64  => Hex64(),
    :Hex125 => Hex125(),
    :Hex216 => Hex216(),
    # Tetrahedra
    :Tet4  => Tet4(),
    :Tet10 => Tet10(),
    :Tet20 => Tet20(),
    :Tet35 => Tet35(),
    :Tet56 => Tet56(),
    # Quadrilaterals
    :Quad4  => Quad4(),
    :Quad9  => Quad9(),
    :Quad16 => Quad16(),
    :Quad25 => Quad25(),
    :Quad36 => Quad36(),
    # Pyramids
    :Pyr5  => Pyr5(),
    :Pyr14 => Pyr14(),
    :Pyr29 => Pyr29(), 
    :Pyr50 => Pyr50(), 
    # Lines
    :Lin2 => Lin2(),
    :Lin3 => Lin3(),
    :Lin4 => Lin4(),
    :Lin5 => Lin5(),
    :Lin6 => Lin6(),
    # Prisms
    :Pri6  => Pri6(),
    :Pri18 => Pri18(),
    :Pri36 => Pri36(), 
    :Pri56 => Pri56(), 
    :Pri78 => Pri78()  
)
function shape_functions(element_type::Symbol, args...)
    # Look up the element type instance from the map.
    # The `get` function provides a way to handle cases where the symbol is not found.
    element_instance = get(ELEMENT_SYMBOL_MAP, element_type, nothing)

    if element_instance === nothing
        error("""
        Unsupported element type: $element_type
        Please ensure it is one of the supported types and is defined in the ELEMENT_SYMBOL_MAP.
        """)
    end

    # This call is now type-stable because Julia's dispatch system takes over.
    return shape_functions(element_instance, args...)
end