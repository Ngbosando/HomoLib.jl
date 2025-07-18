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
function shape_functions(element_type::Symbol, args...)
    if element_type == :Tri3
        shape_functions(Tri3(), args...)
    elseif element_type == :Tri6
        shape_functions(Tri6(), args...)
    elseif element_type == :Tri10
        shape_functions(Tri10(), args...)
    elseif element_type == :Tri15
        shape_functions(Tri15(), args...)
    elseif element_type == :Tri21
        shape_functions(Tri21(), args...)
    elseif element_type == :Hex8
        shape_functions(Hex8(), args...)
    elseif element_type == :Hex27
        shape_functions(Hex27(), args...)
    elseif element_type == :Hex64
        shape_functions(Hex64(), args...)
    elseif element_type == :Hex125
        shape_functions(Hex125(), args...)
    elseif element_type == :Hex216
        shape_functions(Hex216(), args...)    
    elseif element_type == :Tet4
        shape_functions(Tet4(), args...)
    elseif element_type == :Tet10
        shape_functions(Tet10(), args...)
    elseif element_type == :Tet20
        shape_functions(Tet20(), args...)
    elseif element_type == :Tet35
        shape_functions(Tet35(), args...)
    elseif element_type == :Tet56
        shape_functions(Tet56(), args...)
    elseif element_type == :Quad4
        shape_functions(Quad4(), args...)
    elseif element_type == :Quad8
        shape_functions(Quad8(), args...)
    elseif element_type == :Quad9
        shape_functions(Quad9(), args...)
    elseif element_type == :Quad16
        shape_functions(Quad16(), args...)
    elseif element_type == :Quad25
        shape_functions(Quad25(), args...)
    elseif element_type == :Quad36        
        shape_functions(Quad36(), args...)        
    elseif element_type == :Pyr5
        shape_functions(Pyr5(), args...)
    elseif element_type == :Pyr14
        shape_functions(Pyr14(), args...)
    elseif element_type == :Pyr29
        shape_functions(Pyr29(), args...)
    elseif element_type == :Pyr50
        shape_functions(Pyr50(), args...)
    elseif element_type == :Pyr77
        shape_functions(Pyr77(), args...)
    elseif element_type == :Lin1
        shape_functions(Lin1(), args...)
    elseif element_type == :Lin2
        shape_functions(Lin2(), args...)
    elseif element_type == :Lin3
        shape_functions(Lin3(), args...)
    elseif element_type == :Lin4
        shape_functions(Lin4(), args...)
    elseif element_type == :Lin5
        shape_functions(Lin5(), args...)
    elseif element_type == :Lin6
        shape_functions(Lin6(), args...)
    elseif element_type == :Pri6
        shape_functions(Pri6(), args...)
    elseif element_type == :Pri18
        shape_functions(Pri18(), args...)
    elseif element_type == :Pri36
        shape_functions(Pri36(), args...)
    elseif element_type == :Pri56
        shape_functions(Pri56(), args...)
    elseif element_type == :Pri78
        shape_functions(Pri78(), args...)
    else
        error("""
        Unsupported element type: $element_type
        Valid options: 
        - Triangular: :Tri3, :Tri6, :Tri10, :Tri15, :Tri21
        - Hexahedral: :Hex8, :Hex27, :Hex64, :Hex125, :Hex216
        - Tetrahedral: :Tet4, :Tet10, :Tet20, :Tet35, :Tet56
        - Quadrilateral: :Quad4, :Quad9, :Quad16, :Quad25, :Quad36
        - Pyramid: :Pyr5, :Pyr14, :Pyr29, :Pyr50, :Pyr77
        - Prism: :Pri6, :Pri18, :Pri36, :Pri56, :Pri78
        - Linear: :Lin2, :Lin3, :Lin4, :Lin5, :Lin6
        """)
    end
end
  