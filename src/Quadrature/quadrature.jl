
include("Pyr_quadrature.jl")
include("Pri_quadrature.jl")
include("Legendre_quadrature.jl")
include("Tri_quadrature.jl")
include("Tet_quadrature.jl")


# =============================================
# General elements integration rule 
# =============================================

function integration_rule(element_type::Symbol, order::Int)
    if element_type == :Tri3
        integration_rule(Tri3(), order)
    elseif element_type == :Tri6
        integration_rule(Tri6(), order)
    elseif element_type == :Tri10
        integration_rule(Tri10(), order)
    elseif element_type == :Tri15
        integration_rule(Tri15(), order)
    elseif element_type == :Tri21
        integration_rule(Tri21(), order)
    elseif element_type == :Hex8
        integration_rule(Hex8(), order)
    elseif element_type == :Hex27
        integration_rule(Hex27(), order)
    elseif element_type == :Hex64
        integration_rule(Hex64(), order)
    elseif element_type == :Hex125
        integration_rule(Hex125(), order)
    elseif element_type == :Hex216
        integration_rule(Hex216(), order)
    elseif element_type == :Tet4
        integration_rule(Tet4(), order)
    elseif element_type == :Tet10
        integration_rule(Tet10(), order)
    elseif element_type == :Tet20
        integration_rule(Tet20(), order)
    elseif element_type == :Tet35
        integration_rule(Tet35(), order)
    elseif element_type == :Tet56
        integration_rule(Tet56(), order)
    elseif element_type == :Quad4
        integration_rule(Quad4(), order)
    elseif element_type == :Quad8
        integration_rule(Quad8(), order)
    elseif element_type == :Quad9
        integration_rule(Quad9(), order)
    elseif element_type == :Quad16
        integration_rule(Quad16(), order)
    elseif element_type == :Quad25
        integration_rule(Quad25(), order)
    elseif element_type == :Quad36
        integration_rule(Quad36(), order)
    elseif element_type == :Pyr5
        integration_rule(Pyr5(), order)
    elseif element_type == :Pyr14
        integration_rule(Pyr14(), order)
    elseif element_type == :Lin2
        integration_rule(Lin2(), order)
    elseif element_type == :Lin3
        integration_rule(Lin3(), order)
    elseif element_type == :Lin4
        integration_rule(Lin4(), order)
    elseif element_type == :Lin5
        integration_rule(Lin5(), order)
    elseif element_type == :Lin6
        integration_rule(Lin6(), order)
    elseif element_type == :Pri6
        integration_rule(Pri6(), order)
    elseif element_type == :Pri18
        integration_rule(Pri18(), order)
    elseif element_type == :Pri36
        integration_rule(Pri36(), order)
    elseif element_type == :Pri56
        integration_rule(Pri56(), order)
    elseif element_type == :Pri78
        integration_rule(Pri78(), order)
    els
        error("""
        Unsupported element type: $element_type
        Valid options: 
        - Triangular: :Tri3, :Tri6, :Tri10, :Tri15, :Tri21
        - Hexahedral: :Hex8, :Hex27, :Hex64, :Hex125, :Hex216
        - Tetrahedral: :Tet4, :Tet10, :Tet20, :Tet35, :Tet56
        - Quadrilateral: :Quad4, :Quad9, :Quad16, :Quad25, :Quad36
        - Pyramid: :Pyr5, :Pyr14, :Pyr29, :Pyr50, :Pyr77
        - Prism: :Pri6, :Pri18, :Pri36, :Pri56, :Pri78
        - Line: :Lin2, :Lin3, :Lin4, :Lin5, :Lin6
        """)
    end
end