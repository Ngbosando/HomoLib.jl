include("Pyr_quadrature.jl")
include("Pri_quadrature.jl")
include("Legendre_quadrature.jl")
include("Tri_quadrature.jl")
include("Tet_quadrature.jl")


# =============================================
# General elements integration rule 
# =============================================
"""
    HomoLib.NumericalIntegration

    Module providing Gaussian quadrature rules for finite element analysis,
    implementing optimal integration schemes for various element types and orders.

    # Supported Element Families
    - Triangular elements (2D):
    - Linear (Tri3), Quadratic (Tri6), Cubic (Tri10), Quartic (Tri15), Quintic (Tri21)
    - Quadrilateral elements (2D):
    - Bilinear (Quad4), Serendipity (Quad8/Quad9), Lagrangian (Quad16, Quad25, Quad36)
    - Tetrahedral elements (3D):
    - Linear (Tet4), Quadratic (Tet10), Cubic (Tet20), Quartic (Tet35), Quintic (Tet56)
    - Hexahedral elements (3D):
    - Trilinear (Hex8), Quadratic (Hex27), Cubic (Hex64), Quartic (Hex125), Quintic (Hex216)
    - Prismatic elements (3D):
    - Linear (Pri6), Quadratic (Pri18), Cubic (Pri36), Quartic (Pri56), Quintic (Pri78)
    - Pyramid elements (3D):
    - Linear (Pyr5), Quadratic (Pyr14), Cubic (Pyr29), Quartic (Pyr50), Quintic (Pyr77)
    - Line elements (1D):
    - 2-6 point Gauss-Legendre rules (Lin2-Lin6)

    # Integration Rules
    For each element type, provides:
    - Optimal Gauss point locations in reference coordinates
    - Corresponding weights for accurate numerical integration
    - Rules supporting up to quintic polynomial exactness

    # Usage
    ```julia
    # Get integration rule for quadratic triangle
    points, weights = integration_rule(:Tri6, 3)

    # Get rule for cubic hexahedron
    points, weights = integration_rule(:Hex64, 4)
    Implementation Details
    Triangular rules based on Dunavant (1985) and Wandzura (2003)
    Tetrahedral rules from Keast (1986) and Zhang (2009)
    Hexahedral rules use tensor-product Gauss-Legendre
    Prismatic rules combine triangular and 1D Gauss rules

    References
    Dunavant, D.A. (1985). High degree efficient symmetrical Gaussian quadrature rules
    for the triangle. Int. J. Numer. Meth. Eng.
    Keast, P. (1986). Moderate-degree tetrahedral quadrature formulas.
    Comput. Methods Appl. Mech. Engrg.
"""

function lookup_element(sym::Symbol)
elem = get(ELEMENT_SYMBOL_MAP, sym, nothing)
elem === nothing && error("Unsupported element type: $sym")
return elem
end

eldim(::LinearElement)        = Val(1)
eldim(::QuadrilateralElement) = Val(2)
eldim(::TriangularElement)    = Val(2)
eldim(::HexahedralElement)    = Val(3)
eldim(::PrismaticElement)     = Val(3)
eldim(::PyramidElement)       = Val(3)
eldim(::TetrahedralElement)   = Val(3)


#  (Vector{NTuple{D,Float64}}, Vector{Float64})
function integration_rule(element_type::Symbol, order::Int, ::Val{D}) where {D}
    elem = lookup_element(element_type)
    pts_any, w_any = int_rule(elem, order)  # your per-element methods

    n   = length(w_any)
    pts = Vector{NTuple{D,Float64}}(undef, n)
    w   = Vector{Float64}(undef, n)

    @inbounds for i in 1:n
        w[i] = Float64(w_any[i])
        p    = pts_any[i]                     # Tuple/SVector/Vector, length D
        pts[i] = ntuple(j -> Float64(p[j]), D)
    end
    return pts, w
end
