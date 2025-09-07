# =============================================
# Tetraheral integration rule
# =============================================

function int_rule(::TetrahedralElement, order::Int)
    @assert order ≥ 1 && order ≤ 5 "Supported integration orders: 1–5 for Tetrahedron"

    if order == 1
        # 1-point rule 
        pts = [(1/4, 1/4, 1/4)]
        weights = [1.0]

    elseif order == 2
        # 4-point rule 
        a = 0.5854101966249685
        b = 0.1381966011250105
        pts = [
            (b, b, b),
            (a, b, b),
            (b, a, b),
            (b, b, a)
        ]
        weights = fill(0.25, 4) 

    elseif order == 3
        # 5-point rule (degree 3) 
        pts = [
            (1/4, 1/4, 1/4),
            (1/2, 1/6, 1/6),
            (1/6, 1/2, 1/6),
            (1/6, 1/6, 1/2),
            (1/6, 1/6, 1/6)
        ]
        weights = [-4/5, 9/20, 9/20, 9/20, 9/20]  

    elseif order == 4
        # 8-point rule (degree 4)
        a = 0.771642902067237
        b = 0.0761190326442543
        c = 0.315701149778202
        d = 0.108103018168070
        pts = [
            (a, b, b), (b, a, b), (b, b, a), (b, b, b),
            (c, c, d), (c, d, c), (d, c, c), (c, c, c)
        ]
        base_weights = [
            0.00665379170969446,
            0.00665379170969446,
            0.00665379170969446,
            0.00665379170969446,
            0.00167953517588677,
            0.00167953517588677,
            0.00167953517588677,
            0.00167953517588677
        ]
        total = sum(base_weights)
        scale = 1.0 / total
        weights = base_weights .* scale

    elseif order == 5
        # 15-point rule (degree 5) 
        α = 1/4
        β1 = (7 - sqrt(15)) / 34
        β2 = (7 + sqrt(15)) / 34
        γ1 = (13 + 3*sqrt(15)) / 34
        γ2 = (13 - 3*sqrt(15)) / 34
        δ = (5 - sqrt(15)) / 20
        ε = (5 + sqrt(15)) / 20
        pts = [
            (α, α, α),
            (β1, β1, β1),
            (β1, β1, γ1),
            (β1, γ1, β1),
            (γ1, β1, β1),
            (β2, β2, β2),
            (β2, β2, γ2),
            (β2, γ2, β2),
            (γ2, β2, β2),
            (δ, δ, ε),
            (δ, ε, δ),
            (ε, δ, δ),
            (δ, ε, ε),
            (ε, δ, ε),
            (ε, ε, δ)
        ]
        base_weights = [
            0.1884185567365411,
            0.0670385837260428,
            0.0670385837260428,
            0.0670385837260428,
            0.0670385837260428,
            0.0452855923632739,
            0.0452855923632739,
            0.0452855923632739,
            0.0452855923632739,
            0.0330588093338715,
            0.0330588093338715,
            0.0330588093338715,
            0.0250000000000000,
            0.0250000000000000,
            0.0250000000000000
        ]
        total = sum(base_weights)
        scale = 1.0 / total
        weights = base_weights .* scale

    else
        error("Tetrahedron integration rule for order $order not implemented.")
    end

    return pts, weights
end
