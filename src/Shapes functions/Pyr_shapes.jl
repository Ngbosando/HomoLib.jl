# =============================================
# Pyramid Shape Functions (3D Coordinates)
# =============================================

function shapeFunctions_PyramidN(ξ::Real, η::Real, ζ::Real, p::Int)
    # Map to unit pyramid coordinates
    x̂ = (ξ + 1) / 2
    ŷ = (η + 1) / 2
    ẑ = (ζ + 1) / 2

    # Generate base 1D nodes
    nodes = collect(range(0.0, 1.0, length = p + 1))

    # Build Lagrange basis for 1D
    function lagrange_basis(i, nodes, x)
        prod([(x - nodes[j]) / (nodes[i] - nodes[j]) for j in 1:length(nodes) if j != i])
    end
    function lagrange_deriv(i, nodes, x)
        sum([
            prod([(x - nodes[k]) / (nodes[i] - nodes[k]) for k in 1:length(nodes) if k != i && k != j]) /
            (nodes[i] - nodes[j]) for j in 1:length(nodes) if j != i
        ])
    end

    # Pyramid node coordinates in [0,1]^3
    coords = []
    for k in 0:p
        h = k / p
        nxy = p - k
        if nxy == 0
            push!(coords, (0.5, 0.5, h))
            continue
        end
        for j in 0:nxy
            for i in 0:nxy
                push!(coords, (i / nxy, j / nxy, h))
            end
        end
    end

    # Shape functions
    N = Float64[]
    dN_dξ = Float64[]
    dN_dη = Float64[]
    dN_dζ = Float64[]

    for (x, y, z) in coords
        i = findfirst(isequal(x), nodes)
        j = findfirst(isequal(y), nodes)
        k = findfirst(isequal(z), nodes)

        if i === nothing || j === nothing || k === nothing
            error("Invalid local coordinate ($x, $y, $z) for order $p")
        end

        ℓx = lagrange_basis(i, nodes, x̂)
        ℓy = lagrange_basis(j, nodes, ŷ)
        ℓz = lagrange_basis(k, nodes, ẑ)

        dℓx = lagrange_deriv(i, nodes, x̂) * 0.5
        dℓy = lagrange_deriv(j, nodes, ŷ) * 0.5
        dℓz = lagrange_deriv(k, nodes, ẑ) * 0.5

        push!(N, ℓx * ℓy * ℓz)
        push!(dN_dξ, dℓx * ℓy * ℓz)
        push!(dN_dη, ℓx * dℓy * ℓz)
        push!(dN_dζ, ℓx * ℓy * dℓz)
    end

    return N, dN_dξ, dN_dη, dN_dζ
end

function shape_functions(::Pyr5, ξ::T, η::T, ζ::T) where T<:Real
    return shapeFunctions_PyramidN(ξ, η, ζ, 1)
end

function shape_functions(::Pyr14, ξ::T, η::T, ζ::T) where T<:Real
    return shapeFunctions_PyramidN(ξ, η, ζ, 2)
end

function shape_functions(::Pyr29, ξ::T, η::T, ζ::T) where T<:Real
    return shapeFunctions_PyramidN(ξ, η, ζ, 3)
end

function shape_functions(::Pyr50, ξ::T, η::T, ζ::T) where T<:Real
    return shapeFunctions_PyramidN(ξ, η, ζ, 4)
end

function shape_functions(::Pyr77, ξ::T, η::T, ζ::T) where T<:Real
    return shapeFunctions_PyramidN(ξ, η, ζ, 5)
end
