using HomoLib
using CairoMakie


# Géométrie + Matériau

mutable struct MeshData
    nodes::Matrix{Float64}
    elements::Matrix{Int}
    type_elem::Vector{Int}
    Nₓ::Vector{Float64}
    Nᵧ::Vector{Float64}
    boundary::Vector{Int}
    boundary_element
    boundary_inclusion::Vector{Int}
    master::Vector{Int}
    slave::Vector{Int}
end

mutable struct ElemData
    type::Symbol
    int_order::Int
    order::Int
end

function setup_mesh(; width, height, volume_fraction,   
                    n_inclusions, Elem::ElemData,
                    node_divisions, shape, output_file, voids, rdn, show_gui, to_rotate)  # Renamed activate -> voids

    ind_G, ind_D, ind_B, ind_H, ind_C,
    elements, Nₓ, Nᵧ, type_elem, _, boundary_element =
        generate_transfinite_plate_with_inclusions(
            width, height, volume_fraction,
            output_file, shape, n_inclusions,
            Elem.type, Elem.order,
            node_divisions[1], node_divisions[2];
            voids=voids,  # Pass voids parameter directly
            show_gui=show_gui,
            rdn=rdn,  # Disable randomization for reproducibility
            to_rotate=to_rotate
        )
  
    boundary_element = boundary_element
    boundary = unique(sort(vcat(ind_G, ind_D, ind_B, ind_H)))
    boundary_inclusion = ind_C
    master   = vcat(ind_G, ind_B)
    slave    = vcat(ind_D, ind_H)
    nodes    = hcat(Nₓ, Nᵧ)

    MeshData(nodes, elements, type_elem, Nₓ, Nᵧ, boundary, boundary_element, boundary_inclusion, master, slave)
end

Elem = ElemData(:Tri6, 2, 2)
mesh = setup_mesh(;
    width=1.0, height=1.0,
    volume_fraction=0.01, n_inclusions=1,
    Elem, node_divisions=(90, 90), shape=:circle,
    output_file="kirsch.msh",
    voids=true,
    rdn=false,
    show_gui=false,
    to_rotate=false
)

mat = ElasticMaterial(2, E=10.0, ν=1/3, stress=true)
dim = 2
GeoData = precompute_geometric_data(
    Elem.type, Elem.int_order, dim,
    mesh.elements, mesh.nodes, mat
)


# Conditions aux limites

border_elem = mesh.boundary_element
bc_left  = DirichletBC(unique(border_elem[:left]),  [true,false], [0.0,0.0], 2)
bc_right = DirichletBC(unique(border_elem[:right]), [true,false], [0.01,0.0], 2) # ux=0.01
bc_bottom= DirichletBC(unique(border_elem[:bottom]),[false,true],[0.0,0.0], 2)

K, _, _ = assemble_KMF(mesh.elements, mesh.nodes, mat, 2, GeoData)
F = zeros(size(K,1))
u = solve!(K, F, F, [bc_left, bc_right, bc_bottom])


# Récupération contraintes

res = recover_field_values(mat, mesh.elements, mesh.nodes, (U=u,), GeoData)
σxx, σyy, σxy = res.stress[:,1], res.stress[:,2], res.stress[:,3]


# Extraction bord du trou

idx = mesh.boundary_inclusion
xc, yc = 0.5, 0.5
xr = mesh.nodes[idx,1] .- xc
yr = mesh.nodes[idx,2] .- yc

θ = atan.(yr, xr)
θ .= map(x -> x < 0 ? x+2π : x, θ)

sinθ, cosθ = sin.(θ), cos.(θ)
σθθ_fem = similar(θ)
for i in eachindex(θ)
    s, c = sinθ[i], cosθ[i]
    σθθ_fem[i] = σxx[idx[i]]*s^2 + σyy[idx[i]]*c^2 - 2*σxy[idx[i]]*s*c
end

# Calcul des erreurs

min_fem, max_fem = minimum(y_fem), maximum(y_fem)

println("===== Résultats Kirsch FEM vs Analytique =====")
println("Min(FEM_norm)   ≈ ", min_fem, "  (théorique ≈ -1)")
println("Max(FEM_norm)   ≈ ", max_fem, "  (théorique ≈ +3)")
println("Sigma_inf (LS)  ≈ ", σ∞̂)

# Normalisation vs Kirsch

g = 1 .- 2 .* cos.(2 .* θ)
σ∞̂ = sum(σθθ_fem .* g) / sum(g .* g)

σθθ_fem_norm = σθθ_fem ./ σ∞̂
σθθ_kirsch_norm = g

perm = sortperm(θ)
θ_deg = θ[perm] .* 180/pi
y_fem = σθθ_fem_norm[perm]
y_ref = σθθ_kirsch_norm[perm]


# Tracé

fig = Figure(resolution=(1000,400))
ax = Axis(fig[1,1], xlabel="θ (deg)", ylabel="σθθ/σ∞",
          title="Trou circulaire sous traction (Kirsch)")
lines!(ax, θ_deg, y_ref, color=:blue, label="Kirsch")
lines!(ax, θ_deg, y_fem, color=:red, linestyle=:dash, label="FEM")
axislegend(ax, position=:rb)
save("kirsch_compare.png", fig)
