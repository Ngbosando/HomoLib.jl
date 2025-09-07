using HomoLib: ElasticMaterial, precompute_geometric_data,force_computation
using HomoLib: plaque, assemble_KMF, solve!, DirichletBC
using CairoMakie
using LinearAlgebra, SparseArrays, Statistics
using Test
using StaticArrays


function compute_energy(nx, ny, element_type)

    order = 1
    b = 100.0
    h = 20.0
    nodes, connectivity, border_tag = plaque(b, h, 0.0, nx, ny, "patch2D", order, element_type; show_gui=false)
    
    # Create material using new constructor
    mat = ElasticMaterial(2, E=1000.0, ν=0.25, stress=true)  
    
    el_shape = element_type == :triangle ? :Tri3 : :Quad4
    int_shape = element_type == :triangle ? 1 : 2

    Geometric_Data = precompute_geometric_data(
        el_shape, int_shape, mat.dim,
        connectivity, nodes, mat
    );
    
    # Define force only on right edge (parabolic shear)
     NodalForces = Dict(
    :left => (fᵥ = nothing,
          fₜ = force_computation((x,y) -> [0.0, 6 * (1 - ((y-10)/10)^2)])),

    # Right edge (x and y components)
    :right => (fᵥ = nothing,
            fₜ = force_computation((x,y) -> [-3*80*x*(y-10)/(2*10^3), -6 * (1 - ((y-10)/10)^2)]))
    )
    BoundaryFace = Dict(
        :left  => (element_border = border_tag[:left],  element_type = :Lin2, dim = 2, nodes = nodes, int_order = 1),
        :right => (element_border = border_tag[:right], element_type = :Lin2, dim = 2, nodes = nodes, int_order = 1)
    )
   
    # Assemble stiffness matrix and force vector
    K_elasticity,_, F = assemble_KMF(
        connectivity,
        nodes,
        mat,
        2,
        Geometric_Data;
        NodalForces = NodalForces,
        BoundaryFace = BoundaryFace
    );
 
    function center_node(nodes, bnodes)
        coords = nodes[bnodes, :]
        centre = mean(coords, dims=1)
        dists = [norm(c' - centre) for c in eachrow(coords)]
        return bnodes[argmin(dists)]
    end
    ce = center_node(nodes, unique(border_tag[:right]))

    bc1 = DirichletBC([ce], [true,true], [0.0,0.0], 2) 
    bc2 = DirichletBC([2,3], [true,false], [0.0,0.0], 2)

    u = zeros(size(K_elasticity,1))
    U = solve!(K_elasticity, F, u, [bc1,bc2])
    
    # Compute strain energy
    energy = 0.5 * dot(U, K_elasticity * U)
    return 2 * energy
end

# Mesh resolutions and exact energy
sizes = [(11, 5),    (21, 9),    (41, 17),   (81, 33),   (161, 65),  (321, 129)]
exact_energy = 3296.0

# Collect normalized size and normalized error
function convergence_data(elem_sym)
    h_norm = Float64[]
    eta_E = Float64[]
    b = 100.0
    base_nx = 11
    h0 = b / base_nx
    for (nx, ny) in sizes
        W = compute_energy(nx, ny, elem_sym)
        h = b / nx
        push!(h_norm, h / h0)
        push!(eta_E, abs(exact_energy - W) / exact_energy)
    end
    return h_norm, eta_E
end

h_tri, err_tri = convergence_data(:triangle)
h_quad, err_quad = convergence_data(:quadrilateral)

# Plot normalized convergence
fig = Figure(size = (800, 600))
ax = Axis(fig[1,1],
    xlabel = "Normalized element size h/h₀",
    ylabel = "Energy error η_E",
    xscale = log10, 
    yscale = log10,
    title = "Convergence in normalized energy error"
)

lines!(ax, h_tri, err_tri, label = "Tri3", linewidth = 2)
lines!(ax, h_quad, err_quad, label = "Quad4", linewidth = 2)
axislegend(ax; position = :rt)
save("normalized_energy_convergence.png", fig)
display(fig)

# Reference: O.C. Zienkiewicz,R.L. Taylor,J.Z. Zhu
# The Finite Element Method: Its Basis and Fundamentals 
@testset "Table 2.1 convergence behavior" begin
    for elem_sym in (:triangle, :quadrilateral)
        W_vals = [compute_energy(nx, ny, elem_sym) for (nx, ny) in sizes]
        @test isapprox(W_vals[end], exact_energy; rtol=1e-3)
    end
end