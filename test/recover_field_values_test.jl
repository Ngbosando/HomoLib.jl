using Revise
using Test
using HomoLib: recover_field_values,material_def,Material,BoundaryCondition, solve!,getBoundaryNodes,force_computation
using LinearAlgebra
using CairoMakie
using Statistics
using HomoLib: plaque, assemble_global_matrix,Transform_boundary
using BenchmarkTools

@testset "Recover Field Values - Analytical Checks" begin

    dim = 2
    coords = [0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0]
    conn = [1 2 3; 1 3 4]
    element_type = :Tri3
    order = 1
    elem_phase = [1, 1]

    # Elastic Case: Plane stress tension
    E, ν = 210e9, 0.3
    mat_el = material_def([:elastic], dim, :isotropic; E=E, ν=ν, plane_stress=true)
    u_exact = x -> [1e-3 * x[1], 0.0]  # Linear stretch in x

    nodes_u = [u_exact(coords[i, :]) for i in 1:4]
    Uvec = reduce(vcat, nodes_u)
    Uresult = (U = Uvec,)

    recovered = recover_field_values(conn, coords, mat_el, Uresult, elem_phase, element_type, order, dim)

    avg_stress = mean(recovered.stress[:, 1])  # σ₁₁
    σ_expected = E * 1e-3 

    @test isapprox(avg_stress, σ_expected; rtol=0.1)


    # Thermal Case: Gradient of T = [1, 0]
    k = 10.0
    mat_th = material_def([:thermal], dim, :isotropic; κ=k)
    T = coords[:, 1]  # x → T
    Uresult_T = (T = T,)
    recovered_th = recover_field_values(conn, coords, mat_th, Uresult_T, elem_phase, element_type, order, dim)
    avg_flux = mean(recovered_th.flux[:, 1])
    @test isapprox(avg_flux, -k; atol=1e-2)

    # Piezoelectric: Constant ε & E fields
   
    e_tensor = 1e-12 * ones(3, 2)
    ϵ = Matrix(1e-10 * I(2))
 
    mat_pz = material_def([:elastic, :electric], 2, :isotropic, E=100e9, ν=0.25, e=e_tensor,ϵ=ϵ)
    U_mech = reduce(vcat, fill([0.01, 0.0], 8))  # uniform ε₁₁
    ϕ = reduce(vcat, fill(1.0, 4))               # ∇ϕ → zero
    Uresult_pz = (U=U_mech, ϕ=ϕ)
    recovered_pz = recover_field_values(conn, coords, mat_pz, Uresult_pz, elem_phase, element_type, order, dim)

    # @test isapprox(mean(recovered_pz.stress[:, 1]), C[1, 1]*0.01; atol=1e5)
    @test norm(mean(recovered_pz.electric)) ≈ 0.0

end

@testset "1D Heat Conduction – Analytical Thermal" begin
    nodes = [0.0; 0.5; 1.0]
    coords = hcat(nodes, zeros(length(nodes)))
    elements = [1 2; 2 3]

    # Constant thermal conductivity
    k = 10.0
    mat = material_def([:thermal], 1, :isotropic, κ=k)

    # Temperature linearly varying from 100 to 100 across 1m
    T = [100.0, 150.0, 100.0]
    Uresult = (T = T,)

    recovered = recover_field_values(elements, coords[:,1], mat, Uresult, nothing, :Lin2, 1, 1)
    flux = recovered.flux[:,1]

    expected_flux = -k * (100 - 100) / 1.0
    avg_flux = mean(flux)

    @info "Recovered avg_flux = $avg_flux"
    @info "Expected flux = $expected_flux"

    @test isapprox(avg_flux, expected_flux; rtol=1e-10)
end
         
# Geometry and mesh
b, h, lc, lt1, lt2 = 100.0, 10.0, 0.1, 4,3
element_order = 5
element_type = :quadrilateral
filename = "Test_plate.msh"
nodes, elements, border_nodes = plaque(b, h, lc, lt1,lt2, filename, element_order, element_type;show_gui=true);
dim = 2;


# === Elastic Material ===

mat = material_def([:elastic], dim, :isotropic; E=1000, ν=0.25, plane_stress=false);
# Assemble system
NodalForces = Dict(
    :left => (fᵥ = nothing, fₜ = force_computation((x, y)  -> [0.0, 3*80/(4*10) * (1 - (y/10)^2)])),
    :right => (fᵥ = nothing, fₜ = force_computation((x, y) -> [-3*80*x*y/(2*10^3), -3*80/(4*10)*(1 - (y/10)^2)]))
);



BoundaryFace = Dict(
    :right => (
        element_border = border_nodes.right,
        element_type = :Lin6,
        dim = 1,
        nodes = nodes[:,2],
        int_order = 5
    ),
    :left => (
        element_border = border_nodes.left,
        element_type = :Lin6,
        dim = 1,
        nodes = nodes[:,2],
        int_order = 5
    )
);
x_vec = nodes[:,1];
y_vec = nodes[:,2];

y_bl = y_vec[unique(border_nodes.right)];
x_bl = x_vec[unique(border_nodes.right)];

# Compute stresses for boundary nodes
function beam_stresses(P, x_vec::AbstractVector, y_vec::AbstractVector, c)
    """
    Calculate stresses (σ_x, σ_y, τ_xy) for multiple nodes in a beam under a point load.

    Parameters:
    - P: Applied load (force, scalar)
    - x_vec: Vector of x-coordinates
    - y_vec: Vector of y-coordinates (same length as x_vec)
    - c: Half-depth of the beam (scalar, height = 2c)

    Returns:
    - Tuple of vectors (σ_x, σ_y, τ_xy), each of same length as x_vec
    """
    
    # Input validation
    @assert length(x_vec) == length(y_vec) "x and y vectors must have the same length"
    
    # Compute stresses (broadcast over vectors)
    σ_x = @. -(3/2) * (P * x_vec * y_vec) / c^3
    σ_y = zeros(length(x_vec))  # σ_y is always zero
    τ_xy = @. -(3P / (4c)) * (1 - (y_vec/c)^2)
    
    return (σ_x, σ_y, τ_xy)
end

σ_x, σ_y, τ_xy = beam_stresses(80, x_bl, y_bl, 10);
indd = Transform_boundary(unique(border_nodes.right));
indg = Transform_boundary(unique(border_nodes.left));
# === Stiffness assembly ===
K_elasticity,F = assemble_global_matrix(
elements,
nodes,
:Quad36,    # Element type
2*5,         # Integration order
mat,
dim,          # Dimension
nothing;
NodalForces,
BoundaryFace 
);
function check_stiffness_matrix(K)
    n = size(K, 1)
    
    # 1. Check for zero diagonals
    zero_diags = findall(i -> K[i,i] == 0, 1:n)
    if !isempty(zero_diags)
        @warn "Zero diagonals at DOFs: $zero_diags"
    end
    
    # 2. Check minimum stiffness
    non_zero_diags = [K[i,i] for i in 1:n if K[i,i] != 0]
    if isempty(non_zero_diags)
        error("All diagonals are zero!")
    end
    
    min_diag = minimum(non_zero_diags)
    max_diag = maximum(non_zero_diags)
    @info "Stiffness matrix diagonal range: $min_diag to $max_diag"
    
    # 3. Check condition number (expensive but useful)
    if n < 5000
        cond_num = cond(Matrix(K))
        @info "Condition number: $cond_num"
        if cond_num > 1e12
            @warn "Extremely ill-conditioned matrix!"
        end
    end
end

# Before solving:
check_stiffness_matrix(K_elasticity)



f_global = zeros(size(K_elasticity, 1));
mid = unique(border_nodes.right)[ceil(Int, length(unique(border_nodes.right))/2)];
dofs_per_node = 2
leftNodes,rightNodes, topNodes, bottomNodes =  getBoundaryNodes(nodes, b);
dirichlet_bc1 = BoundaryCondition(:dirichlet, [1,4], [true, false], [0.0,0.0], dofs_per_node); # x= L ,y=c & y=-c
dirichlet_bc2 = BoundaryCondition(:dirichlet, [mid], [true, true],  [0.0,0.0], dofs_per_node); # x= L ,y=0


# Solve
U = solve!(K_elasticity, F, f_global, [dirichlet_bc1,dirichlet_bc2]);
El =  U' * K_elasticity * U
ue = F'*U
function beam_displacement(P, x_vec::AbstractVector, y_vec::AbstractVector, L, c, E, I, G, ν)
    """
    Compute displacements (u, v) for multiple nodes in a beam under a point load.

    Parameters:
    - P: Applied load (force, scalar)
    - x_vec, y_vec: Vectors of x and y coordinates (same length)
    - L: Beam length (scalar)
    - c: Half-depth of beam (scalar, height = 2c)
    - E: Young's modulus (scalar)
    - I: Area moment of inertia (scalar)
    - G: Shear modulus (scalar)
    - ν: Poisson's ratio (scalar)

    Returns:
    - Tuple (u, v): Vectors of horizontal and vertical displacements
    """

    @assert length(x_vec) == length(y_vec) "x and y vectors must have the same length"

    # Precompute common terms
    EI = E * I
    GI = G * I

    # Compute u components (broadcast over vectors)
    term1_u = @. -P * (x_vec^2 - L^2) * y_vec / (2 * EI)
    term2_u = @. -ν * P * y_vec * (y_vec^2 - c^2) / (6 * EI)
    term3_u = @. P * y_vec * (y_vec^2 - c^2) / (6 * GI)
    u = term1_u + term2_u + term3_u

    # Compute v components
    term1_v = @. ν * P * x_vec * y_vec^2 / (2 * EI)
    term2_v = @. P * (x_vec^3 - L^3) / (6 * EI)
    term3_v = @. (P * L^2 / (2 * EI) + ν * P * c^2 / (6 * EI) + P * c^2 / (3 * GI)) * (x_vec - L)
    v = term1_v + term2_v - term3_v

    return (u, v)
end

# Define material and geometric properties
P = 80.0  # Load in N
L = 100.0    # Length in m
c = 10.0    # Half-depth in m
E = 1000    # Young's modulus in Pa (steel)
ν = 0.25    # Poisson's ratio
G = E/(2*(1+ν)) # Shear modulus
I = 2*c^3/3  # Moment of inertia for square cross-section

# Evaluate at a point (x,y)

u, v = beam_displacement(P, nodes[:,1], nodes[:,2], L, c, E, I, G, ν)

findall(x -> x == 0, sum(K_elasticity .!= 0, dims=2))
# Recover
magn=0.0;
Ux = U[1:2:end];
Uy = U[2:2:end];
Nx₂ = nodes[:,1] + Ux * magn;
Ny₂ = nodes[:,2] + Uy * magn ;
nodes2 = hcat(Nx₂, Ny₂);
Uresult = (U = U,)
field_vals = recover_field_values(
    elements, nodes2, mat, Uresult,
    nothing, :Quad4, 2, 2
)
σ₁₁, σ₂₂, σ₁₂ = field_vals.stress[:,1], field_vals.stress[:,2], field_vals.stress[:,3];

# Analytical strain = [1 0; 0 0], stress = C * ε
E = 2.1e11
ν = 0.3
C = (E / (1 - ν^2)) * [
    1.0  ν   0.0;
    ν    1.0 0.0;
    0.0  0.0 (1 - ν) / 2
]
ε = [1.0, 0.0, 0.0]
σ = C * ε

# Numerical average
avg_stress = mean(field_vals.stress, dims=1)
avg_strain = mean(field_vals.strain, dims=1)

# Compute L2 residuals
stress_residual = norm(avg_stress[:] - σ) / norm(σ)
strain_residual = norm(avg_strain[:] - ε) / norm(ε)

println("Stress residual: ", stress_residual)
println("Strain residual: ", strain_residual)
avg_strain = [mean(field_vals.strain[:, 1]) mean(field_vals.strain[:, 2]) mean(field_vals.strain[:, 3])]  
println("avg_strain: ", avg_strain)
println(" ε: ",  ε')
# Visualisation displacement and strain/stress
using DelaunayTriangulation
using HomoLib: plot_champ_scalaire
            magn=0.0;
            Ux = U[1:2:end];
            Uy = U[2:2:end];
            Nx₂ = nodes[:,1] + Ux * magn;
            Ny₂ = nodes[:,2] + Uy * magn ;
            nodes2 = hcat(Nx₂, Ny₂);
            # Displacement
                points = Point2f.(nodes[:,1], nodes[:,2]);
                # On fait une triangulation Delaunay sur ces points
                tri = triangulate(points);
                fig, ax, _ = triplot(
                tri;                   # Triangulation existante
                show_points    = false,
                triangle_color = :transparent,
                strokecolor    = :blue,
                strokewidth    = 0.8
                );
                n_major = 5;
                x_min, x_max = minimum(nodes[:,1]), maximum(nodes[:,1]);
                y_min, y_max = minimum(nodes[:,2]), maximum(nodes[:,2]);
                xt = range(x_min, x_max; length = n_major);
                yt = range(y_min, y_max; length = n_major);

                # labels arrondis à 1 décimale
                xtl = string.(round.(xt; digits = 1));
                ytl = string.(round.(yt; digits = 1));

                ax.xticks = (xt, xtl);
                ax.yticks = (yt, ytl);
                fig

σ₁₁, σ₂₂, σ₁₂ = field_vals.stress[:,1], field_vals.stress[:,2], field_vals.stress[:,3];
                ϵ₁₁, ϵ₂₂, ϵ₁₂ = field_vals.strain[:,1], field_vals.strain[:,2], field_vals.strain[:,3];
                σ₁₁_plot = plot_champ_scalaire( nodes2, elements, σ₁₁, :Tri6)
                σ₂₂_plot = plot_champ_scalaire( nodes2, elements, σ₂₂, :Tri6)
                σ₁₂_plot = plot_champ_scalaire( nodes2, elements, σ₁₂, :Tri6)
                ϵ₁₁_plot = plot_champ_scalaire( nodes2, elements, ϵ₁₁, :Tri6)
                ϵ₂₂_plot = plot_champ_scalaire( nodes2, elements, ϵ₂₂, :Tri6)
                ϵ₁₂_plot = plot_champ_scalaire( nodes2, elements, ϵ₁₂, :Tri6)