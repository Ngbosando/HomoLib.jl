using Revise 

using HomoLib: getBoundaryNodes,plaque
using SparseArrays, Statistics, CairoMakie, Test


using HomoLib: generate_irregular_quad_patch_mesh
# Define mesh parameters
b = 1;
h = 1;
lc = 0.1; # mesh size
lt = 5;
filename = "irregular_patch.msh";
element_order = 1;
element_type = :quadrilateral;
perturb=0.3

nodes, elements, border_nodes = generate_irregular_quad_patch_mesh(element_order,:quadrilateral; nx=2, ny=2, perturb, filename);
using HomoLib: getBoundaryNodes,plaque
nodes1, elements1, border_nodes1= plaque(b, h, lc,lt,lt, filename, element_order, element_type;);
# Define exact solutions for displacement and stress
function u_exact(x, y)
    E = 210e9
    ν = 0.3
    σ_xx = 210e9
    u_x = (σ_xx / E) * x
    u_y = -ν * (σ_xx / E) * y
    return [u_x, u_y]
end

function stress_exact()
    σ_xx = 210e9
    σ_yy = 0.0
    σ_xy = 0.0
    return [σ_xx, σ_yy, σ_xy]
end


# Define boundary conditions and forces
using HomoLib: force_computation
NodalForces = Dict(
    :right => (fᵥ = nothing, fₜ = force_computation( [210e9, 0])),
    :left => (fᵥ = nothing, fₜ = force_computation( [-210e9, 0]))
)

BoundaryFace = Dict(
    :right => (
        element_border = border_nodes1.right,
        element_type = :Lin2,
        dim = 1,
        nodes = nodes1[:,2],
        int_order = 1
    ),
    :left => (
        element_border = border_nodes1.left,
        element_type = :Lin2,
        dim = 1,
        nodes = nodes1[:,2],
        int_order = 1
    )
)

using HomoLib: material_def
elastic_mat = material_def([:elastic], 2, :isotropic, 
                                E=210e9, ν=0.3, plane_stress=true);

# Assemble global stiffness matrix and force vector
using HomoLib: assemble_global_matrix
K_elasticity = assemble_global_matrix(
    elements1,
    nodes1,
    :Quad4,
    2,
    elastic_mat,
    2,
    nothing;
  
);
@show(unique(K_elasticity-K_elasticity')');
# Verification checks
function check_stiffness_matrix(K)
    n = size(K, 1)
    zero_diags = findall(i -> K[i,i] == 0, 1:n)
    !isempty(zero_diags) && @warn "Zero diagonals at DOFs: $zero_diags"
    
    non_zero_diags = [K[i,i] for i in 1:n if K[i,i] != 0]
    isempty(non_zero_diags) && error("All diagonals are zero!")
    
    min_diag, max_diag = extrema(non_zero_diags)
    @info "Stiffness matrix diagonal range: $min_diag to $max_diag"
    
    if n < 5000
        cond_num = cond(Matrix(K))
        @info "Condition number: $cond_num"
        cond_num > 1e12 && @warn "Extremely ill-conditioned matrix!"
    end
end
isapprox(K_elasticity, K_elasticity', atol=1e-12)

@test size(K_elasticity, 1) == size(K_elasticity, 2)
@test issymmetric(Matrix(K_elasticity))
@test norm(K_elasticity) > 1e-12
@test all(isfinite, K_elasticity)
# Test matrix properties
@test issymmetric(K_elasticity) 
@test all(diag(K_elasticity) .>= 0)

# Test expected matrix size
n_nodes = size(nodes, 1)
expected_size = n_nodes * 2
@test size(K_elasticity) == (expected_size, expected_size)

# Rigid body motion checks
Δx, Δy = 1.0, 1.0
u_translation = repeat([Δx, Δy], n_nodes)
θ = 0.05
uθ  = Float64[];
    for (x, y) in eachrow(nodes)
        push!(uθ , -θ*y)
        push!(uθ ,  θ*x)
    end
    u_rotation  = collect(uθ ) ; # make i

@test norm(K_elasticity * u_translation) < 1e-12
@test norm(K_elasticity * u_rotation) < 1e-12

# Apply boundary conditions
f_global = zeros(size(K_elasticity, 1))
leftBoundaryNodes, rightBoundaryNodes, top, bot = getBoundaryNodes(nodes, b)

# Fix left boundary in x, bottom boundary in y
dirichlet_bc_left = BoundaryCondition(
    :dirichlet,
    leftBoundaryNodes,
    [true, false],  # Constrain x-direction
    [0.0, 0.0],
    2
)

dirichlet_bc_bottom = BoundaryCondition(
    :dirichlet,
    bot,
    [false, true],  # Constrain y-direction
    [0.0, 0.0],
    2
)

# Solve system
U = solve!(K_elasticity, F, f_global, [dirichlet_bc_left, dirichlet_bc_bottom])
Uresult = (U = U,)

# Recover stress values
field_vals = recover_field_values(
    elements, nodes, elastic_mat, Uresult,
    nothing, :Quad4, 2, 2
)
σ₁₁, σ₂₂, σ₁₂ = field_vals.stress[:,1], field_vals.stress[:,2], field_vals.stress[:,3]

# =============================================
# Patch Test Verification
# =============================================
exact_stress = stress_exact()

# Check stress values
@testset "Stress Field Verification" begin
    for i in eachindex(σ₁₁)
        @test isapprox(σ₁₁[i], exact_stress[1], rtol=1e-4)
        @test isapprox(σ₂₂[i], exact_stress[2], rtol=1e-4)
        @test isapprox(σ₁₂[i], exact_stress[3], rtol=1e-4)
    end
end

# Check displacement solutions
@testset "Displacement Field Verification" begin
    for (i, (x, y)) in enumerate(eachrow(nodes))
        exact_disp = u_exact(x, y)
        node_dof_x = 2i - 1
        node_dof_y = 2i
        
        @test isapprox(U[node_dof_x], exact_disp[1], rtol=1e-4)
        @test isapprox(U[node_dof_y], exact_disp[2], rtol=1e-4)
    end
end

# Check residual forces
@testset "Residual Force Check" begin
    residual = K_elasticity * U - F
    @test norm(residual) < 1e-8
end