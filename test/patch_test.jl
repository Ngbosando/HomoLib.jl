using Revise
using HomoLib: cube,material_def,assemble_global_matrix,plaque,extract_border_nodes_from_elements
using LinearAlgebra

# ---------------------------------------------------------------------
# 3D Elastic patch test with a single Hex8 element
# ---------------------------------------------------------------------

# Generate simple cube mesh with one hexahedral element
nodes3D, elems3D, bottom, back, front, left, top, right = cube("patch3D.msh",
    1, 1.0, 1.0, 1.0, 1, :Hex8; show_gui=false)

# Material definition (linear isotropic elasticity)
mat3D = material_def([:elastic], 3, :isotropic; E=1.0e9, ν=0.3)

# Assemble stiffness matrix
K3D = assemble_global_matrix(elems3D, nodes3D, :Hex8, 2, mat3D, 3, nothing)

# Exact linear displacement producing constant strain ε_xx
ε = 1.0e-3
U3D = Float64[]
for (x, y, z) in eachrow(nodes3D)
    push!(U3D, ε * x)   # u_x
    push!(U3D, 0.0)     # u_y
    push!(U3D, 0.0)     # u_z
end

# Residual of K * U for patch test
res3D = K3D * U3D
println("3D elastic Hex8 patch residual norm = ", norm(res3D))

# ---------------------------------------------------------------------
# 2D Piezoelectric patch test with Quad4 elements
# ---------------------------------------------------------------------

# Create a 1×1 quadrilateral patch (2×2 nodes)
nodes2D, elems2D, _ = plaque(1.0, 1.0, 0.1, 2, 2, "patch2D", 1, :quadrilateral; show_gui=false)

# Piezoelectric material (plane stress)
mat2D = material_def([:elastic, :electric], 2, :isotropic;
    E=1.0e9, ν=0.3, plane_stress=true,
    e=zeros(2,3), ϵ=Matrix(0.1I(2)))

# Assemble coupled stiffness matrix
K2D = assemble_global_matrix(elems2D, nodes2D, :Quad4, 2, mat2D, 2, nothing)

# Linear mechanical displacement and electric potential
α = 1.0e-3                     # mechanical strain
ϕgrad = [1.0, 0.0]             # constant electric field Ex = 1
U_mech = Float64[]
for (x, y) in eachrow(nodes2D)
    push!(U_mech, α * x)       # u_x
    push!(U_mech, 0.0)         # u_y
end
ϕ = [ϕgrad[1]*x + ϕgrad[2]*y for (x, y) in eachrow(nodes2D)]
U2D = vcat(U_mech, ϕ)

# Residual check
res2D = K2D * U2D
println("2D piezoelectric Quad4 patch residual norm = ", norm(res2D))
 