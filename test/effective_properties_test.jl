using Revise 
using Test
using HomoLib: material_def, compute_effective_property
 

# @testset "Effective Property (Thermal)" begin
    # Single triangular element
    nodes = [
        0.0 0.0;
        1.0 0.0;
        0.0 1.0
    ]
    elements = [1 2 3]
 
# using HomoLib:shape_data,build_B_matrices,jacobian_data
# Define inputs
Uᵣ₁ = rand(100)  # Example nodal solution for load case 1
Uᵣ₂ = rand(100)  # Example nodal solution for load case 2
connectivity = [1 2 3; 4 5 6]  # Example connectivity matrix
Prop = [1.0 0.0; 0.0 1.0], [2.0 0.0; 0.0 2.0]  # Material properties
type_elem = [1, 2]  # Element types
physical_group = [(1, 1)]  # Physical groups
nodes_per_element = 3
element_order = 1
nodes = rand(100, 2)  # Example nodal coordinates (2D)
element_type = :Tri3  # Element type (triangle)
    mat = material_def([:thermal], 2, :isotropic; κ=1.0)
    solver_results = (U = (zeros(3), zeros(3)),)
 
# Compute effective conductivity

    result = compute_effective_property(
        mat,
        elements,
        nodes,
        nothing,
        solver_results,
        :Tri3,
        1,
        2,
    )

    @test size(result.K) == (2, 2)

# end
