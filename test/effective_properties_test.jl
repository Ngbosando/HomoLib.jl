using Revise 

using HomoLib : compute_effective_properties

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

# Compute effective conductivity
k_eff = effective_properties_conductivity(:UTG, Uᵣ₁, Uᵣ₂, connectivity, Prop, type_elem, physical_group,
                                         nodes_per_element, element_order, nodes, element_type)