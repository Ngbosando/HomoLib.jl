using Revise
using HomoLib
using Test
function generate_plate_with_inclusion(L=1.0, R=0.5, N=5)
    # Generate node coordinates
    x_coords = range(-L, L, length=N)
    y_coords = range(-L, L, length=N)
    
    # Create nodes matrix [x, y]
    nodes = Matrix{Float64}(undef, N*N, 2)
    k = 1
    for j in 1:N
        y = y_coords[j]
        for i in 1:N
            x = x_coords[i]
            nodes[k, 1] = x
            nodes[k, 2] = y
            k += 1
        end
    end
    
    # Initialize element containers
    elements = Vector{Vector{Int}}()
    type_elem = Int[]
    Prop = Float64[]
    
    # Create triangular elements from grid
    for j in 1:N-1
        for i in 1:N-1
            ll = (j-1)*N + i      # Lower-left node index
            lr = ll + 1           # Lower-right node index
            ur = lr + N           # Upper-right node index
            ul = ll + N           # Upper-left node index
            
            # Create two triangles per square
            push!(elements, [ll, lr, ul])  # First triangle
            push!(elements, [lr, ur, ul])  # Second triangle
        end
    end
    
    # Classify elements and assign properties
    for (e, elem) in enumerate(elements)
        # Get element centroid coordinates
        n1 = nodes[elem[1], :]
        n2 = nodes[elem[2], :]
        n3 = nodes[elem[3], :]
        centroid = [(n1[1] + n2[1] + n3[1])/3, (n1[2] + n2[2] + n3[2])/3]
        
        # Calculate distance from origin
        distance = sqrt(centroid[1]^2 + centroid[2]^2)
        
        # Assign element type and properties
        if distance <= R
            push!(type_elem, 1)
            push!(Prop, 1.0)       # k1 = 1 for inclusion
        else
            push!(type_elem, 2)
            push!(Prop, 5.0)       # k2 = 5 for matrix
        end
    end
    
    # Define physical_group as a tuple of pairs
    physical_group = ((1, 1.0), (2, 5.0))  # (group_id, property_value)
    
    # Mesh parameters
    nodes_per_element = 3
    element_order = 1
    
    return (nodes = nodes,
            elements = elements,
            Prop = Prop,
            type_elem = type_elem,
            physical_group = physical_group,
            nodes_per_element = nodes_per_element,
            element_order = element_order)
end

# Example usage
mesh_data = generate_plate_with_inclusion()

 
function mat_properties(kα,kᵦ)
    # Propriétés thermiques de l'inclusion
    kₘ = [kᵦ 0
           0 kᵦ]  # conductivité thermique (isotrope)
    kᵢ =  [kα 0
            0 kα]  # conductivité thermique (isotrope)
    Prop = [kₘ, kᵢ]
    return Prop
end
# inclusion
kα = 5
# matrice 
kᵦ = 1
Prop = mat_properties(kα, kᵦ)


@testset "HomoLib.jl" begin
    K = fem_th(mesh_data[2],
                mesh_data[1],
                Prop,
                mesh_data[4],
                mesh_data[5],
                mesh_data[6],
                mesh_data[7]) 
                
    
end
