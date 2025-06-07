

"""
    subdivide_element(elements, element_type)

Subdivide elements based on their type using predefined subdivision patterns.
Returns a vector of element vectors for consistent downstream processing.
"""
function subdivide_element(elements, element_type)
    if element_type == :Tri3
        return subdivide_node(Tri3(), elements)
    elseif element_type == :Tri6
        return subdivide_node(Tri6(), elements)
    elseif element_type == :Tri10
        return subdivide_node(Tri10(), elements)
    elseif element_type == :Tri15
        return subdivide_node(Tri15(), elements)
    elseif element_type == :Tri21
        return subdivide_node(Tri21(), elements)
    elseif element_type == :Quad4
        return subdivide_node(Quad4(), elements)
    elseif element_type == :Quad9
        return subdivide_node(Quad9(), elements)
    elseif element_type == :Quad16
        return subdivide_node(Quad16(), elements)
    elseif element_type == :Quad25
        return subdivide_node(Quad25(), elements)
    elseif element_type == :Quad36
        return subdivide_node(Quad36(), elements)
    else
        error("Unsupported element type: $element_type")
    end
end

# Triangular element subdivisions
function subdivide_node(::Tri3, elements)
    return [collect(row) for row in eachrow(elements)]
end

# Triangular element subdivisions
function subdivide_node(::Tri6,elements)
    n_elements = size(elements, 1)
    subdivided = Vector{Vector{Int}}(undef, 4*n_elements)
    for (i, element) in enumerate(eachrow(elements))
        n = element
        subdivided[4i-3] = [n[1], n[4], n[6]]
        subdivided[4i-2] = [n[4], n[2], n[5]]
        subdivided[4i-1] = [n[6], n[5], n[3]]
        subdivided[4i]   = [n[4], n[5], n[6]]
    end
    return subdivided
end

function subdivide_node(::Tri10,elements)
    n_elements = size(elements, 1)
    subdivided = Vector{Vector{Int}}(undef, 9*n_elements)
    for (i, element) in enumerate(eachrow(elements))
        n = element
        offset = 9i-8
        subdivided[offset+0] = [n[1], n[4], n[9]]
        subdivided[offset+1] = [n[4], n[9], n[10]]
        subdivided[offset+2] = [n[4], n[5], n[10]]
        subdivided[offset+3] = [n[5], n[10], n[6]]
        subdivided[offset+4] = [n[5], n[2], n[6]]
        subdivided[offset+5] = [n[9], n[10], n[8]]
        subdivided[offset+6] = [n[10], n[8], n[7]]
        subdivided[offset+7] = [n[10], n[7], n[6]]
        subdivided[offset+8] = [n[8], n[7], n[3]]
    end
    return subdivided
end

function subdivide_node(::Tri15,elements)
    n_elements = size(elements, 1)
    subdivided = Vector{Vector{Int}}(undef, 16*n_elements)
    for (i, element) in enumerate(eachrow(elements))
        n = element
        offset = 16i-15
        subdivided[offset+0]  = [n[1], n[4], n[12]]
        subdivided[offset+1]  = [n[4], n[12], n[13]]
        subdivided[offset+2]  = [n[4], n[5], n[13]]
        subdivided[offset+3]  = [n[5], n[13], n[14]]
        subdivided[offset+4]  = [n[5], n[6], n[14]]
        subdivided[offset+5]  = [n[6], n[14], n[7]]
        subdivided[offset+6]  = [n[6], n[2], n[7]]
        subdivided[offset+7]  = [n[12], n[13], n[11]]
        subdivided[offset+8]  = [n[13], n[11], n[15]]
        subdivided[offset+9]  = [n[13], n[14], n[15]]
        subdivided[offset+10] = [n[14], n[15], n[8]]
        subdivided[offset+11] = [n[14], n[7], n[8]]
        subdivided[offset+12] = [n[11], n[15], n[10]]
        subdivided[offset+13] = [n[15], n[10], n[9]]
        subdivided[offset+14] = [n[15], n[9], n[8]]
        subdivided[offset+15] = [n[10], n[9], n[3]]
    end
    return subdivided
end

function subdivide_node(::Tri21, elements)
    n_elements = size(elements, 1)
    subdivided = Vector{Vector{Int}}(undef, 25*n_elements)
    
    for (i, element) in enumerate(eachrow(elements))
        n = element
        offset = 25i - 24
        
        # First layer
        subdivided[offset+0]  = [n[1], n[4], n[19]]
        subdivided[offset+1]  = [n[4], n[5], n[20]]
        subdivided[offset+2]  = [n[5], n[6], n[21]]
        subdivided[offset+3]  = [n[6], n[2], n[7]]
        subdivided[offset+4]  = [n[19], n[20], n[22]]
        
        # Second layer
        subdivided[offset+5]  = [n[20], n[21], n[23]]
        subdivided[offset+6]  = [n[21], n[7], n[8]]
        subdivided[offset+7]  = [n[22], n[23], n[24]]
        subdivided[offset+8]  = [n[23], n[8], n[9]]
        subdivided[offset+9]  = [n[24], n[9], n[10]]
        
        # Third layer
        subdivided[offset+10] = [n[19], n[22], n[25]]
        subdivided[offset+11] = [n[22], n[24], n[26]]
        subdivided[offset+12] = [n[24], n[10], n[11]]
        subdivided[offset+13] = [n[25], n[26], n[27]]
        subdivided[offset+14] = [n[26], n[11], n[12]]
        
        # Fourth layer
        subdivided[offset+15] = [n[27], n[12], n[13]]
        subdivided[offset+16] = [n[25], n[27], n[28]]
        subdivided[offset+17] = [n[27], n[13], n[14]]
        subdivided[offset+18] = [n[28], n[14], n[15]]
        subdivided[offset+19] = [n[28], n[15], n[16]]
        
        # Fifth layer
        subdivided[offset+20] = [n[25], n[28], n[17]]
        subdivided[offset+21] = [n[28], n[16], n[18]]
        subdivided[offset+22] = [n[17], n[18], n[3]]
        subdivided[offset+23] = [n[18], n[3], n[1]]
        subdivided[offset+24] = [n[17], n[1], n[19]]
    end
    
    return subdivided
end

# Quadrilateral element subdivisions
function subdivide_node(::Quad4, elements)
   return [collect(row) for row in eachrow(elements)]
end
function subdivide_node(::Quad9,elements)
    n_elements = size(elements, 1)
    subdivided = Vector{Vector{Int}}(undef, 4*n_elements)
    for (i, element) in enumerate(eachrow(elements))
        n = element
        offset = 4i-3
        subdivided[offset+0] = [n[1], n[5], n[9], n[8]]
        subdivided[offset+1] = [n[5], n[2], n[6], n[9]]
        subdivided[offset+2] = [n[9], n[6], n[3], n[7]]
        subdivided[offset+3] = [n[8], n[9], n[7], n[4]]
    end
    return subdivided
end

function subdivide_node(::Quad16,elements)
    n_elements = size(elements, 1)
    subdivided = Vector{Vector{Int}}(undef, 9*n_elements)
    for (i, element) in enumerate(eachrow(elements))
        n = element
        offset = 9i-8
        subdivided[offset+0] = [n[1], n[5], n[13], n[12]]  # q1
        subdivided[offset+1] = [n[5], n[6], n[14], n[13]]  # q2
        subdivided[offset+2] = [n[6], n[2], n[7], n[14]]   # q3
        subdivided[offset+3] = [n[12], n[13], n[16], n[11]]# q4
        subdivided[offset+4] = [n[13], n[14], n[15], n[16]]# q5
        subdivided[offset+5] = [n[14], n[7], n[8], n[15]]  # q6
        subdivided[offset+6] = [n[11], n[16], n[10], n[4]] # q7
        subdivided[offset+7] = [n[16], n[15], n[9], n[10]] # q8
        subdivided[offset+8] = [n[15], n[8], n[3], n[9]]   # q9
    end
    return subdivided
end

function subdivide_node(::Quad25,elements)
    n_elements = size(elements, 1)
    subdivided = Vector{Vector{Int}}(undef, 16*n_elements)
    for (i, element) in enumerate(eachrow(elements))
        n = element
        offset = 16i-15
        subdivided[offset+0]  = [n[1], n[5], n[17], n[16]]   # q1
        subdivided[offset+1]  = [n[5], n[6], n[21], n[17]]   # q2
        subdivided[offset+2]  = [n[6], n[7], n[18], n[21]]   # q3
        subdivided[offset+3]  = [n[7], n[2], n[8], n[18]]    # q4
        subdivided[offset+4]  = [n[16], n[17], n[24], n[15]] # q5
        subdivided[offset+5]  = [n[17], n[21], n[25], n[24]] # q6
        subdivided[offset+6]  = [n[21], n[18], n[22], n[25]] # q7
        subdivided[offset+7]  = [n[18], n[8], n[9], n[22]]   # q8
        subdivided[offset+8]  = [n[15], n[24], n[20], n[14]] # q9
        subdivided[offset+9]  = [n[24], n[25], n[23], n[20]] # q10
        subdivided[offset+10] = [n[25], n[22], n[19], n[23]] # q11
        subdivided[offset+11] = [n[22], n[9], n[10], n[19]]  # q12
        subdivided[offset+12] = [n[14], n[20], n[13], n[4]]  # q13
        subdivided[offset+13] = [n[20], n[23], n[12], n[13]] # q14
        subdivided[offset+14] = [n[23], n[19], n[11], n[12]] # q15
        subdivided[offset+15] = [n[19], n[10], n[3], n[11]]  # q16
    end
    return subdivided
end

function subdivide_node(::Quad36, elements)
    n_elements = size(elements, 1)
    subdivided = Vector{Vector{Int}}(undef, 25*n_elements)
    
    for (i, element) in enumerate(eachrow(elements))
        n = element
        offset = 25i - 24
        
        # First row of subdivisions
        subdivided[offset+0]  = [n[1], n[5], n[21], n[20]]
        subdivided[offset+1]  = [n[5], n[6], n[25], n[21]]
        subdivided[offset+2]  = [n[6], n[7], n[26], n[25]]
        subdivided[offset+3]  = [n[7], n[8], n[22], n[26]]
        subdivided[offset+4]  = [n[8], n[2], n[9], n[22]]

        # Second row
        subdivided[offset+5]  = [n[20], n[21], n[32], n[19]]
        subdivided[offset+6]  = [n[21], n[25], n[33], n[32]]
        subdivided[offset+7]  = [n[25], n[26], n[34], n[33]]
        subdivided[offset+8]  = [n[26], n[22], n[27], n[34]]
        subdivided[offset+9]  = [n[22], n[9], n[10], n[27]]

        # Third row
        subdivided[offset+10] = [n[19], n[32], n[31], n[18]]
        subdivided[offset+11] = [n[32], n[33], n[36], n[31]]
        subdivided[offset+12] = [n[33], n[34], n[35], n[36]]
        subdivided[offset+13] = [n[34], n[27], n[28], n[35]]
        subdivided[offset+14] = [n[27], n[10], n[11], n[28]]

        # Fourth row
        subdivided[offset+15] = [n[18], n[31], n[24], n[17]]
        subdivided[offset+16] = [n[31], n[36], n[30], n[24]]
        subdivided[offset+17] = [n[36], n[35], n[29], n[30]]
        subdivided[offset+18] = [n[35], n[28], n[23], n[29]]
        subdivided[offset+19] = [n[28], n[11], n[12], n[23]]

        # Fifth row
        subdivided[offset+20] = [n[17], n[24], n[16], n[4]]
        subdivided[offset+21] = [n[24], n[30], n[15], n[16]]
        subdivided[offset+22] = [n[30], n[29], n[14], n[15]]
        subdivided[offset+23] = [n[29], n[23], n[13], n[14]]
        subdivided[offset+24] = [n[23], n[12], n[3], n[13]]
    end
    
    return subdivided
end

"""
    plot_champ_scalaire(nodes, elements, field, element_type)

Plot a scalar field on a mesh using Makie. Automatically handles element subdivision
and face type detection based on element type prefix (Tri* or Quad*).
"""
function plot_champ_scalaire(nodes, elements, field, element_type)
    subdivided_elements = subdivide_element(elements, element_type)
    
    # Create geometry objects
    points = [Point2f(row...) for row in eachrow(nodes)]
    
    # Determine face type based on element type prefix
    is_tri = startswith(string(element_type), "Tri")
    face_type = is_tri ? TriangleFace : QuadFace
    faces = [face_type(element...) for element in eachrow(subdivided_elements)]
    
    # Create and plot mesh
    mesh_obj = GeometryBasics.Mesh(points, faces)
    fig, ax, plt = mesh(mesh_obj, color = field, colormap = :jet)
    Colorbar(fig[1, 2], plt, label = "Temperature")
    
    return fig
end
