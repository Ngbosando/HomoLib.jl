function polygons_intersect(poly1::Vector{Tuple{Float64,Float64}}, poly2::Vector{Tuple{Float64,Float64}})

    axes1 = get_axes(poly1)
    axes2 = get_axes(poly2)
    axes = vcat(axes1, axes2)

    for axis in axes
        proj1 = [p[1]*axis[1] + p[2]*axis[2] for p in poly1]
        proj2 = [p[1]*axis[1] + p[2]*axis[2] for p in poly2]
        if maximum(proj1) < minimum(proj2) || maximum(proj2) < minimum(proj1)
            return false
        end
    end
    return true
end

function get_rotated_rectangle(x::Float64, y::Float64, width::Float64, height::Float64, θ::Float64)
    half_w = width / 2
    half_h = height / 2
    
    local_corners = [(-half_w, -half_h), ( half_w, -half_h),
                     ( half_w,  half_h), (-half_w,  half_h)]
    
    rotated_corners = [( x + cos(θ)*pt[1] - sin(θ)*pt[2],
                          y + sin(θ)*pt[1] + cos(θ)*pt[2] ) for pt in local_corners]
    return rotated_corners
end

function get_rotated_triangle(x::Float64, y::Float64, s::Float64, θ::Float64)
   
    h = s * sqrt(3)/2
    local_triangle = [(0.0, 2*h/3), (-s/2, -h/3), (s/2, -h/3)]
    rotated_triangle = [( x + cos(θ)*pt[1] - sin(θ)*pt[2],
                           y + sin(θ)*pt[1] + cos(θ)*pt[2] ) for pt in local_triangle]
    return rotated_triangle
end

function get_rotated_ellipse(x::Float64, y::Float64, a::Float64, b::Float64, θ::Float64; n::Int=20)
    points = Tuple{Float64,Float64}[]
    for i in 0:(n-1)
        t = 2*pi*i/n
        # Point on the ellipse in local coordinates
        local_point = (a*cos(t), b*sin(t))
        # Rotate and translate
        rotated_point = ( x + cos(θ)*local_point[1] - sin(θ)*local_point[2],
                          y + sin(θ)*local_point[1] + cos(θ)*local_point[2] )
        push!(points, rotated_point)
    end
    return points
end

"""
    create_inclusion(
        shape::Symbol, 
        x::Float64, 
        y::Float64, 
        size::Tuple{Float64, Float64}, 
        θ,
        existing_inclusions,
        voids
    )

    Create a single inclusion with collision checking.

    # Arguments
    - "shape": Shape of inclusion (":circle", ":square", ":ellipse", ":triangle")
    - "x": x-coordinate of center
    - "y": y-coordinate of center
    - "size": Tuple of size parameters (interpretation depends on shape)
    - "θ": Rotation angle in radians
    - "existing_inclusions": List of already placed inclusions
    - "voids": Whether to create void (true) or material inclusion (false)

    # Returns
    Tuple containing curve loop, surface, and boundary tags, or nothing if collision detected

    # Notes
    - Uses Separating Axis Theorem for precise collision detection
    - Supports rotated shapes
    - Delegates to shape-specific creation functions after collision check
"""
function create_inclusion(
    shape::Symbol, 
    x::Float64, 
    y::Float64, 
    size::Tuple{Float64, Float64}, 
    θ,
    existing_inclusions,
    voids,
    to_rotate)
    
    for inc in existing_inclusions
        if length(inc) == 5
            existing_x, existing_y, size1, size2, existing_θ = inc
        else
            existing_x, existing_y, size1, size2 = inc
            existing_θ = 0.0
        end

        if shape == :circle
            r_new = size[1]
            r_existing = size1
            dist = sqrt((x - existing_x)^2 + (y - existing_y)^2)
            if dist < (r_new + r_existing)
                return nothing
            end
        elseif shape == :square
            # Compute corners of the rotated rectangles.
            new_poly = get_rotated_rectangle(x, y, size[1], size[2], θ)
            existing_poly = get_rotated_rectangle(existing_x, existing_y, size1, size2, existing_θ)
            if polygons_intersect(new_poly, existing_poly)
                return nothing
            end
        elseif shape == :ellipse
            # Approximate ellipse boundary with polygon.
            new_poly = get_rotated_ellipse(x, y, size[1], size[2], θ)
            existing_poly = get_rotated_ellipse(existing_x, existing_y, size1, size2, existing_θ)
            if polygons_intersect(new_poly, existing_poly)
                return nothing
            end
        elseif shape == :triangle
            # Compute vertices of the rotated triangle.
            new_poly = get_rotated_triangle(x, y, size[1], θ)
            existing_poly = get_rotated_triangle(existing_x, existing_y, size1, existing_θ)
            if polygons_intersect(new_poly, existing_poly)
                return nothing
            end
        else
            error("Unknown shape: $shape")
        end
    end

    if shape == :circle
        return create_circular_inclusion(x, y, size[1], voids)
    elseif shape == :square
        return create_square_inclusion(x, y, size[1], size[2], θ, voids, to_rotate)
    elseif shape == :triangle
        return create_triangular_inclusion(x, y, size[1], θ, voids, to_rotate)
    elseif shape == :ellipse
        return create_elliptical_inclusion(x, y, size[1], size[2], θ, voids, to_rotate)
    else
        error("Unknown shape: $shape")
    end
end

function get_axes(poly)
    axes = Tuple{Float64,Float64}[]
    n = length(poly)
    for i in 1:n
        j = i % n + 1  
        edge = (poly[j][1] - poly[i][1], poly[j][2] - poly[i][2])
        # Compute normal (perpendicular vector)
        normal = (-edge[2], edge[1])
        mag = sqrt(normal[1]^2 + normal[2]^2)
        if mag != 0
            normal = (normal[1]/mag, normal[2]/mag)
        end
        push!(axes, normal)
    end
    return axes
end
