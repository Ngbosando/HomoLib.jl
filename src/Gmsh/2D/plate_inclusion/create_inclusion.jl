

# Utility function: Polygon intersection test using the Separating Axis Theorem (SAT)
function polygons_intersect(poly1::Vector{Tuple{Float64,Float64}}, poly2::Vector{Tuple{Float64,Float64}})
    # Compute the axes (normals to edges) to test for a given polygon.
    axes1 = get_axes(poly1)
    axes2 = get_axes(poly2)
    axes = vcat(axes1, axes2)

    # For each axis, project both polygons and check if there's a separating gap.
    for axis in axes
        proj1 = [p[1]*axis[1] + p[2]*axis[2] for p in poly1]
        proj2 = [p[1]*axis[1] + p[2]*axis[2] for p in poly2]
        if maximum(proj1) < minimum(proj2) || maximum(proj2) < minimum(proj1)
            # There is a separating axis → no overlap
            return false
        end
    end
    return true
end

# Returns the corners of a rotated rectangle (or square) centered at (x, y) with angle θ
function get_rotated_rectangle(x::Float64, y::Float64, width::Float64, height::Float64, θ::Float64)
    half_w = width / 2
    half_h = height / 2
    # Local corners (centered at origin)
    local_corners = [(-half_w, -half_h), ( half_w, -half_h),
                     ( half_w,  half_h), (-half_w,  half_h)]
    # Apply rotation and translation
    rotated_corners = [( x + cos(θ)*pt[1] - sin(θ)*pt[2],
                          y + sin(θ)*pt[1] + cos(θ)*pt[2] ) for pt in local_corners]
    return rotated_corners
end

# Returns the vertices of an equilateral triangle of side s, centered at (x, y), rotated by θ
function get_rotated_triangle(x::Float64, y::Float64, s::Float64, θ::Float64)
    # Height of the equilateral triangle
    h = s * sqrt(3)/2
    # To center it at the centroid:
    # v1 = (0, 2h/3), v2 = (-s/2, -h/3), v3 = (s/2, -h/3)
    local_triangle = [(0.0, 2*h/3), (-s/2, -h/3), (s/2, -h/3)]
    rotated_triangle = [( x + cos(θ)*pt[1] - sin(θ)*pt[2],
                           y + sin(θ)*pt[1] + cos(θ)*pt[2] ) for pt in local_triangle]
    return rotated_triangle
end

# Approximates the boundary of an ellipse with n rotated points
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

# Main function to create an inclusion with overlap detection, taking into account the angle
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
    to_rotate
)
    for inc in existing_inclusions
        # Extract parameters of the existing inclusion.
        # If angle is not present, assume 0.0.
        if length(inc) == 5
            existing_x, existing_y, size1, size2, existing_θ = inc
        else
            existing_x, existing_y, size1, size2 = inc
            existing_θ = 0.0
        end

        if shape == :circle
            # Rotation doesn't matter for circles.
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

    # If no overlap is detected, create the inclusion based on its shape.
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

# Returns the normals (axes) to the edges of a polygon (for SAT test)
function get_axes(poly)
    axes = Tuple{Float64,Float64}[]
    n = length(poly)
    for i in 1:n
        j = i % n + 1  # next index (wrap-around)
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
