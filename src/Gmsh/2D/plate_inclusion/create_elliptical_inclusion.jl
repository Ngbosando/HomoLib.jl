"""
    create_elliptical_inclusion(x, y, a, b, θ, voids)

    Create an elliptical inclusion in GMSH.

    # Arguments
    - "x": x-coordinate of center
    - "y": y-coordinate of center
    - "a": Semi-major axis length
    - "b": Semi-minor axis length
    - "θ": Rotation angle in radians
    - "voids": Whether to create void (true) or material (false)

    # Returns
    Tuple containing:
    - "loop": Curve loop tag
    - "surface": Surface tag (nothing if voids=true)
    - "boundary": List of elliptical arc tags

    # Notes
    - Creates ellipse using four elliptical arcs
    - Applies rotation after creation
"""

function create_elliptical_inclusion(
    x, y, a, b, θ, voids, to_rotate
)

    # Creation of points and arcs
   
    center = gmsh.model.geo.addPoint(x, y, 0.0)
    point_top = gmsh.model.geo.addPoint(x, y + b, 0.0)
    point_right = gmsh.model.geo.addPoint(x + a, y, 0.0)
    point_bottom = gmsh.model.geo.addPoint(x, y - b, 0.0)
    point_left = gmsh.model.geo.addPoint(x - a, y, 0.0)

    C1 = gmsh.model.geo.addEllipseArc(point_top, center, point_right, point_right)
    C2 = gmsh.model.geo.addEllipseArc(point_right, center, point_bottom, point_bottom)
    C3 = gmsh.model.geo.addEllipseArc(point_bottom, center, point_left, point_left)
    C4 = gmsh.model.geo.addEllipseArc(point_left, center, point_top, point_top)

    C = [C1, C2, C3, C4]
    
    loop = gmsh.model.geo.addCurveLoop([C1, C2, C3, C4])
    surface = voids ? nothing : gmsh.model.geo.addPlaneSurface([loop])
    
    # Rotate the square around its center
    if surface !== nothing
        rotate = [(2, surface)]                    # rotate the surface
    else
        rotate = [(1, l) for l in C]           # rotate the 4 lines
        # —or— rotate the points instead:
        # to_rotate = [(0, p) for p in pts]
    end
    if to_rotate
        gmsh.model.geo.rotate([rotate], x, y, 0.0, 0.0, 0.0, 1.0, θ)
    end
    
    return loop, surface, C
end