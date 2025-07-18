"""
    create_circular_inclusion(x, y, radius, voids)

    Create a circular inclusion in GMSH.

    # Arguments
    - "x": x-coordinate of center
    - "y": y-coordinate of center
    - "radius": Radius of circle
    - "voids": Whether to create void (true) or material (false)

    # Returns
    Tuple containing:
    - "loop": Curve loop tag
    - "surface": Surface tag (nothing if voids=true)
    - "boundary": List of circular arc tags

    # Notes
    - Creates circle using four circular arcs
    - For voids, no surface is created
"""

function create_circular_inclusion(x, y, radius, voids)

    # Creation of points and arcs

    center = gmsh.model.geo.addPoint(x, y, 0.0)
    point_top = gmsh.model.geo.addPoint(x, y + radius, 0.0)
    point_right = gmsh.model.geo.addPoint(x + radius, y, 0.0)
    point_bottom = gmsh.model.geo.addPoint(x, y - radius, 0.0)
    point_left = gmsh.model.geo.addPoint(x - radius, y, 0.0)

    C1 = gmsh.model.geo.addCircleArc(point_top, center, point_right)
    C2 = gmsh.model.geo.addCircleArc(point_right, center, point_bottom)
    C3 = gmsh.model.geo.addCircleArc(point_bottom, center, point_left)
    C4 = gmsh.model.geo.addCircleArc(point_left, center, point_top)
    C = [C1, C2, C3, C4]

    loop = gmsh.model.geo.addCurveLoop([C1, C2, C3, C4])
    surface = voids ? nothing : gmsh.model.geo.addPlaneSurface([loop])


    return loop, surface, C
end