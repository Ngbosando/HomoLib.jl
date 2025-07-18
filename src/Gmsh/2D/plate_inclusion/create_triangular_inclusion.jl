"""
    create_triangular_inclusion(x, y, side_length, θ, voids)

    Create an equilateral triangular inclusion in GMSH.

    # Arguments
    - "x": x-coordinate of center
    - "y": y-coordinate of center
    - "side_length": Length of triangle sides
    - "θ": Rotation angle in radians
    - "voids": Whether to create void (true) or material (false)

    # Returns
    Tuple containing:
    - "loop": Curve loop tag
    - "surface": Surface tag (nothing if voids=true)
    - "boundary": List of line tags

    # Notes
    - Creates equilateral triangle centered at (x,y) and applies rotation
"""

function create_triangular_inclusion(x, y, side_length, θ, voids)
    height = sqrt(3) / 2 * side_length  # Calculate the height of the equilateral triangle
   
    # Creation of points and lines
    
    point1 = gmsh.model.geo.addPoint(x, y + 2 * height / 3, 0.0)
    point1 = gmsh.model.geo.addPoint(x, y + 2 * height / 3, 0.0)
    point2 = gmsh.model.geo.addPoint(x - side_length / 2, y - height / 3, 0.0)
    point3 = gmsh.model.geo.addPoint(x + side_length / 2, y - height / 3, 0.0)

    L1 = gmsh.model.geo.addLine(point1, point2)
    L2 = gmsh.model.geo.addLine(point2, point3)
    L3 = gmsh.model.geo.addLine(point3, point1)
    L = [L1, L2, L3]

    loop = gmsh.model.geo.addCurveLoop([L1, L2, L3])
    surface = voids ? nothing : gmsh.model.geo.addPlaneSurface([loop])
   
    gmsh.model.geo.rotate([(2, surface)], x, y, 0.0, 0.0, 0.0, 1.0, θ)
    
    return loop, surface, L
end