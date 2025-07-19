"""
    create_square_inclusion(x, y, size1, size2, θ, voids)

    Create a rotated square inclusion in GMSH.

    # Arguments
    - "x": x-coordinate of center
    - "y": y-coordinate of center
    - "size1": Length of side 1
    - "size2": Length of side 2
    - "θ": Rotation angle in radians
    - "voids": Whether to create void (true) or material (false)

    # Returns
    Tuple containing:
    - "loop": Curve loop tag
    - "surface": Surface tag (nothing if voids=true)
    - "boundary": List of line tags

    # Notes
    - Creates square centered at (x,y) and applies rotation
"""

function create_square_inclusion(x, y, size1, size2, θ,voids, to_rotate)
   
 
    a = size1 / 2
    b = size2 / 6
    
    # Creation of points and lines
    plate_points = [
        gmsh.model.geo.addPoint(x - a, y - b, 0.0),
        gmsh.model.geo.addPoint(x + a, y - b, 0.0),
        gmsh.model.geo.addPoint(x + a, y + b, 0.0),
        gmsh.model.geo.addPoint(x - a, y + b, 0.0)
    ]

    plate_lines = [
        gmsh.model.geo.addLine(plate_points[1], plate_points[2]),
        gmsh.model.geo.addLine(plate_points[2], plate_points[3]),
        gmsh.model.geo.addLine(plate_points[3], plate_points[4]),
        gmsh.model.geo.addLine(plate_points[4], plate_points[1])
    ]


    plate_loop = gmsh.model.geo.addCurveLoop(plate_lines)
    surface = voids ? nothing : gmsh.model.geo.addPlaneSurface([plate_loop])
    
    # Rotate the square around its center
    if surface !== nothing
        rotate = [(2, surface)]                    # rotate the surface
    else
        rotate = [(1, l) for l in plate_lines]           # rotate the 4 lines
        # —or— rotate the points instead:
        # to_rotate = [(0, p) for p in pts]
    end
    if to_rotate
        gmsh.model.geo.rotate(rotate, x, y, 0.0, 0.0, 0.0, 1.0, θ)
    end
    return plate_loop, surface, plate_lines
end