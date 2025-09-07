"""
    create_plate(plate_width, plate_height)

    Create a rectangular plate geometry in GMSH.

    # Arguments
    - "plate_width": Width of the plate (x-direction)
    - "plate_height": Height of the plate (y-direction)

    # Returns
    Tuple containing:
    - "plate_loop": Curve loop tag of the plate boundary
    - "plate_lines": List of line tags for plate edges

    # Notes
    - Creates a rectangle with corners at (0,0) to (width,height)
    - Returns geometry tags for further operations
"""

function create_plate(
    plate_width, 
    plate_height)

    plate_points = [
        gmsh.model.geo.addPoint(0.0, 0.0, 0.0),
        gmsh.model.geo.addPoint(plate_width, 0.0, 0.0),
        gmsh.model.geo.addPoint(plate_width, plate_height, 0.0),
        gmsh.model.geo.addPoint(0.0, plate_height, 0.0)
    ]

    plate_lines = [
        gmsh.model.geo.addLine(plate_points[1], plate_points[2]),
        gmsh.model.geo.addLine(plate_points[2], plate_points[3]),
        gmsh.model.geo.addLine(plate_points[3], plate_points[4]),
        gmsh.model.geo.addLine(plate_points[4], plate_points[1])
    ]

    plate_loop = gmsh.model.geo.addCurveLoop(plate_lines)

    return plate_loop, plate_lines
end