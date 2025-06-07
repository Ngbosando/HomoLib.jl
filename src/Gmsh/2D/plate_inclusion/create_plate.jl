
function create_plate(
    plate_width, 
    plate_height
)
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