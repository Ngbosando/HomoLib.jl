
function create_square_inclusion(x, y, size1, size2, θ)
   
 
    a = size1 / 2
    b = size2 / 2
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
    surface = gmsh.model.geo.addPlaneSurface([plate_loop])
 
    gmsh.model.geo.rotate([(2, surface)], x, y, 0.0, 0.0, 0.0, 1.0, θ)

    return plate_loop, surface, plate_lines
end