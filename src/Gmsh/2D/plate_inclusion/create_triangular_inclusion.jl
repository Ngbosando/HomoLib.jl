

function create_triangular_inclusion(x, y, side_length, θ)
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
    surface = gmsh.model.geo.addPlaneSurface([loop])
   
    gmsh.model.geo.rotate([(2, surface)], x, y, 0.0, 0.0, 0.0, 1.0, θ)
    
    return loop, surface, L
end