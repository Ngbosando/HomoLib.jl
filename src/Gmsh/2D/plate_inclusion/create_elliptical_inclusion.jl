
function create_elliptical_inclusion(
    x, y, a, b, θ
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
    surface = gmsh.model.geo.addPlaneSurface([loop])
   
    gmsh.model.geo.rotate([(2, surface)], x, y, 0.0, 0.0, 0.0, 1.0, θ)
    
    return loop, surface, C
end