# Fonction utilitaire : Test d'intersection de deux polygones avec la méthode des axes séparants (SAT)
function polygons_intersect(poly1::Vector{Tuple{Float64,Float64}}, poly2::Vector{Tuple{Float64,Float64}})
    # Calcule les axes (normales aux arêtes) à tester pour un polygone donné.

    axes1 = get_axes(poly1)
    axes2 = get_axes(poly2)
    axes = vcat(axes1, axes2)
    # Pour chacun des axes, on projette les polygones et on vérifie s'il y a séparation.
    for axis in axes
        proj1 = [p[1]*axis[1] + p[2]*axis[2] for p in poly1]
        proj2 = [p[1]*axis[1] + p[2]*axis[2] for p in poly2]
        if maximum(proj1) < minimum(proj2) || maximum(proj2) < minimum(proj1)
            # Il existe un axe séparant les deux polygones → pas de recouvrement
            return false
        end
    end
    return true
end

# Renvoie les coins d'un rectangle (ou carré) tourné d'un angle θ autour de son centre (x,y)
function get_rotated_rectangle(x::Float64, y::Float64, width::Float64, height::Float64, θ::Float64)
    half_w = width / 2
    half_h = height / 2
    # Coins locaux (centre en (0,0))
    local_corners = [(-half_w, -half_h), ( half_w, -half_h),
                     ( half_w,  half_h), (-half_w,  half_h)]
    # Rotation et translation
    rotated_corners = [( x + cos(θ)*pt[1] - sin(θ)*pt[2],
                          y + sin(θ)*pt[1] + cos(θ)*pt[2] ) for pt in local_corners]
    return rotated_corners
end

# Renvoie les sommets d'un triangle équilatéral de côté s, centré en (x,y), tourné de θ.
function get_rotated_triangle(x::Float64, y::Float64, s::Float64, θ::Float64)
    # Hauteur du triangle équilatéral
    h = s * sqrt(3)/2
    # Pour centrer le triangle au niveau de son centroïde, on peut utiliser :
    # v1 = (0, 2h/3), v2 = (-s/2, -h/3), v3 = (s/2, -h/3)
    local_triangle = [(0.0, 2*h/3), (-s/2, -h/3), (s/2, -h/3)]
    rotated_triangle = [( x + cos(θ)*pt[1] - sin(θ)*pt[2],
                           y + sin(θ)*pt[1] + cos(θ)*pt[2] ) for pt in local_triangle]
    return rotated_triangle
end

# Approxime le contour d'une ellipse par n points
function get_rotated_ellipse(x::Float64, y::Float64, a::Float64, b::Float64, θ::Float64; n::Int=20)
    points = Tuple{Float64,Float64}[]
    for i in 0:(n-1)
        t = 2*pi*i/n
        # Point sur l'ellipse dans le repère local
        local_point = (a*cos(t), b*sin(t))
        # Rotation et translation
        rotated_point = ( x + cos(θ)*local_point[1] - sin(θ)*local_point[2],
                          y + sin(θ)*local_point[1] + cos(θ)*local_point[2] )
        push!(points, rotated_point)
    end
    return points
end

# Fonction principale de création d'inclusion avec détection de chevauchement prenant en compte l'angle.
function create_inclusion(shape::Symbol, 
    x::Float64, 
    y::Float64, 
    size::Tuple{Float64, Float64}, 
    θ,
    existing_inclusions
)
    for inc in existing_inclusions
        # Extraction des paramètres de l'inclusion existante.
        # Si l'angle n'est pas présent, on suppose qu'il vaut 0.0.
        if length(inc) == 5
            existing_x, existing_y, size1, size2, existing_θ = inc
        else
            existing_x, existing_y, size1, size2 = inc
            existing_θ = 0.0
        end

        if shape == :circle
            # Pour les cercles, la rotation n'a pas d'importance.
            r_new = size[1]
            r_existing = size1
            dist = sqrt((x - existing_x)^2 + (y - existing_y)^2)
            if dist < (r_new + r_existing)
                return nothing
            end
        elseif shape == :square
            # On calcule les coins des rectangles tournés.
            new_poly = get_rotated_rectangle(x, y, size[1], size[2], θ)
            existing_poly = get_rotated_rectangle(existing_x, existing_y, size1, size2, existing_θ)
            if polygons_intersect(new_poly, existing_poly)
                return nothing
            end
        elseif shape == :ellipse
            # Approximation par polygone du contour de l'ellipse.
            new_poly = get_rotated_ellipse(x, y, size[1], size[2], θ)
            existing_poly = get_rotated_ellipse(existing_x, existing_y, size1, size2, existing_θ)
            if polygons_intersect(new_poly, existing_poly)
                return nothing
            end
        elseif shape == :triangle
            # Calcul des sommets du triangle tourné.
            new_poly = get_rotated_triangle(x, y, size[1], θ)
            existing_poly = get_rotated_triangle(existing_x, existing_y, size1, existing_θ)
            if polygons_intersect(new_poly, existing_poly)
                return nothing
            end
        else
            error("Unknown shape: $shape")
        end
    end
    
    # Si aucun chevauchement n'est détecté, on crée l'inclusion selon la forme.
    if shape == :circle
        return create_circular_inclusion(x, y, size[1])
    elseif shape == :square
        return create_square_inclusion(x, y, size[1], size[2], θ)
    elseif shape == :triangle
        return create_triangular_inclusion(x, y, size[1], θ)
    elseif shape == :ellipse
        return create_elliptical_inclusion(x, y, size[1], size[2], θ)
    else
        error("Unknown shape: $shape")
    end
end
function get_axes(poly)
    axes = Tuple{Float64,Float64}[]
    n = length(poly)
    for i in 1:n
        j = i % n + 1  # indice suivant (boucle circulaire)
        edge = (poly[j][1] - poly[i][1], poly[j][2] - poly[i][2])
        # La normale (on peut choisir l'orientation opposée)
        normal = (-edge[2], edge[1])
        mag = sqrt(normal[1]^2 + normal[2]^2)
        if mag != 0
            normal = (normal[1]/mag, normal[2]/mag)
        end
        push!(axes, normal)
    end
    return axes
end