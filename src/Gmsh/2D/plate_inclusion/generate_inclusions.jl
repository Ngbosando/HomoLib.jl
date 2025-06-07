function generate_inclusions(volume_fraction, 
    min_size, 
    max_size, 
    plate_width, 
    plate_height, 
    shape, 
    N_inclu
)

inclusion_loops = Int[]
inclusion_surfaces = Int[]
inclusion_borders = []
existing_inclusions = []

volume_total_inclusions = volume_fraction
attempt_count = 0
max_attempts = 1e4

# Détermination des tailles de base en fonction de la forme
if shape == :circle
    r = plate_width * sqrt(volume_fraction / (N_inclu * π))   
    b = r    
elseif shape == :square
    k = 0.5
    r = plate_width * sqrt(volume_fraction / (N_inclu * k))   
    b = k * r   
elseif shape == :ellipse
    k = 0.5
    r = plate_width * sqrt(volume_fraction / (N_inclu * k * π))   
    b = k * r * π 
elseif shape == :triangle
    k = 0.5
    r = plate_width * sqrt(volume_fraction / (N_inclu * k * 2))   
    b = k * r * 2  
else 
    error("Unknown shape: $shape")
    gmsh.finalize()
end

remaining_volume = volume_total_inclusions
attempt_inclu = 0
num_inclusions = 0

function triangle_margins(r, θ)
    v1 = (-r/2, -r/3)
    v2 = ( r/2, -r/3)
    v3 = (0, 2r/3)
    # Fonction de rotation
    rotate(v, θ) = (cos(θ)*v[1] - sin(θ)*v[2], sin(θ)*v[1] + cos(θ)*v[2])
    rv1 = rotate(v1, θ)
    rv2 = rotate(v2, θ)
    rv3 = rotate(v3, θ)
    margin_x = maximum(abs.([rv1[1], rv2[1], rv3[1]]))
    margin_y = maximum(abs.([rv1[2], rv2[2], rv3[2]]))
return margin_x, margin_y
end

while remaining_volume >= 1e-3 && attempt_count < 1e6 

    attempt_count += 1
    θ = 2*pi*rand() - pi

    if shape == :circle

        margin_x = r
        margin_y = r
    elseif shape == :square

        half_side = r/2
        margin_x = half_side*(abs(cos(θ)) + abs(sin(θ)))
        margin_y = margin_x

    elseif shape == :ellipse

        margin_x = sqrt((r*cos(θ))^2 + (b*sin(θ))^2)
        margin_y = sqrt((r*sin(θ))^2 + (b*cos(θ))^2)

    elseif shape == :triangle
        margin_x, margin_y = triangle_margins(r, θ)
    else
    error("Unknown shape: $shape")
    end

    # Si les marges dépassent la moitié des dimensions de la plaque,
    # on passe cette tentative (cela peut arriver si la forme est trop grande).
    if margin_x > plate_width/2 || margin_y > plate_height/2
    continue
    end


    # xmin = margin_x
    # xmax = plate_width - margin_x
    # ymin = margin_y
    # ymax = plate_height - margin_y

    # x = xmin + (xmax - xmin)*rand()
    # y = ymin + (ymax - ymin)*rand()
    x = 1/2
    y = 1/2

    result = create_inclusion(shape,
            x, 
            y, 
            (r, b),
            θ,
            existing_inclusions
            )

    if result !== nothing
        loop, surface, border = result
        
        # On stocke aussi θ pour pouvoir vérifier la position et la rotation ultérieurement
        push!(existing_inclusions, (x, y, r, b, θ))
        push!(inclusion_loops, loop)
        push!(inclusion_surfaces, surface)
        push!(inclusion_borders, border)
        num_inclusions += 1 

        if shape == :circle
            remaining_volume -= (π * r^2) / plate_width^2
        elseif shape == :square
            remaining_volume -= (r * b) / plate_width^2
        elseif shape == :triangle
            remaining_volume -= ((r * b) / 2) / plate_width^2
        elseif shape == :ellipse
            remaining_volume -= (π * r * b) / plate_width^2
        end
        attempt_inclu = 0 
    else
        attempt_inclu += 1
    end

    if attempt_inclu == max_attempts
        N_inclu += 1
        if shape == :circle
        r = plate_width * sqrt(remaining_volume / (N_inclu * π))       
        elseif shape == :square
        k = 0.5
        r = plate_width * sqrt(remaining_volume / (N_inclu * k))   
        b = k * r
        elseif shape == :ellipse
        k = 0.5
        r = plate_width * sqrt(remaining_volume / (N_inclu * k * π))   
        b = k * r * π
        elseif shape == :triangle
        k = 0.5
        r = plate_width * sqrt(remaining_volume / (N_inclu * k * 2))   
        b = k * r * 2
    else 
        println("Forme inconnue : $shape")
    end
    attempt_inclu = 0  
    end
end

return inclusion_loops, inclusion_surfaces, inclusion_borders, num_inclusions
end
