
"""
    resolve_material_type(physics::Vector{Symbol}) -> Symbol

Resolves a composite material `physics` vector (e.g. `[:elastic, :thermal]`) 
into a legacy type symbol used for dispatch.
"""
function resolve_material_type(materials::Union{Material,Vector{Material}})
 
    physics = materials.type
    s = Set(physics)

    if s == Set([:elastic]) && !haskey(materials.properties, :α_p) 
        return :elastic
    elseif s == Set([:elastic, :electric])
        return :piezoelectric
    elseif s == Set([:thermal])
        return :thermal
    elseif s == Set([:elastic, :thermal])
        return :thermoelastic
    elseif s == Set([:elastic, :pressure]) || haskey(materials.properties, :α_p)
        return :poroelastic
    elseif s == Set([:fluid, :pressure])
        return :stokes
    else
        error("Unknown material combination: $physics")
    end
end

