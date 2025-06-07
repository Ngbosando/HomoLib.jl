
"""
    resolve_material_type(physics::Vector{Symbol}) -> Symbol

Resolves a composite material `physics` vector (e.g. `[:elastic, :thermal]`) 
into a legacy type symbol used for dispatch.
"""
function resolve_material_type(physics::Vector{Symbol})
    s = Set(physics)

    if s == Set([:elastic])
        return :elastic
    elseif s == Set([:elastic, :electric])
        return :piezoelectric
    elseif s == Set([:thermal])
        return :thermal
    elseif s == Set([:elastic, :thermal])
        return :thermoelastic
    elseif s == Set([:elastic, :thermal, :electric])
        return :thermo_piezoelectric
    elseif s == Set([:strain])
        return :viscoelastic
    elseif s == Set([:elastic, :pressure])
        return :poroelastic
    else
        error("Unknown material combination: $physics")
    end
end

