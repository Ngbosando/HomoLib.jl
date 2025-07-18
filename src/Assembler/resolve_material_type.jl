"""
    resolve_material_type(material)

    Determine the physics type of a material based on its properties.

    # Arguments
    - "material": A "Material" instance or vector of "Material" instances (for composites)

    # Returns
    Symbol representing the resolved material physics type:
    - ":elastic": Pure elasticity
    - ":piezoelectric": Coupled elastic-electric behavior
    - ":thermal": Thermal conduction
    - ":thermoelastic": Coupled thermal-elastic behavior
    - ":poroelastic": Poroelasticity (Biot model)
    - ":stokes": Stokes flow (fluid mechanics)

    # Material Type Resolution Logic
    1. "Pure Elastic": 
    - Physics = [:elastic]
    - No Biot coefficient (:α_p) present
    2. "Piezoelectric":
    - Physics = [:elastic, :electric]
    3. "Thermal":
    - Physics = [:thermal]
    4. "Thermoelastic":
    - Physics = [:elastic, :thermal]
    5. "Poroelastic":
    - Physics includes [:elastic, :pressure] OR 
    - Material has Biot coefficient (:α_p)
    6. "Stokes Flow":
    - Physics = [:fluid, :pressure]
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

