# """
#     resolve_material_type(material)

#     Determine the physics type of a material based on its properties.

#     # Arguments
#     - "material": A "Material" instance or vector of "Material" instances (for composites)

#     # Returns
#     Symbol representing the resolved material physics type:
#     - ":elastic": Pure elasticity
#     - ":piezoelectric": Coupled elastic-electric behavior
#     - ":thermal": Thermal conduction
#     - ":thermoelastic": Coupled thermal-elastic behavior
#     - ":poroelastic": Poroelasticity (Biot model)
#     - ":stokes": Stokes flow (fluid mechanics)

#     # Material Type Resolution Logic
#     1. "Pure Elastic": 
#     - Physics = [:elastic]
#     - No Biot coefficient (:α_p) present
#     2. "Piezoelectric":
#     - Physics = [:elastic, :electric]
#     3. "Thermal":
#     - Physics = [:thermal]
#     4. "Thermoelastic":
#     - Physics = [:elastic, :thermal]
#     5. "Poroelastic":
#     - Physics includes [:elastic, :pressure] OR 
#     - Material has Biot coefficient (:α_p)
#     6. "Stokes Flow":
#     - Physics = [:fluid, :pressure]
# """
# # function resolve_material_type(materials::Union{Material,Vector{Material}})
    
# #     # For composites, use the first material as the reference for determining the physics type.
# #     ref_material = materials isa Vector ? materials[1] : materials
    
# #     types = ref_material.type
# #     props = ref_material.properties
# #     num_types = length(types)

# #     # Use direct boolean checks on the `types` vector instead of creating a Set.
# #     # This is much more efficient when called inside a loop.
# #     has_elastic = :elastic in types
# #     has_electric = :electric in types
# #     has_thermal = :thermal in types
# #     has_pressure = :pressure in types
# #     has_fluid = :fluid in types

# #     if num_types == 1 && has_elastic && !haskey(props, :α_p) 
# #         return :elastic
# #     elseif num_types == 2 && has_elastic && has_electric
# #         return :piezoelectric
# #     elseif num_types == 1 && has_thermal
# #         return :thermal
# #     elseif num_types == 2 && has_elastic && has_thermal
# #         return :thermoelastic
# #     elseif (num_types == 2 && has_elastic && has_pressure) || haskey(props, :α_p)
# #         return :poroelastic
# #     elseif num_types == 2 && has_fluid && has_pressure
# #         return :stokes
# #     else
# #         error("Unknown material combination: $types")
# #     end
# # end

#  resolve_material_type(m::Material) = begin
#     # Decide your tag using m.B_types (fast, no strings)
#     bt = m.B_types
#     if (:strain in bt) && (:electric_field in bt)
#         :piezoelectric
#     elseif (:strain in bt)
#         :elastic
#     elseif (:temperature_gradient in bt)
#         :thermal
#     elseif (:velocity_gradient in bt) && (:pressure in bt)
#         :stokes
#     else
#         error("unsupported material B_types")
#     end
# end

# @noinline mismatch_error() = error("all element materials must have the same resolved type")

#  function resolve_material_type(ms::Vector{Material})
#     t = resolve_material_type(ms[1])
#     @inbounds for i in 2:length(ms)
#         resolve_material_type(ms[i]) === t || mismatch_error()
#     end
#     t
# end
