using Revise
using HomoLib
using Test

include("quadrature_test.jl")
include("precompute_data_test.jl")
include("global_assembly_test.jl")
include("analytical_bound_test.jl")
include("effective_properties_test.jl")
include("force_test.jl")
include("homogenization_test.jl")
include("material_model_test.jl")
include("patch_test.jl")
include("recover_field_values_test.jl")
include("shapes_function_test.jl")
include("test_global_assembly.jl")
include("SIMP_test.jl")
