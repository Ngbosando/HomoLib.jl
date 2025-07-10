module HomoLib

using SparseArrays, LinearAlgebra, Statistics
using Reexport, ForwardDiff, Parameters
using Base.Threads, StaticArrays, ProgressMeter
using NLsolve, IterativeSolvers, LoopVectorization
using Gmsh: gmsh
using Random, Tensors, Combinatorics
using GeometryBasics, CairoMakie, DelaunayTriangulation
import FastGaussQuadrature.gausslegendre, MAT
using StatsBase,BenchmarkTools



struct Material 
    type::Vector{Symbol}
    dim::Int
    symmetry::Symbol
    properties::Dict{Symbol, Union{Number, Vector, Matrix}}
    tensors::Vector{AbstractMatrix{<:Real}}
    B_types::Vector{Symbol}
    mass_properties::Union{Nothing, Dict{Symbol, Any}}
end

abstract type TriangularElement   end
struct Tri3  <: TriangularElement end    # 3-nodes Linear 
struct Tri6  <: TriangularElement end    # 6-nodes Quadratic 
struct Tri10 <: TriangularElement end    # 10-nodes Lagrangian
struct Tri15 <: TriangularElement end    # 15-nodes quartic 
struct Tri21 <: TriangularElement end    # 21-nodes quintic 

abstract type HexahedralElement     end
struct Hex8  <: HexahedralElement end   # 8-node linear
struct Hex27 <: HexahedralElement end  # 27-node quadratic
struct Hex64 <: HexahedralElement end  # 64-node cubic
struct Hex125 <: HexahedralElement end # 125-node quartic
struct Hex216 <: HexahedralElement end # 216-node quintic

abstract type PrismaticElement      end
struct Pri6  <: PrismaticElement end  # 6-node linear prism
struct Pri18 <: PrismaticElement end # 18-node quadratic prism
struct Pri36 <: PrismaticElement end # 36-node cubic prism
struct Pri56 <: PrismaticElement end # 56-node quartic prism
struct Pri78 <: PrismaticElement end # 78-node quintic prism

abstract type PyramidElement        end
struct Pyr5  <: PyramidElement end   # 5-nodes Linear pyramid 
struct Pyr14 <: PyramidElement end  # 14-nodes Quadratic pyramid 
struct Pyr29 <: PyramidElement end  # 29-nodes Lagrangian pyramid
struct Pyr50 <: PyramidElement end  # 50-nodes cubic pyramid
struct Pyr77 <: PyramidElement end  # 77-nodes quartic pyramid

abstract type QuadrilateralElement  end
struct Quad4 <: QuadrilateralElement end  # 4-node bilinear
struct Quad9 <: QuadrilateralElement end  # 8-node serendipity
struct Quad16 <: QuadrilateralElement end  # 16-node Lagrangian
struct Quad25 <: QuadrilateralElement end  # 25-node Lagrangian
struct Quad36 <: QuadrilateralElement end  # 36-node Lagrangian

abstract type TetrahedralElement    end
struct Tet4  <: TetrahedralElement end  # 4-node linear 
struct Tet10 <: TetrahedralElement end # 10-node quadratic 
struct Tet20 <: TetrahedralElement end # 20-node cubic
struct Tet35 <: TetrahedralElement end # 35-node quartic
struct Tet56 <: TetrahedralElement end # 56-node quintic

abstract type LinearElement         end
struct Lin1 <: LinearElement end  # 1 points d'intégration
struct Lin2 <: LinearElement end  # 2 points d'intégration
struct Lin3 <: LinearElement end  # 3 points d'intégration
struct Lin4 <: LinearElement end  # 4 points d'intégration
struct Lin5 <: LinearElement end  # 5 points d'intégration
struct Lin6 <: LinearElement end  # 6 points d'intégration
# =============================================
# Material System
# =============================================




include("Analytical model/analytical_bound.jl")
include("Shapes functions/shapes_interpolation.jl")
include("Quadrature/quadrature.jl")
include("Assembler/assemble.jl")
include("Solve/solver.jl")
include("Compute properties/effective_properties.jl")
include("Gmsh/2D/plate_inclusion/generate_transfinite_plate_with_inclusions.jl")
include("Gmsh/2D/Simple plate/plaque.jl")
include("Gmsh/2D/Simple plate/quad_patch.jl")
include("Gmsh/2D/Simple plate/get_boundary.jl")
include("Plots/2D/plot_champ_scalaire.jl")
include("Plots/2D/plot_structure_with_centroids.jl")
include("Plots/3D/build_periodic_node_pairs.jl")
include("Plots/3D/visualize_periodic_pairs.jl")


include("exports.jl")



end
