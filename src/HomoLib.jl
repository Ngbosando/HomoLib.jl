"""
    HomoLib
    
    A Julia package for computational homogenization and finite element analysis.

    # Modules
    - "Analytical Models": Contains analytical bounds for material properties
    - "Shape Functions": Provides interpolation functions for various element types
    - "Quadrature": Numerical integration methods
    - "Assembler": Finite element matrix assembly routines
    - "Solver": Numerical solvers for FEM systems
    - "Effective Properties": Computes effective material properties
    - "Gmsh Integration": Utilities for working with Gmsh meshes (2D and 3D)
    - "Visualization": Plotting utilities for 2D and 3D results

    # Element Types
    The package supports a wide range of finite elements:
    - Triangular elements (Tri3, Tri6, Tri10, Tri15, Tri21)
    - Hexahedral elements (Hex8, Hex27, Hex64, Hex125, Hex216)
    - Prismatic elements (Pri6, Pri18, Pri36, Pri56, Pri78)
    - Pyramid elements (Pyr5, Pyr14, Pyr29, Pyr50, Pyr77)
    - Quadrilateral elements (Quad4, Quad9, Quad16, Quad25, Quad36)
    - Tetrahedral elements (Tet4, Tet10, Tet20, Tet35, Tet56)
    - Linear elements (Lin1-Lin6)

    # Material System
    The `Material` struct represents material properties with:
    - Type (e.g., :elastic, :poroelastic)
    - Dimension
    - Symmetry type
    - Mechanical properties
    - Tensor representations
    - Boundary condition types
    - Mass properties

    # Dependencies
    The package relies on several Julia packages including:
    - LinearAlgebra, SparseArrays, Statistics (standard library)
    - ForwardDiff, Parameters, StaticArrays
    - NLsolve, IterativeSolvers
    - GeometryBasics, Tensors
    - CairoMakie for visualization
    - Gmsh for mesh generation

    # Examples
    ```julia
    using HomoLib

    # Create a material
    mat = Material([:elastic], 3, :isotropic, Dict(:E => 210e9, :ν => 0.3), ...)

    # Generate a mesh
    nodes, elements = generate_plate_with_inclusions(...)

    # Run analysis
    results = solve_fem_problem(nodes, elements, mat)
"""
module HomoLib

using SparseArrays 
using LinearAlgebra
using Statistics
using Reexport 
using ForwardDiff
using Parameters
using Base.Threads 
using StaticArrays
using ProgressMeter
using NLsolve 
using IterativeSolvers
using LoopVectorization
using Gmsh: gmsh
using Random 
using Tensors
using Combinatorics
using GeometryBasics
using Dates
using CairoMakie
import FastGaussQuadrature.gausslegendre
using MAT
using StatsBase 
using DelaunayTriangulation
using BenchmarkTools
using WriteVTK


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
include("Gmsh/3D/simple_struct/build_reference_tetrahedron.jl")
include("Plots/2D/plot_champ_scalaire.jl")
include("Plots/2D/plot_structure_with_centroids.jl")
include("Plots/3D/build_periodic_node_pairs.jl")
include("Plots/3D/visualize_periodic_pairs.jl")
include("exports.jl")



end
