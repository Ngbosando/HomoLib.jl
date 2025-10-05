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
using Base.Threads
using StaticArrays
using ForwardDiff
using ProgressMeter
using Random
using Gmsh: gmsh
import FastGaussQuadrature.gausslegendre
using LoopVectorization
using GeometryBasics
using CairoMakie

abstract type AbstractMaterial end
abstract type AbstractElementType end

abstract type LinearElement        <: AbstractElementType end
abstract type TriangularElement    <: AbstractElementType end
abstract type QuadrilateralElement <: AbstractElementType end
abstract type TetrahedralElement   <: AbstractElementType end
abstract type HexahedralElement    <: AbstractElementType end
abstract type PrismaticElement     <: AbstractElementType end
abstract type PyramidElement       <: AbstractElementType end

# Triangles
struct Tri3  <: TriangularElement end
struct Tri6  <: TriangularElement end
struct Tri10 <: TriangularElement end
struct Tri15 <: TriangularElement end
struct Tri21 <: TriangularElement end

# Quadrilaterals
struct Quad4  <: QuadrilateralElement end
struct Quad9  <: QuadrilateralElement end
struct Quad16 <: QuadrilateralElement end
struct Quad25 <: QuadrilateralElement end
struct Quad36 <: QuadrilateralElement end

# Tetrahedra
struct Tet4  <: TetrahedralElement end
struct Tet10 <: TetrahedralElement end
struct Tet20 <: TetrahedralElement end
struct Tet35 <: TetrahedralElement end
struct Tet56 <: TetrahedralElement end

# Hexahedra
struct Hex8    <: HexahedralElement end
struct Hex27   <: HexahedralElement end
struct Hex64   <: HexahedralElement end
struct Hex125  <: HexahedralElement end
struct Hex216  <: HexahedralElement end

# Prisms
struct Pri6  <: PrismaticElement end
struct Pri18 <: PrismaticElement end
struct Pri36 <: PrismaticElement end
struct Pri56 <: PrismaticElement end
struct Pri78 <: PrismaticElement end

# Pyramids
struct Pyr5  <: PyramidElement end
struct Pyr14 <: PyramidElement end
struct Pyr29 <: PyramidElement end
struct Pyr50 <: PyramidElement end
struct Pyr77 <: PyramidElement end

# 1D (Gauss-Legendre)
struct Lin1 <: LinearElement end
struct Lin2 <: LinearElement end
struct Lin3 <: LinearElement end
struct Lin4 <: LinearElement end
struct Lin5 <: LinearElement end
struct Lin6 <: LinearElement end

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

include("Contact/Formulation for unilateral contac/spectral_strain_decomposition.jl")

include("exports.jl")

end