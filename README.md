[![Build Status](https://github.com/leycrimson/HomoLib.jl/actions/workflows/CI.yml/badge.svg?branch=main)]
# HomoLib

**HomoLib.jl** is a Julia library for multi-physics computational homogenization and topology optimization based on the finite element method (FEM).

The project focuses on the development and validation of a modular FEM solver for heterogeneous materials and structural optimization problems.

---

## Current capabilities
- Finite element homogenization of heterogeneous materials
- Multi-physics support:
- elasticity
- thermal conduction
- poroelasticity
- piezoelectricity
- SIMP topology optimization with adaptive penalization
- Statistical convergence analysis for representative volume elements
- Validation against analytical benchmarks (Kirsch solution)
- Computation of effective material properties
- VTK export and visualization tools
- Automated test suite and CI integration
---
## Work in progress
- Fluid–structure interaction coupling
- Extended multi-physics coupling strategies
- Additional performance optimizations
---
## Installation
```julia
pkg> add HomoLib
```
---
## Usage
```julia
using HomoLib
# Homogenization analysis
C_eff = compute_effective_C(volume_fraction, n_inclusions, order, shape, element_type)
# Topology optimization
ρ_opt, history = topology_optimization_CB()
```
---
## Examples
- https://github.com/Ngbosando/HomoLib.jl/blob/main/Example/sensitivity_analysis.md
- https://github.com/Ngbosando/HomoLib.jl/blob/main/Example/topology_optimization_guide.md
- https://github.com/Ngbosando/HomoLib.jl/blob/main/Example/stress_homogenization_example.md
---
## Project goals
HomoLib is primarily developed as a research and training framework to deepen understanding of computational homogenization, composite materials modeling, and multi-physics finite element coupling.

The library serves as a testbed for validating numerical methods and exploring advanced FEM strategies in a modular Julia environment.
---
## License
MIT
