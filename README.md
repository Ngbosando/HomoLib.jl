[![Build Status](https://github.com/leycrimson/HomoLib.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/leycrimson/HomoLib.jl/actions/workflows/CI.yml?query=branch%3Amain)
# HomoLib

Computational homogenization and topology optimization in Julia.

## Installation

```julia
pkg> add HomoLib
```

## Usage

```julia
using HomoLib

# Homogenization analysis
C_eff = compute_effective_C(volume_fraction, n_inclusions, order, shape, element_type)

# Topology optimization  
ρ_opt, history = topology_optimization_CB()
```

## Examples
 * [Statistical analysis of composite materials with convergence studies.](https://github.com/Ngbosando/HomoLib.jl/blob/main/Example/sensitivity_analysis.md)
- `examples/topology_optimization_cb.jl` - SIMP method for cantilever beam optimization (based on Sigmund's 99-line code)
- `examples/stress_homogenization_example.md` – Validation of stress analysis (Kirsch solution for a perforated plate) and poroelastic homogenization with comparison to analytical estimates

Run examples to see full capabilities and generated plots.

## Features

- Finite element homogenization of heterogeneous materials
- SIMP topology optimization with adaptive penalization
- Statistical convergence analysis
- Stress analysis validation against analytical references (Kirsch)
- Effective poroelastic properties with numerical vs analytical comparison
- Visualization and VTK export

