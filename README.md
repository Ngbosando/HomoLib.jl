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

- `examples/composite_homogenization.jl` - Statistical analysis of composite materials with convergence studies
- `examples/topology_optimization_cb.jl` - SIMP method for cantilever beam optimization (based on Sigmund's 99-line code)

Run examples to see full capabilities and generated plots.

## Features

- Finite element homogenization of heterogeneous materials
- SIMP topology optimization with adaptive penalization
- Statistical convergence analysis
- Visualization and VTK export

