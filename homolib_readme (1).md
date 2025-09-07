# HomoLib.jl

A Julia package for analyzing composite materials and optimizing structures using finite elements.

## What it does

HomoLib.jl helps you:
- Find effective properties of composite materials
- Optimize structures to be stronger and lighter
- Analyze materials with holes, inclusions, or different phases
- Create visualizations and export results

## Installation

```julia
using Pkg
Pkg.add("HomoLib")
```

## Quick Example

```julia
using HomoLib

# Analyze a composite with 40% inclusions
C_effective = compute_effective_C(0.4, 10, 2, :square, :Tri6, 10, 20)
```

## Examples Included

### 1. Composite Material Analysis
**File**: `examples/composite_homogenization.jl`

Shows how to analyze materials with inclusions (like concrete with aggregates or metals with voids):
- Tests how many samples you need for good results
- Compares different material models
- Creates plots showing material properties vs. inclusion amount

**What you get**:
- `convergence_study.png` - Shows when you have enough samples
- `sensitivity_analysis.png` - How properties change with inclusion content
- `constitutive_law_*.png` - Stress-strain behavior plots

### 2. Structure Optimization
**File**: `examples/topology_optimization_cb.jl`

Optimizes a beam structure to be as stiff as possible with limited material:
- Based on famous "99-line" Matlab code by Sigmund
- Uses modern improvements from COMET-FEniCS tutorials  
- Shows structure evolving during optimization

**What you get**:
- `compliance_history.png` - Optimization progress
- Real-time structure evolution plots
- VTK files for 3D visualization in ParaView

## Main Features

- **Easy to use**: Simple functions for common tasks
- **Well documented**: Each example explains what it does
- **Good graphics**: Built-in plotting with publication-quality figures
- **Export options**: Save results for other software (VTK format)
- **Educational**: Learn topology optimization step-by-step

## Requirements

- Julia 1.6 or newer
- Basic packages (automatically installed): CairoMakie, LinearAlgebra, WriteVTK

## Getting Started

1. Install the package
2. Run the examples to see what it can do
3. Modify parameters to try different cases
4. Use the functions in your own projects

## Help and Contributing

- Report problems or ask questions through GitHub issues
- Suggestions for improvements are welcome
- Documentation and examples can always be improved

---

**Try the examples first** - they show you everything the package can do!