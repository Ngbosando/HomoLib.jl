# Finite Element Analysis of Composite Materials using HomoLib.jl

## Executive Summary

This document presents a comprehensive numerical analysis framework for determining the effective mechanical properties of composite materials using the Finite Element Method (FEM). The analysis employs the HomoLib.jl package to perform homogenization of Representative Volume Elements (RVEs) containing randomly distributed square inclusions.

### Key Features
- **Model Comparison**: Standard elastic model vs. contact interface model
- **Statistical Convergence Analysis**: Determines optimal number of realizations
- **Sensitivity Studies**: Explores volume fraction effects on effective properties
- **Theoretical Validation**: Comparison with Voigt, Hashin-Shtrikman, and Mori-Tanaka bounds

---

## 1. Introduction
The composite is characterized by a 2-phase micro-structure. Each phase is linear elastic isotropic
with different coefficients for each phase. The inclusions are assumed perfectly square, and their
position is distributed according to a uniform probability law with the exclusion of superposition
of inclusions.

### 1.1 Objective
The primary goal is to characterize the effective elastic properties of heterogeneous materials with soft inclusions (or voids) through computational homogenization. This analysis is particularly relevant for:
- Porous materials
- Particle-reinforced composites
- Materials with damage or voids
- Soft inclusion composites

### 1.2 Methodology Overview
We perform two complementary studies:

1. **Convergence Study**: Establishes the minimum number of random RVE realizations required for statistically stable results
2. **Sensitivity Analysis**: Investigates the relationship between inclusion volume fraction and effective stiffness

### 1.3 Physical Models

#### Standard Elastic Model
- Assumes perfect bonding between matrix and inclusions
- Both tension and compression treated identically
- Linear elastic behavior throughout
- Interface are consider perfect

#### Contact Interface Model
- Accounts for potential debonding at matrix-inclusion interfaces
- Asymmetric response in tension vs compression
- More realistic for materials voids

---

## 2. Implementation

### 2.1 Dependencies and Configuration

```julia
using Revise
using HomoLib
using CairoMakie
using LinearAlgebra
using Statistics
using Printf
```

### 2.2 Data Structures

#### MeshData Structure
Stores all mesh-related information:
```julia
mutable struct MeshData
    nodes::Matrix{Float64}           # Nodal coordinates [x, y]
    elements::Matrix{Int}            # Element connectivity
    type_elem::Vector{Int}           # Element type identifier
    boundary::Vector{Int}            # Boundary node indices
    boundary_inclusion::Vector{Int}  # Inclusion boundary nodes
    master::Vector{Int}              # Master nodes for periodic BC
    slave::Vector{Int}               # Slave nodes for periodic BC
end
```

#### ElemData Structure
Defines element properties:
```julia
mutable struct ElemData
    type::Symbol        # Element type (:Tri3, :Tri6, :Quad4, etc.)
    int_order::Int      # Integration order
    order::Int          # Polynomial order
end
```

### 2.3 Core Analysis Functions

#### Mesh Generation
The `setup_mesh` function creates an RVE with randomly distributed inclusions:
- Generates transfinite mesh with specified refinement
- Places inclusions randomly while avoiding overlap
- Applies periodic boundary conditions for homogenization

#### Standard Elastic Analysis
Implements classical FEM with:
- Symmetric stiffness response
- Perfect interface assumption
- Three loading cases for full tensor characterization

#### Contact Interface Analysis
Model featuring:
- Interface elements at inclusion boundaries
- Asymmetric stiffness in tension/compression

---

## 3. Simulation Parameters

### 3.1 Numerical Settings

| Parameter | Value | Description |
|-----------|-------|-------------|
| `n_realizations` | 30 | Number of random RVE configurations per data point |
| `n_samples` | 5 | Number of volume fraction points for sensitivity |
| `max_conv_real` | 50 | Maximum realizations for convergence study |
| `n_conv_steps` | 10 | Sampling points in convergence study |

### 3.2 RVE Configuration

| Parameter | Value | Description |
|-----------|-------|-------------|
| `volume_fraction` | 0.4 & 0.6 | Target porosity (40% and 60%) for detailed analysis |
| `n_inclusions` | 10 | Number of square inclusions per RVE |
| `element_order` | 2 | Quadratic elements (6-node triangles) |
| `shape` | :square | Inclusion geometry |
| `element_type` | :Tri6 | 6-node triangular elements |
| `node_div_inc` | 10 | Mesh divisions around inclusions |
| `node_div_mat` | 20 | Mesh divisions in matrix |

### 3.3 Material Properties

| Material | Young's Modulus (GPa) | Poisson's Ratio | Notes |
|----------|----------------------|-----------------|-------|
| Matrix | 30.0 | 0.3 | Base material |
| Inclusion | 1×10⁻¹⁰ | 1×10⁻¹⁰ | Near-zero stiffness (void-like) |

---

## 4. Results and Analysis

### 4.1 Convergence Study Results

The convergence study demonstrates that:
- **Mean values stabilize** after approximately 15-20 realizations
- **Standard deviations converge** at the same, requiring 15+ realizations for stability
- **Contact model** shows slightly higher variability than standard model

### 4.2 Sensitivity Analysis Results

The sensitivity analysis reveals:
- **Linear degradation** of stiffness with increasing porosity
- Results fall within theoretical bounds (Hashin-Shtrikman)
- Good agreement with Mori-Tanaka model for spherical voids

### 4.3 Constitutive Behavior Comparison
The Contact model incorporates interface mechanics that cause its stiffness predictions to differ from the Standard model, with the effects depending on porosity and load type.

- At 40% Porosity : The model predicts a softer response under compression, and reduces shear and off-axis stiffness by 4-5% 

-At 55% Porosity : The main differences occur under tension. The Contact model shows that interface opening softens the material (~2.6%), but the resulting structural rearrangement causes a significant shear stiffening of over 7%.

### 4.4 Summary Table 
### 4.4.1 40% Porosity
| Component | Standard (GPa) | Contact (GPa) | Difference (%) |
|-----------|---------------|---------------|----------------|
| C₁₁₁₁ | 10.8714 | 10.8899 | +0.17 |
| C₂₂₂₂ | 11.0928 | 10.5339 | -5.04 |
| C₁₁₂₂ | 4.0507 | 3.8834 | -4.13 |
| C₁₂₁₂ | 2.4894 | 2.3707 | -4.77 |
### 4.4.2 55% Porosity
| Component | Standard (GPa) | Contact (GPa) | Difference (%) |
|-----------|---------------|---------------|----------------|
| C₁₁₁₁ | 10.1334 | 9.8662 | -2.64 |
| C₂₂₂₂ | 9.9611 | 10.0315 | +0.71 |
| C₁₁₂₂ | 3.5346 | 3.5707 | +1.02 |
| C₁₂₁₂ | 1.9675 | 2.1185 | +7.67 |
---

## 5. Theoretical Bounds Comparison

### 5.1 Implemented Bounds

1. **Voigt (Upper) Bound**: assumes uniform strain
2. **Hashin-Shtrikman Bounds**: Tighter bounds based on variational principles
3. **Mori-Tanaka Model**: Mean-field approximation for dilute inclusions

### 5.2 Validation Results

The numerical results successfully:
- Remain within Voigt bounds (as required theoretically)
- Fall close to Hashin-Shtrikman bounds 
- Show excellent agreement with Mori-Tanaka for moderate volume fractions

---

## 6. Visualization Outputs

The analysis generates three key figures:

### Figure 1: Convergence Study
![Convergence Study: Stiffness evolution and standard deviation convergence](/docs/images/sensitivity_analysis/convergence_study.png)
- Top row: Evolution of mean stiffness components
- Bottom row: Convergence of standard deviations
- Comparison between standard (blue) and contact (red) models

### Figure 2: Sensitivity Analysis
![Sensitivity Analysis: Stiffness components with 95% confidence intervals](/docs/images/sensitivity_analysis/sensitivity_analysis.png)
- Four subplots for each stiffness component
- Error bars showing 95% confidence intervals
- Theoretical bounds overlaid for validation

### Figure 3: Constitutive Law
![Constitutive Law: Stress-strain relationships for 40% porosity](/docs/images/sensitivity_analysis/constitutive_law_40percent.png)
![Constitutive Law: Stress-strain relationships for 55% porosity](/docs/images/sensitivity_analysis/constitutive_law_60percent.png)
- Stress-strain relationships for 40% & 55% porosity
- Demonstrates asymmetric behavior in contact model
- Includes uniaxial, shear, and Poisson effects

---

## 7. Conclusions

1. **Statistical Convergence**: 30 realizations provide adequate statistical stability for engineering purposes
2. **Model Comparison**: Contact model predicts 5-10% reduction/increase in stiffness versus standard model
depending on the volume fraction consider
3. **Theoretical Agreement**: Numerical results align well with analytical bounds
4. **Asymmetric Response**: Contact model successfully captures tension-compression asymmetry

---

## 8. Code Repository and Usage

### Installation
```julia
using Pkg
Pkg.add("HomoLib")
```

### Basic Usage Example
```julia
# Define material properties
E_matrix, ν_matrix = 30.0, 0.3
E_inclusion, ν_inclusion = 1e-10, 1e-10

# Run analysis
C_effective = compute_effective_C(
    volume_fraction, n_inclusions, element_order,
    shape, element_type, node_div_inc, node_div_mat;
    analysis_type=:contact,
    material_props=(E_matrix, ν_matrix, E_inclusion, ν_inclusion)
)
```

### Output Files
- `convergence_study.png`: Statistical convergence analysis
- `sensitivity_analysis.png`: Volume fraction sensitivity
- `constitutive_law_$X$percent.png`: Asymmetric stress-strain behavior

---

## References

For more information on the theoretical background and implementation details, please refer to:
- HomoLib.jl documentation
- Yvonnet, Computational Homogenization of Heterogeneous Materials with Finite Elements (2019), Periodic BCs, Fine Mesh
- C. Miehe and M. Lambrecht. Algorithms for computation of stresses and elasticity moduli in terms of seth-hill’s family of generalized strain tensors.
---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*