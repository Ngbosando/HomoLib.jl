# Analyse par éléments finis des matériaux composites avec HomoLib.jl

## Résumé

Ce document présente un cadre numérique pour calculer les propriétés effectives de matériaux composites à l’aide de la méthode des éléments finis (MEF). L’analyse s’appuie sur la bibliothèque HomoLib.jl pour réaliser l’homogénéisation de volume élémentaire représentatifs (VER) contenant des inclusions rectangulaire réparties aléatoirement.

### Points clés
- **Comparaison de modèles**: modèle élastique standard vs. modèle avec interface de contact
- **Analyse de convergence statistique**: pour déterminer le nombre optimal de simulations nécessaires
- **Études de sensibilité**: effet de la fraction volumique sur les propriétés effectives
- **Validation théorique**: comparaison avec les bornes de Voigt, Hashin-Shtrikman et Mori-Tanaka

---

## 1. Introduction
Le composite étudié est constitué d’une microstructure biphasée. Chaque phase est isotrope, linéairement élastique, avec des propriétés différentes. Les inclusions sont supposées parfaitement rectangulaire, et leur position suit une loi de probabilité uniforme, sans chevauchement.

### 1.1 Objectif
Déterminer les propriétés élastiques effectives de matériaux hétérogènes avec inclusions molles (ou vides) par homogénéisation numérique.

### 1.2 Méthodologie
Deux études complémentaires :

1. **Convergence statistique**: nombre minimal de réalisations RVE nécessaires pour obtenir des résultats fiables
2. **Analyse de sensibilité**: lien entre la fraction volumique des inclusions et la rigidité effective

### 1.3 Modèles physiques

#### Modèle élastique
- Parfait collage matrice/inclusion
- Comportement identique en traction et compression
- Elasticité linéaire
- Interfaces supposées parfaites

#### Contact Interface Model
- Réponse différente en traction et en compression
- Plus réaliste pour les matériaux poreux

---

## 2. Implémentation

### 2.1 Dépendances

```julia
using Revise
using HomoLib
using CairoMakie
using LinearAlgebra
using Statistics
using Printf
```

### 2.2 Structures

#### MeshData 
Stocke toutes les infos du maillage:
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

#### ElemData 
Définit les propriétés des éléments:
```julia
mutable struct ElemData
    type::Symbol        # Element type (:Tri3, :Tri6, :Quad4, etc.)
    int_order::Int      # Integration order
    order::Int          # Polynomial order
end
```

### 2.3 Fonctions

- **Génération de maillage** : création d’un VER avec inclusions aléatoires, maillage raffiné, conditions périodiques (GMSH)
- **Analyse élastique standard** : MEF, interfaces parfaites, trois chargements pour caractérisation complète 
- **Analyse avec contact** : Décomposition spectrale et énergétique

---

## 3. Paramètres de simulation

### 3.1 Données statistique

| Paramètre        | Valeur | Description                            |
| ---------------- | ------ | -------------------------------------- |
| `n_realizations` | 30     | Nombre de configurations RVE par point |
| `n_samples`      | 5      | Points de fraction volumique           |
| `max_conv_real`  | 50     | Réalisations max pour convergence      |
| `n_conv_steps`   | 10     | Points d’échantillonnage               |


### 3.2 Configuration VER

| Paramètre         | Valeur    | Description                       |
| ----------------- | --------- | --------------------------------- |
| `volume_fraction` | 0.4 & 0.6 | Porosité cible (40 % et 60 %)     |
| `n_inclusions`    | 10        | Nombre d’inclusions               |
| `element_order`   | 2         | Éléments quadratiques             |
| `shape`           | :square   | Forme des inclusions              |
| `element_type`    | :Tri6     | Triangles à 6 nœuds               |
| `node_div_inc`    | 10        | Raffinement autour des inclusions |
| `node_div_mat`    | 20        | Raffinement matrice               |


### 3.3 Propriétés matériaux

| Matériau  | E (GPa) | ν       | Remarques                   |
| --------- | ------- | ------- | --------------------------- |
| Matrice   | 30.0    | 0.3     | Matériau de base            |
| Inclusion | 1×10⁻¹⁰ | 1×10⁻¹⁰ | Rigidité quasi nulle (vide) |

---

## 4. Résultats et analyse

### 4.1 Convergence

- Moyennes stables dès 15–20 réalisations
- Écarts-types convergent au même rythme
- Le modèle contact montre un peu plus de variabilité

### 4.2 Sensitivity Analysis Results

- Rigidité diminue quasi linéairement avec la porosité
- Accord avec Mori-Tanaka en fraction volumique moyen (0.2-0.4)

### 4.3 Constitutive Behavior Comparison

- **À 40 % porosité** : réponse plus molle en compression, rigidité cisaillement réduite (~5 %)

- **À 55 % porosité** : ouverture d’interface adoucit la traction (~2.6 %), mais provoque un renforcement cisaillement >7 %

### 4.4 Résumés
### 4.4.1 40% porosité
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

## 5. Comparaison avec les bornes théoriques

### 5.1 Implémentation 

1. **Voigt (borne sup)**
2. **Hashin-Shtrikman**
3. **Mori-Tanaka**

### 5.2 Résultats

The numerical results successfully:
- Restent en dessous de la borne de Voigt (logique)
- Proches de Hashin-Shtrikman, à certaines porosités.
- Excellente cohérence avec Mori-Tanaka (fractions modérées).
---

## 6. Visualisations


### Figure 1: Convergence
![Convergence Study: Stiffness evolution and standard deviation convergence](/docs/images/sensitivity_analysis/convergence_study.png)
- Haut : évolution de la rigidité moyenne
- Bas : convergence des écarts-types
- Comparaison entre modèle standard (bleu) et modèle contact (rouge)

### Figure 2: Sensibilité
![Sensitivity Analysis: Stiffness components with 95% confidence intervals](/docs/images/sensitivity_analysis/sensitivity_analysis.png)
- Quatre sous-graphes pour chaque composant de rigidité
- Barres d’erreur avec intervalles de confiance à 95 %
- Bornes théoriques tracées en référence

### Figure 3: Loi de comportement
![Constitutive Law: Stress-strain relationships for 40% porosity](/docs/images/sensitivity_analysis/constitutive_law_40percent.png)
![Constitutive Law: Stress-strain relationships for 55% porosity](/docs/images/sensitivity_analysis/constitutive_law_60percent.png)
- Courbes contrainte–déformation pour porosités 40 % et 55 %
- Mise en évidence de l’asymétrie
- Cas étudiés : traction, cisaillement et effet de Poisson

---

## 7. Conclusions

1. 30 réalisations suffisent pour une bonne stabilité statistique
2. Le modèle contact modifie la rigidité de 5–10 % selon la porosité
depending on the volume fraction consider
3. Les résultats numériques collent bien aux bornes théoriques
4. L’asymétrie traction/compression est correctement capturée

---

## 8. Utilisation du code

### Installation
```julia
using Pkg
Pkg.add("HomoLib")
```

### Exemple minimal
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

### Fichiers générés
- `convergence_study.png`: Statistical convergence analysis
- `sensitivity_analysis.png`: Volume fraction sensitivity
- `constitutive_law_$X$percent.png`: Asymmetric stress-strain behavior

---

## Références

For more information on the theoretical background and implementation details, please refer to:
- HomoLib.jl documentation
- Yvonnet, Computational Homogenization of Heterogeneous Materials with Finite Elements (2019), Periodic BCs, Fine Mesh
- C. Miehe and M. Lambrecht. Algorithms for computation of stresses and elasticity moduli in terms of seth-hill’s family of generalized strain tensors.
---

