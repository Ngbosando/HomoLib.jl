# Optimisation topologique SIMP

## Objectif

La méthode SIMP pour l’optimisation topologique est présentée ici, comme un exercice pratique avant un cours plus avancé sur le sujet. 

---

## 1. Fondement 


L’implémentation repose sur deux sources :

1. **Sigmund, O. (2001)**: "A 99 line topology optimization code written in Matlab" - *Structural and Multidisciplinary Optimization*
   - Fournit un cadre de base de la méthode SIMP
   - Montre comment faire de l’optimisation topologique avec très peu de code
   - Pose les bases mathématiques de l’optimisation par densité

2. **COMET-FEniCS Tutorial**: [SIMP Topology Optimization](https://comet-fenics.readthedocs.io/en/latest/demo/topology_optimization/simp_topology_optimization.html)
   - Techniques d'implémentation et améliorations
   - Méthodes de filtrage et stratégies de convergence
---

## 2. Notions clés

### 2.1 Fondement mathématique

**Variable de conception**: EDensités des éléments ρₑ ∈ [0,1]

**Interpolation du matériau**: 
```
E(ρₑ) = Eₘᵢₙ + (E₀ - Eₘᵢₙ) × ρₑᵖ
```

**Problème d’optimisation**:
```
minimize:    C = UᵀKU  (souplesse)
subject to:  V/V₀ ≤ Vf (volume max)
             KU = F    (équilibre)
```
---

## 3. Implémentation 

### 3.1 Adptation du code 99-Line de Matlab à HomoLib.jl

#### Version Matlab:
```matlab
% Interpolation du matériau
E = Emin + (E0-Emin)*x.^penal;

% Analyse de sensibilité
dc = -penal*(E0-Emin)*x.^(penal-1).*Ue;
```

#### Version HomoLib.jl:
```julia
# Interpolation avec pénalisation adaptative
E_vec = E * (ρ_min .+ (1 - ρ_min) .* ρ_vec.^p_current)

# Sensibilité avec MEF complet 
dC[e] = -p_current * (1 - ρ_min) * ρ_vec[e]^(p_current - 1) * element_energy
```

### 3.2 Améliorations

#### Pénalisation adaptative
```julia
gray_level = 4 * sum((ρ_vec .- ρ_min) .* (1 .- ρ_vec) .* areas) / total_volume
if gray_level < threshold && iterations_condition
    p_current = min(p_current * (1 + 0.2 * (1 - gray_level)), pmax)
end
```

#### Filtrage des densités
```julia
for i in 1:n_elem
    for j in 1:n_elem
        dist = norm(elem_centers[i] - elem_centers[j])
        weight = max(0, filter_radius - dist)
    end
end
```
---

## 4. Cas : Cantilever Beam (poutre en porte-à-faux)

### Données

**Géométrie**: poutre 6.0 × 1.0
**Charge**: force au coin supérieur droit
**Appuis**: encastrement au coin inférieur gauche
**Contrainte**: volume limité à 40 %

---

## 5. Structure du code

### 5.1 Modular Design 

```julia
# Main learning loop
for iter in 1:max_iter
    # 1. Résolution MEF
    K, F = assemble_system(materials, geometry)
    U = solve_equilibrium(K, F, boundary_conditions)
    
    # 2. Calcul de l’objectif 
    compliance = compute_compliance(U, F)
    
    # 3. Sensibilités
    sensitivities = compute_sensitivities(U, materials, ρ_vec)
    
    # 4. Filtrage
    ρ_filtered = apply_density_filter(ρ_vec, filter_radius)
    
    # 5.  Mise à jour des densités
    ρ_vec = update_densities(ρ_filtered, sensitivities, constraints)
    
    # 6. Visualisation
    visualize_progress(iter, ρ_vec, compliance)
end
```

### 5.2 Visualization
### Figure 1: Convergence
![Convergence history](/docs/images/topology_optimization_cantilever/compliance_history.png)

### Figure 2: Sensibilité
![Final structure](/docs/images/topology_optimization_cantilever/Final_state.png)

```julia

plot_density(nodes, elems, ρ_vec, "Iteration $iter")

plot_compliance_history(compliance_history)

write_vtk_results(iter, nodes, elems, ρ_vec)
```

## 6. Results and Validation

### 6.1 Résultats de la convergence

| Indicateur | Début | Fin |
|--------|---------|-------|
| Souplesse | ~145-50 | ~120-25 | 
| p-value | 1.0 | 4.0 |
| Gray level | ~0.25 | <0.01 |
| Itérations | - | 61 |

### 6.2 Observations

#### Evolution des densités 

- **Redistribution du matériau** following structural logic
- **Forme Treillis** qui apparaissent progressivement
- **Amélioration de la convergence** au fil des itérations

---

## 7.Fonctions principales

```julia

function topology_optimization_CB()
function update_density(ρ, dC, areas, ρ_min, V_f, total_volume)
function plot_density(nodes, connect, ρ, title_str)
function plot_compliance_history(compliance_history)
function write_vtk_results(iter, nodes, elems, ρ_vec)

```
---

## 8. Conclusion: 

Cette implémentation est un **bon point de départ pratique** pour l’optimisation topologique ::

- Compréhension des bases de la méthode SIMP
- Intuition sur le comportement et la convergence des algorithmes
- Savoir-faire pour coder, visualiser et analyser des résultats
- Préparation à explorer des méthodes avancées et à suivre des cours plus poussés

**Prochaine étape**: un enseignement plus théorique et la découverte de nouvelles approches.

---

## Références

1. **Sigmund, O. (2001)**. "A 99 line topology optimization code written in Matlab." *Structural and Multidisciplinary Optimization*, 21(2), 120-127.

2. **COMET-FEniCS Project**. "SIMP Topology Optimization Tutorial." Available at: https://comet-fenics.readthedocs.io/en/latest/demo/topology_optimization/simp_topology_optimization.html

3. **HomoLib.jl Documentation**. Package Julia pour l’homogénéisation et l’analyse EF.