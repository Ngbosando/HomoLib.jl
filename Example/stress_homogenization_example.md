# Exemple : Stress Analysis et Homogénéisation Poroélastique

## 1. Analyse de contraintes autour d’un trou circulaire (Kirsch)

### Configuration
- **Géométrie** : plaque carrée (1×1) avec inclusion circulaire (rayon = 0.2).  
- **Chargement** : traction uniaxiale (ε₁₁ imposé).  
- **Conditions aux limites** :  
  - Déplacement imposé sur le bord droit (ux = 0.01).  
  - Bord gauche bloqué en x, base bloquée en y.  
- **Matériau** : isotrope, E = 10, ν = 1/3.  
- **Maillage** : éléments quadratiques, raffinement local autour du trou.  

La référence analytique est donnée par la **solution de Kirsch** :

σθθ(r=a,θ) = σ∞ (1 - 2 cos(2θ))

---

### Résultats
- Minimum FEM ≈ **-0.983** (théorique -1)  
- Maximum FEM ≈ **3.011** (théorique +3)  

#### Figure
![Comparaison Kirsch FEM](/docs/images/stress_homogenization_example/kirsch_compare.png)  

→ La distribution FEM reproduit fidèlement les pics de traction et compression.  

---

## 2. Homogénéisation poroélastique

### Configuration
- **Géométrie** : cellule carrée avec inclusion circulaire (volume de vide ≈ 20%).  
- **Chargements appliqués** :  
  - ε₁₁ = 1, ε₂₂ = 0, ε₁₂ = 0  
  - ε₁₁ = 0, ε₂₂ = 1, ε₁₂ = 0  
  - ε₁₁ = 0, ε₂₂ = 0, ε₁₂ = 1  
- **Matériau solide** : E = 10, ν = 1/3.  
- **Méthode** : FEM (Quad36, ordre 2, intégration 5).  
- **Vérification énergétique** : W_int ≈ W_ext.  

---

### Résultats numériques
- **Tenseur effectif C** :  
```
C = [
  10.3009   3.9914   ~0
   3.9914  10.3009   ~0
   ~0       ~0       3.0622
]
```

- **Vecteur effectif B** :  
```
B = [0.492, 0.492, ~0]
```

- **Coefficient Biot (numérique)** : solid_biot ≈ 0.026.  

### Comparaison avec modèles analytiques
- **Dilute** : b = 0.60, s_biot = 0.040  
- **Refined** : b = 0.478, s_biot = 0.028  
- **Numérique** : b_eff ≈ 0.492, s_biot ≈ 0.026  

→ Les résultats FEM sont en **très bon accord** avec le modèle raffiné.  

---

### Figures

#### Contraintes locales (f = 0.2)  
![σ₁₁](/docs/images/stress_homogenization_example/σ₁₁_poro.png)  
![σ₂₂](/docs/images/stress_homogenization_example/σ₂₂_poro.png)  
![σ₁₂](/docs/images/stress_homogenization_example/σ₁₂_poro.png)  

#### Déformations locales (f = 0.2)  
![ε₁₁](/docs/images/stress_homogenization_example/ϵ₁₁_poro.png)  
![ε₂₂](/docs/images/stress_homogenization_example/ϵ₂₂_poro.png)  
![ε₁₂](/docs/images/stress_homogenization_example/ϵ₁₂_poro.png)  

---

## Conclusion
- La simulation de la plaque perforée reproduit correctement la solution analytique de Kirsch, avec un facteur de concentration de contraintes ≈ 3 et un minimum ≈ −1.
- L’homogénéisation poroélastique fournit des coefficients effectifs cohérents avec les estimations analytiques.
