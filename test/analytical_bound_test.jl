
using Revise
using HomoLib: voigt_conductivity, reuss_conductivity,hashin_shtrikman_bounds,maxwell_bounds
                

# Define the conductivities
k1 = 1.0
k2 = 5.0

# Define the volume fraction range for k2
v2_range = range(1e-6, 0.7, length=10)

# Initialize arrays to store results
voigt_results = []
reuss_results = []
hs_lower_results = []
hs_upper_results = []
mw_lower_results = []
mw_upper_results = []

# Compute the bounds for each volume fraction
for v2 in v2_range
    v1 = 1.0 - v2  # Volume fraction of k1
    
    # Voigt and Reuss bounds
    push!(voigt_results, voigt_conductivity([v1, v2], [k1, k2]))
    push!(reuss_results, reuss_conductivity([v1, v2], [k1, k2]))
    
    # Hashin-Shtrikman bounds
    hs = hashin_shtrikman_bounds(k1, v1, k2, v2)
    push!(hs_lower_results, hs.lower)
    push!(hs_upper_results, hs.upper)
    
    # Maxwell bounds (same as Hashin-Shtrikman for two-phase materials)
    mw = maxwell_bounds(k1, v1, k2, v2)
    push!(mw_lower_results, mw.lower)
    push!(mw_upper_results, mw.upper)
end

# thermal bound ok 

using HomoLib: mori_tanaka_elastic,self_consistent_elastic,voigt_elastic,reuss_elastic,hashin_shtrikman_bounds_elastic

# Material properties: Soft matrix (K1=5, G1=2) and stiff inclusions (K2=50, G2=30)
K1, G1 = 5.0, 2.0  # Matrix: Bulk=5 GPa, Shear=2 GPa
K2, G2 = 50.0, 30.0  # Inclusion: Bulk=50 GPa, Shear=30 GPa

# Volume fraction range for inclusion (v2)
v2_range = range(1e-6, 0.7, length=10)

# Initialize arrays to store results
voigt_K = []
voigt_G = []
reuss_K = []
reuss_G = []
hs_K_lower = []
hs_K_upper = []
hs_G_lower = []
hs_G_upper = []

# Compute bounds for each volume fraction
for v2 in v2_range
    v1 = 1.0 - v2  # Volume fraction of matrix
    
    # Voigt and Reuss bounds
    voigt_result = voigt_elastic(K1, G1, v1, K2, G2, v2)
    reuss_result = reuss_elastic(K1, G1, v1, K2, G2, v2)
    push!(voigt_K, voigt_result.K)
    push!(voigt_G, voigt_result.G)
    push!(reuss_K, reuss_result.K)
    push!(reuss_G, reuss_result.G)
    
    # Hashin-Shtrikman bounds
    hs_result = hashin_shtrikman_bounds_elastic(K1, G1, v1, K2, G2, v2)
    push!(hs_K_lower, hs_result.K_lower)
    push!(hs_K_upper, hs_result.K_upper)
    push!(hs_G_lower, hs_result.G_lower)
    push!(hs_G_upper, hs_result.G_upper)
end


