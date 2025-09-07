#' #'Finite Element Analysis of Composite Materials
#'
#' ##'Introduction
#'
#' This document presents a numerical analysis to determine the effective mechanical
#' properties of composite materials. We use the Finite Element Method (FEM)
#' to homogenize a Representative Volume Element (RVE) containing inclusions.
#' The analysis compares a standard elastic model with a model that accounts for
#' contact at the material interfaces.
#'
#' We perform two main studies:
#' 1. **Convergence Study:** To determine how many random realizations are
#'    needed to obtain stable statistical results.
#' 2. **Sensitivity Analysis:** To observe the evolution of the effective stiffness
#'    as a function of the inclusion volume fraction.
#'
#' The results are compared against the theoretical Voigt and Hashin-Shtrikman bounds.
#'
#'--
#' Imports and Configuration
#'--
using Revise
using HomoLib
using CairoMakie
using LinearAlgebra
using Statistics
using Printf
#'
#'--
#'Data Structures
#'--
mutable struct MeshData
    nodes::Matrix{Float64}
    elements::Matrix{Int}
    type_elem::Vector{Int}
    boundary::Vector{Int}
    boundary_inclusion::Vector{Int}
    master::Vector{Int}
    slave::Vector{Int}
end
#'
mutable struct ElemData
    type::Symbol
    int_order::Int
    order::Int
end
#'
#'--
#'Core Analysis Functions
#--
#'
function setup_mesh(; width, height, volume_fraction,
                      n_inclusions, elem_data::ElemData,
                      node_divisions, shape, output_file, voids, rdn, show_gui, to_rotate)

    ind_G, ind_D, ind_B, ind_H, ind_C,
    elements, N_x, N_y, type_elem, _, _ =
        generate_transfinite_plate_with_inclusions(
            width, height, volume_fraction,
            output_file, shape, n_inclusions,
            elem_data.type, elem_data.order,
            node_divisions[1], node_divisions[2];
            voids=voids,
            show_gui=show_gui,
            rdn=rdn,
            to_rotate=to_rotate
        )
    
    nodes = hcat(N_x, N_y)
    boundary = unique(sort(vcat(ind_G, ind_D, ind_B, ind_H)))
    master = vcat(ind_G, ind_B)
    slave = vcat(ind_D, ind_H)

    MeshData(nodes, elements, type_elem, boundary, ind_C, master, slave)
end
#'
function run_standard_elastic_analysis(mesh::MeshData, elem_data::ElemData, E_matrix, ν_matrix, E_inclusion, ν_inclusion)
    material_matrix = ElasticMaterial(2, E=E_matrix, ν=ν_matrix)
    material_inclusion = ElasticMaterial(2, E=E_inclusion, ν=ν_inclusion)

    geometric_data = precompute_geometric_data(
        elem_data.type, elem_data.int_order, material_matrix.dim,
        mesh.elements, mesh.nodes, material_matrix
    )

    K_assembly = assemble_KMF(
        mesh.elements, mesh.nodes,
        [material_matrix, material_inclusion],
        material_matrix.dim, mesh.type_elem, geometric_data
    )
    K_matrix = K_assembly.K

    bc1 = PeriodicHomogenizationBC([1.0 0.0; 0.0 0.0], (mesh.master, mesh.slave), mesh.nodes, material_matrix.dim, [true, true])
    bc2 = PeriodicHomogenizationBC([0.0 0.0; 0.0 1.0], (mesh.master, mesh.slave), mesh.nodes, material_matrix.dim, [true, true])
    bc3 = PeriodicHomogenizationBC([0.0 1/2; 1/2 0.0], (mesh.master, mesh.slave), mesh.nodes, material_matrix.dim, [true, true])

    F_global = zeros(size(K_matrix, 1))
    U_initial = zeros(size(K_matrix, 1))

    U1 = solve!(K_matrix, F_global, U_initial, [bc1])
    U2 = solve!(K_matrix, F_global, U_initial, [bc2])
    U3 = solve!(K_matrix, F_global, U_initial, [bc3])
    
    solver_results = (U = (U1, U2, U3),)

    effective_tensor, _ = compute_effective_property(
        [material_matrix, material_inclusion], mesh.elements, mesh.type_elem,
        solver_results, material_matrix.dim, geometric_data
    )

    return effective_tensor.C
end
#'
#'
function run_contact_elastic_analysis(mesh::MeshData, elem_data::ElemData, E_matrix, ν_matrix, E_inclusion, ν_inclusion)
    material_matrix = ElasticMaterial(2, E=E_matrix, ν=ν_matrix)
    material_inclusion = ElasticMaterial(2, E=E_matrix, ν=ν_matrix)
    material_inclusion.properties[:cracks] = true

    geometric_data = precompute_geometric_data(
        elem_data.type, elem_data.int_order, material_matrix.dim,
        mesh.elements, mesh.nodes, material_matrix
    )

    K_assembly = assemble_KMF(
        mesh.elements, mesh.nodes,
        [material_matrix, ElasticMaterial(2, E=E_inclusion, ν=ν_inclusion)],
        material_matrix.dim, mesh.type_elem, geometric_data
    )
    K_matrix = K_assembly.K

    bc1 = PeriodicHomogenizationBC([1.0 0.0; 0.0 0.0], (mesh.master, mesh.slave), mesh.nodes, material_matrix.dim, [true, true])
    bc2 = PeriodicHomogenizationBC([0.0 0.0; 0.0 1.0], (mesh.master, mesh.slave), mesh.nodes, material_matrix.dim, [true, true])
    bc3 = PeriodicHomogenizationBC([0.0 1/2; 1/2 0.0], (mesh.master, mesh.slave), mesh.nodes, material_matrix.dim, [true, true])

    F_global = zeros(size(K_matrix, 1))
    U_initial = zeros(size(K_matrix, 1))

    U1 = solve!(K_matrix, F_global, U_initial, [bc1])
    U2 = solve!(K_matrix, F_global, U_initial, [bc2])
    U3 = solve!(K_matrix, F_global, U_initial, [bc3])

    solver_results = (U = (U1, U2, U3),)

    effective_tensor, _ = compute_effective_property(
        [material_matrix, material_inclusion], mesh.elements, mesh.type_elem,
        solver_results, material_matrix.dim, geometric_data
    )

    return effective_tensor.C
end
#'
function compute_effective_C(
    volume_fraction::Float64, n_inclusions::Int, element_order::Int, shape::Symbol, element_type::Symbol,
    node_div_inc::Int, node_div_mat::Int; analysis_type::Symbol, material_props)

    is_tri = startswith(string(element_type), "Tri")
    int_order = is_tri ? element_order : element_order * 2
    elem_data = ElemData(element_type, int_order, element_order)

    mesh = setup_mesh(
        width=1.0, height=1.0, volume_fraction=volume_fraction, n_inclusions=n_inclusions,
        elem_data=elem_data, node_divisions=(node_div_inc, node_div_mat), shape=shape,
        output_file="RVE.msh", voids=false, rdn=true, show_gui=false, to_rotate=true
    )

    if analysis_type == :standard
        return run_standard_elastic_analysis(mesh, elem_data, material_props...)
    elseif analysis_type == :contact
        return run_contact_elastic_analysis(mesh, elem_data, material_props...)
    else
        error("Unsupported analysis_type: $analysis_type")
    end
end
#'
function get_mean_effective_C(
    n_realizations::Int, volume_fraction::Float64, n_inclusions::Int, element_order::Int,
    shape::Symbol, element_type::Symbol, node_div_inc::Int, node_div_mat::Int; analysis_type::Symbol, material_props)

    C_tensors = Vector{Matrix{Float64}}(undef, n_realizations)
    for i in 1:n_realizations
        C_tensors[i] = compute_effective_C(
            volume_fraction, n_inclusions, element_order, shape, element_type,
            node_div_inc, node_div_mat; analysis_type=analysis_type, material_props=material_props
        )
    end
    
    mean_C_tensor = sum(C_tensors) / n_realizations
    
    stds = [
        std([c[1, 1] for c in C_tensors]),
        std([c[2, 2] for c in C_tensors]),
        std([c[1, 2] for c in C_tensors]),
        std([c[3, 3] for c in C_tensors])
    ]

    return mean_C_tensor, stds
end
#'
function run_sensitivity_analysis(
    n_samples::Int, n_realizations::Int, n_inclusions::Int, element_order::Int, shape::Symbol,
    element_type::Symbol, node_div_inc::Int, node_div_mat::Int;
    min_vf::Float64, max_vf::Float64, analysis_type::Symbol, material_props)

    vf_values = range(min_vf, max_vf, length=n_samples)
    results = Vector{NamedTuple{(:vf, :mean_C, :stds_C), Tuple{Float64, Matrix{Float64}, Vector{Float64}}}}(undef, n_samples)

    println("Starting Sensitivity Analysis for $(analysis_type)...")
    for i in 1:n_samples
        vf = vf_values[i]
        mean_C, stds_C = get_mean_effective_C(
            n_realizations, vf, n_inclusions, element_order,
            shape, element_type, node_div_inc, node_div_mat;
            analysis_type=analysis_type, material_props=material_props
        )
        results[i] = (vf=vf, mean_C=mean_C, stds_C=stds_C)
    end

    sort!(results, by = x -> x.vf)
    
    #'Convert to the format expected by the plotting functions
    vf_sorted = [r.vf for r in results]
    means_sorted = (
        c1111=[r.mean_C[1,1] for r in results],
        c2222=[r.mean_C[2,2] for r in results],
        c1122=[r.mean_C[1,2] for r in results],
        c1212=[r.mean_C[3,3] for r in results]
    )
    stds_sorted = (
        c1111=[r.stds_C[1] for r in results],
        c2222=[r.stds_C[2] for r in results],
        c1122=[r.stds_C[3] for r in results],
        c1212=[r.stds_C[4] for r in results]
    )

    return vf_sorted, means_sorted, stds_sorted
end
#'
function run_convergence_study(
    max_realizations::Int, n_steps::Int, volume_fraction::Float64,
    n_inclusions::Int, element_order::Int, shape::Symbol, element_type::Symbol,
    node_div_inc::Int, node_div_mat::Int; analysis_type::Symbol, material_props)

    all_C_tensors = [compute_effective_C(
        volume_fraction, n_inclusions, element_order, shape, element_type,
        node_div_inc, node_div_mat; analysis_type=analysis_type, material_props=material_props
    ) for _ in 1:max_realizations]

    step_size = floor(Int, max_realizations / n_steps)
    if step_size == 0; step_size = 1; end
    realization_counts = unique(vcat(collect(1:step_size:max_realizations), max_realizations))
    
    means = [(
        c1111=mean([c[1,1] for c in @view all_C_tensors[1:n]]),
        c2222=mean([c[2,2] for c in @view all_C_tensors[1:n]]),
        c1122=mean([c[1,2] for c in @view all_C_tensors[1:n]]),
        c1212=mean([c[3,3] for c in @view all_C_tensors[1:n]])
    ) for n in realization_counts]

    stds = [(
        c1111= n > 1 ? std([c[1,1] for c in @view all_C_tensors[1:n]]) : 0.0,
        c2222= n > 1 ? std([c[2,2] for c in @view all_C_tensors[1:n]]) : 0.0,
        c1122= n > 1 ? std([c[1,2] for c in @view all_C_tensors[1:n]]) : 0.0,
        c1212= n > 1 ? std([c[3,3] for c in @view all_C_tensors[1:n]]) : 0.0
    ) for n in realization_counts]

    return realization_counts, means, stds
end
#'
function compute_CI(std_devs, n_realizations::Int; confidence::Float64=0.95)
    z_scores = Dict(0.90=>1.645, 0.95=>1.960, 0.99=>2.576)
    z = get(z_scores, confidence, 1.960)
    return z .* std_devs ./ sqrt(n_realizations)
end
#'
function plot_convergence_study(
    realizations, means_std, stds_std, means_contact, stds_contact)
    
    fig = Figure(size=(1200, 1000), fontsize=18)
    ga = fig[1, 1] = GridLayout()
    gb = fig[2, 1] = GridLayout()
    Label(ga[1, 1:2], "Convergence of Mean Stiffness Components", fontsize=24, font=:bold, padding=(0,0,20,0))
    Label(gb[1, 1:2], "Convergence of Standard Deviation", fontsize=24, font=:bold, padding=(0,0,20,0))

    axes_mean = [Axis(ga[2, c], xlabel="Number of Realizations", ylabel="Mean C (GPa)") for c in 1:2]
    append!(axes_mean, [Axis(ga[3, c], xlabel="Number of Realizations", ylabel="Mean C (GPa)") for c in 1:2])
    axes_std = [Axis(gb[2, c], xlabel="Number of Realizations", ylabel="Std Dev (GPa)") for c in 1:2]
    append!(axes_std, [Axis(gb[3, c], xlabel="Number of Realizations", ylabel="Std Dev (GPa)") for c in 1:2])

    components = [:c1111, :c2222, :c1122, :c1212]
    titles = ["C₁₁₁₁", "C₂₂₂₂", "C₁₁₂₂", "C₁₂₁₂"]

    for i in 1:4
        ax_m, ax_s = axes_mean[i], axes_std[i]
        ax_m.title = titles[i]
        ax_s.title = titles[i]

        ax_m.xticks = realizations
        ax_s.xticks = realizations

        ax_m.xticklabelrotation = π/4
        ax_s.xticklabelrotation = π/4

        lines!(ax_m, realizations, [m[i] for m in means_std], color=:blue, label="Standard")
        scatter!(ax_m, realizations, [m[i] for m in means_std], color=:blue)
        lines!(ax_s, realizations, [s[i] for s in stds_std], color=:blue, label="Standard")
        scatter!(ax_s, realizations, [s[i] for s in stds_std], color=:blue)

        lines!(ax_m, realizations, [m[i] for m in means_contact], color=:red, label="Contact")
        scatter!(ax_m, realizations, [m[i] for m in means_contact], color=:red, marker=:utriangle)
        lines!(ax_s, realizations, [s[i] for s in stds_contact], color=:red, label="Contact")
        scatter!(ax_s, realizations, [s[i] for s in stds_contact], color=:red, marker=:utriangle)
    end
    
    axislegend(axes_mean[1], position=:rt)
    return fig
end
#'
function compute_theoretical_bounds(vf_values, C_matrix)
    #'extract effective bulk (κ) and shear (μ) from your 2D C-matrix convention
    κ_m = C_matrix[1,1] - (4.0/3.0)*C_matrix[3,3]
    μ_m = C_matrix[3,3]
    
    #'Extract ν_matrix for Mori-Tanaka
    ν_matrix = (3*κ_m - 2*μ_m) / (2*(3*κ_m + μ_m))

    #'phase 2 = inclusion / void; 
    κ_i = 0.0
    μ_i = 0.0

    safe_inv(x; tol=1e-12) = abs(x) < tol ? sign(x)/tol : 1.0/x

    #'Voigt & Reuss for bulk and shear
    function voigt_reuss(κ1, μ1, κ2, μ2, f2)
        f1 = 1.0 - f2
        κV = f1*κ1 + f2*κ2
        μV = f1*μ1 + f2*μ2
        invκR = (f1/κ1) + (f2/κ2)
        invμR = (f1/μ1) + (f2/μ2)
        κR = invκR == 0.0 ? 0.0 : 1.0 / invκR
        μR = invμR == 0.0 ? 0.0 : 1.0 / invμR
        return κV, μV, κR, μR
    end

    #' Hashin-Shtrikman (returns κ_lower, μ_lower, κ_upper, μ_upper)
    function hashin_shtrikman_safe(κ1, μ1, κ2, μ2, f2; tol=1e-12)
        f1 = 1.0 - f2

        Δκ = κ2 - κ1
        termK1 = 3*κ1 + 4*μ1
        denom_a = safe_inv(Δκ; tol=tol) + 3*f1/termK1
        a = κ1 + (abs(denom_a) < tol ? 0.0 : f2 / denom_a)

        Δκ_b = κ1 - κ2
        termK2 = 3*κ2 + 4*μ2
        denom_b = safe_inv(Δκ_b; tol=tol) + 3*f2/termK2
        b = κ2 + (abs(denom_b) < tol ? 0.0 : f1 / denom_b)

        Δμ = μ2 - μ1
        safe_mu = x -> abs(x) < tol ? tol * sign(x) : x

        denom_alpha = safe_inv(Δμ; tol=tol) + 6*f1*(κ1 + 2*μ1) / (5*safe_mu(μ1)*(3*κ1 + 4*μ1))
        α = μ1 + (abs(denom_alpha) < tol ? 0.0 : f2 / denom_alpha)

        denom_beta = safe_inv(-Δμ; tol=tol) + 6*f2*(κ2 + 2*μ2) / (5*safe_mu(μ2)*(3*κ2 + 4*μ2))
        β = μ2 + (abs(denom_beta) < tol ? 0.0 : f1 / denom_beta)

        κ_lower, κ_upper = min(a,b), max(a,b)
        μ_lower, μ_upper = min(α,β), max(α,β)
        return κ_lower, μ_lower, κ_upper, μ_upper
    end

 
    C1111 = (κ, μ) -> κ + (4.0/3.0)*μ
    C1122 = (κ, μ) -> κ - (2.0/3.0)*μ
    C1212 = (κ, μ) -> μ


    vf = collect(vf_values)
    n = length(vf)
    voigt = Dict(k => zeros(n) for k in (:c1111,:c1122,:c1212))
    reuss = Dict(k => zeros(n) for k in (:c1111,:c1122,:c1212))
    hs_low = Dict(k => zeros(n) for k in (:c1111,:c1122,:c1212))
    hs_up  = Dict(k => zeros(n) for k in (:c1111,:c1122,:c1212))
    mt = Dict(k => zeros(n) for k in (:c1111,:c1122,:c1212))  

    
    κ2 = κ_i
    μ2 = μ_i
    κ1 = κ_m
    μ1 = μ_m

    for (i, f) in enumerate(vf)
        κV, μV, κR, μR = voigt_reuss(κ1, μ1, κ2, μ2, f)
        kll, mll, kuu, muu = hashin_shtrikman_safe(κ1, μ1, κ2, μ2, f)
        
      
        κ_mt = κ1 * (1 - f) / (1 + 3*(1-ν_matrix)/(2*(1-2*ν_matrix)) * f)
        μ_mt = μ1 * (1 - f) / (1 + (15*(1-ν_matrix))/(7-5*ν_matrix) * f)

        voigt[:c1111][i] = C1111(κV, μV)
        voigt[:c1122][i] = C1122(κV, μV)
        voigt[:c1212][i] = C1212(κV, μV)

        reuss[:c1111][i] = C1111(κR, μR)
        reuss[:c1122][i] = C1122(κR, μR)
        reuss[:c1212][i] = C1212(κR, μR)

        hs_low[:c1111][i] = C1111(kll, mll)
        hs_low[:c1122][i] = C1122(kll, mll)
        hs_low[:c1212][i] = C1212(kll, mll)

        hs_up[:c1111][i] = C1111(kuu, muu)
        hs_up[:c1122][i] = C1122(kuu, muu)
        hs_up[:c1212][i] = C1212(kuu, muu)
        
        #'ADDED: Mori-Tanaka values
        mt[:c1111][i] = C1111(κ_mt, μ_mt)
        mt[:c1122][i] = C1122(κ_mt, μ_mt)
        mt[:c1212][i] = C1212(κ_mt, μ_mt)
    end

    return (voigt=voigt, reuss=reuss, hs_lower=hs_low, hs_upper=hs_up, mori_tanaka=mt)  
end
#'
function plot_sensitivity_results(
    vf_values, means_std, stds_std, means_contact, stds_contact, n_realizations,
    E_matrix, ν_matrix, E_inclusion, ν_inclusion)

    # Calculate matrix and inclusion properties for theoretical bounds
    κ_m = E_matrix * ν_matrix / ((1 + ν_matrix) * (1 - 2ν_matrix))
    μ_m = E_matrix / (2 * (1 + ν_matrix))
    κ_i = E_inclusion * ν_inclusion / ((1 + ν_inclusion) * (1 - 2ν_inclusion))
    μ_i = E_inclusion / (2 * (1 + ν_inclusion))

    # Compute theoretical bounds
    hs_lower_c1111 = zeros(length(vf_values))
    voigt_c1111 = zeros(length(vf_values))
    mt_c1111 = zeros(length(vf_values))  # Mori-Tanaka
    
    for (i, f) in enumerate(vf_values)
        # Hashin-Shtrikman lower bound
        hs_bounds = theorical_bound(:elastic, :hashin_shtrikman, κ_m, μ_m, κ_i, μ_i, f)
        hs_lower_c1111[i] = hs_bounds.κ_lower + (4/3) * hs_bounds.μ_lower
        
        # Voigt upper bound
        voigt_bounds = theorical_bound(:elastic, :voigt, κ_m, μ_m, κ_i, μ_i, f)
        voigt_c1111[i] = voigt_bounds.κ + (4/3) * voigt_bounds.μ
        
        # Mori-Tanaka estimate
        κ_mt = κ_m * (1 - f) / (1 + 3*(1-ν_matrix)/(2*(1-2ν_matrix)) * f)
        μ_mt = μ_m * (1 - f) / (1 + (15*(1-ν_matrix))/(7-5ν_matrix) * f)
        mt_c1111[i] = κ_mt + (4/3) * μ_mt
    end

    fig = Figure(size=(1200, 1000), fontsize=18)
    Label(fig[1, 1:2], "Sensitivity Analysis (N = $(n_realizations))", 
          fontsize=24, font=:bold, padding=(0,0,20,0))
    
    # Create all axes explicitly
    axes = [
        Axis(fig[2, 1], xlabel="Volume Fraction (Porosity) ϕ", ylabel="Effective Stiffness (GPa)"),
        Axis(fig[2, 2], xlabel="Volume Fraction (Porosity) ϕ", ylabel="Effective Stiffness (GPa)"),
        Axis(fig[3, 1], xlabel="Volume Fraction (Porosity) ϕ", ylabel="Effective Stiffness (GPa)"),
        Axis(fig[3, 2], xlabel="Volume Fraction (Porosity) ϕ", ylabel="Effective Stiffness (GPa)")
    ]

    components = [:c1111, :c2222, :c1122, :c1212]
    titles = ["C₁₁₁₁", "C₂₂₂₂", "C₁₁₂₂", "C₁₂₁₂"]

    for i in 1:4
        ax = axes[i]
        ax.title = titles[i]
        comp_key = components[i]

        ci_std = compute_CI(stds_std[comp_key], n_realizations)
        ci_contact = compute_CI(stds_contact[comp_key], n_realizations)

        # Plot numerical results
        scatter!(ax, vf_values, means_std[comp_key], marker=:circle, color=:blue)
        errorbars!(ax, vf_values, means_std[comp_key], ci_std, color=:blue)
        
        scatter!(ax, vf_values, means_contact[comp_key], marker=:utriangle, color=:red)
        errorbars!(ax, vf_values, means_contact[comp_key], ci_contact, color=:red)
        
        # Only plot relevant bounds for C1111 and C2222 (using the same bounds for both)
        if i == 1 || i == 2
            lines!(ax, vf_values, hs_lower_c1111, color=:orange, linestyle=:dash, linewidth=2)
            lines!(ax, vf_values, voigt_c1111, color=:green, linestyle=:dash, linewidth=2)
            lines!(ax, vf_values, mt_c1111, color=:brown, linestyle=:dash, linewidth=2)
        end
    end
    
    elements = [
        MarkerElement(marker=:circle, color=:blue, markersize=15),
        LineElement(color=:blue, linestyle=:solid),
        MarkerElement(marker=:utriangle, color=:red, markersize=15),
        LineElement(color=:red, linestyle=:solid),
        LineElement(color=:orange, linestyle=:dash, linewidth=2),
        LineElement(color=:green, linestyle=:dash, linewidth=2),
        LineElement(color=:brown, linestyle=:dash, linewidth=2)
    ]
    
    labels = [
        "Standard Mean", "Standard 95% CI",
        "Contact Mean", "Contact 95% CI",
        "HS Lower Bound", "Voigt Upper Bound", "Mori-Tanaka"
    ]
    
    Legend(fig[4, 1:2], elements, labels, 
           orientation=:horizontal, tellwidth=false, tellheight=true,
           framevisible=false, nbanks=3)

    return fig
end
#'
function plot_constitutive_law(vf_values, means_std, means_contact, target_vf::Float64)
    
    idx = argmin(abs.(vf_values .- target_vf))
    vf_target = vf_values[idx]
    
    C_std = means_std.c1111[idx], means_std.c2222[idx], means_std.c1122[idx], means_std.c1212[idx]
    C_contact = means_contact.c1111[idx], means_contact.c2222[idx], means_contact.c1122[idx], means_contact.c1212[idx]
    
    C_std_matrix = [C_std[1] C_std[3] 0; C_std[3] C_std[2] 0; 0 0 C_std[4]]
    C_contact_matrix = [C_contact[1] C_contact[3] 0; C_contact[3] C_contact[2] 0; 0 0 C_contact[4]]
    
    fig = Figure(size=(1200, 1000), fontsize=18)
    Label(fig[1, 1:2], "Effective Asymmetric Constitutive Law (ϕ = $(round(vf_target*100, digits=1))% Porosity)", 
          fontsize=24, font=:bold, padding=(0,0,20,0))
    
    axes = [
        Axis(fig[2, 1]),
        Axis(fig[2, 2]),
        Axis(fig[3, 1]),
        Axis(fig[3, 2])
    ]

    titles = [
        "σ₁₁ vs ε₁₁ (Uniaxial Loading)", "σ₂₂ vs ε₂₂ (Uniaxial Loading)",
        "σ₁₂ vs 2ε₁₂ (Shear Loading)", "σ₂₂ vs ε₁₁ (Poisson Effect)"
    ]
    xlabels = ["Strain ε₁₁", "Strain ε₂₂", "Shear Strain 2ε₁₂", "Strain ε₁₁"]
    ylabels = ["Stress σ₁₁ (GPa)", "Stress σ₂₂ (GPa)", "Stress σ₁₂ (GPa)", "Stress σ₂₂ (GPa)"]
    
    strain_range = -0.01:0.0002:0.01
    
    #'C1111 plot
    stress_std_11 = C_std_matrix[1, 1] .* strain_range
    stress_contact_11 = [ε >= 0 ? C_contact_matrix[1, 1] * ε : C_std_matrix[1, 1] * ε for ε in strain_range]
    lines!(axes[1], strain_range, stress_std_11, label="Standard", color=:blue)
    lines!(axes[1], strain_range, stress_contact_11, label="Contact", color=:red)

    #'C2222 plot
    stress_std_22 = C_std_matrix[2, 2] .* strain_range
    stress_contact_22 = [ε >= 0 ? C_contact_matrix[2, 2] * ε : C_std_matrix[2, 2] * ε for ε in strain_range]
    lines!(axes[2], strain_range, stress_std_22, color=:blue)
    lines!(axes[2], strain_range, stress_contact_22, color=:red)

    #'C1212 (shear) plot
    stress_std_12 = C_std_matrix[3, 3] .* strain_range
    stress_contact_12 = [ε >= 0 ? C_contact_matrix[3, 3] * ε : C_std_matrix[3, 3] * ε for ε in strain_range]
    lines!(axes[3], strain_range, stress_std_12, color=:blue)
    lines!(axes[3], strain_range, stress_contact_12, color=:red)

    #'C1122 (Poisson) plot
    stress_std_1122 = C_std_matrix[1, 2] .* strain_range
    stress_contact_1122 = [ε >= 0 ? C_contact_matrix[1, 2] * ε : C_std_matrix[1, 2] * ε for ε in strain_range]
    lines!(axes[4], strain_range, stress_std_1122, color=:blue)
    lines!(axes[4], strain_range, stress_contact_1122, color=:red)
    
    for i in 1:4
        axes[i].title = titles[i]
        axes[i].xlabel = xlabels[i]
        axes[i].ylabel = ylabels[i]
        vlines!(axes[i], [0.0], color=:black, linestyle=:dash)
        hlines!(axes[i], [0.0], color=:black, linestyle=:dash)
    end

    axislegend(axes[1], position=:lt)
    return fig
end
#'
#' 
#' ##'Analysis Execution
#'
#' We now define the simulation parameters and run the studies.
#'
#' 
#' Simulation Parameters
n_realizations = 30
n_samples = 5
max_conv_real = 100 
n_conv_steps = 10
#'
#'RVE Parameters
n_inclusions = 10
element_order = 2
shape = :circle
element_type = :Tri6
node_div_inc = 10
node_div_mat = 20
#'
#'Material Properties: Generic Composite with Soft Inclusions
#'Units are in GPa.
E_matrix, ν_matrix = 30.0, 0.3
E_inclusion, ν_inclusion = 1e-10, 1e-10
material_props = (E_matrix, ν_matrix, E_inclusion, ν_inclusion)
#'
#'
#' ###'1. Convergence Study
#'
#' The convergence study is performed for a fixed inclusion volume fraction of 60%.
#' We observe how the mean and standard deviation of the effective stiffness
#' components stabilize as the number of geometric realizations increases.
#'
#'
println("\n Running Convergence Study ")
realization_counts, means_std_conv, stds_std_conv = run_convergence_study(
    max_conv_real, n_conv_steps, volume_fraction, n_inclusions, element_order,
    shape, element_type, node_div_inc, node_div_mat; analysis_type=:standard, material_props=material_props
)
_, means_contact_conv, stds_contact_conv = run_convergence_study(
    max_conv_real, n_conv_steps, volume_fraction, n_inclusions, element_order,
    shape, element_type, node_div_inc, node_div_mat; analysis_type=:contact, material_props=material_props
)
fig_conv = plot_convergence_study(realization_counts, means_std_conv, stds_std_conv, means_contact_conv, stds_contact_conv)
save("convergence_study.png", fig_conv)
fig_conv
#'
#'
#' ###'2. Sensitivity Analysis
#'
#' The sensitivity analysis explores the impact of the volume fraction,
#' varying from 10% to 70%, on the effective stiffness components.
#' Each data point is the average of 20 realizations.
#'
#'
println("\n Running Sensitivity Analysis ")
min_vf, max_vf = 0.1, 0.70
#'
vf_values, means_std, stds_std = run_sensitivity_analysis(
    n_samples, n_realizations, n_inclusions, element_order, shape, element_type,
    node_div_inc, node_div_mat; min_vf=min_vf, max_vf=max_vf, analysis_type=:standard, material_props=material_props
)
_, means_contact, stds_contact = run_sensitivity_analysis(
    n_samples, n_realizations, n_inclusions, element_order, shape, element_type,
    node_div_inc, node_div_mat; min_vf=min_vf, max_vf=max_vf, analysis_type=:contact, material_props=material_props
)
#'
#'Plot with theoretical bounds including Mori-Tanaka
fig_sens = plot_sensitivity_results(
    vf_values, means_std, stds_std, means_contact, stds_contact, n_realizations,
    E_matrix, ν_matrix, E_inclusion, ν_inclusion
)
save("sensitivity_analysis.png", fig_sens)
fig_sens
#'
#'
#' ###'3. Final Comparison and Conclusion
#'
#' Finally, we compare the effective constitutive law for both models
#' at a fixed volume fraction of 60%, reusing the data computed during the sensitivity study.
#'
#'
println("\n Plotting Constitutive Law using Sensitivity Analysis Results")
#'
#'Plot constitutive law using the averaged results from sensitivity analysis
target_vf = 0.4  #'60% volume fraction
fig_constitutive = plot_constitutive_law(vf_values, means_std, means_contact, target_vf)
save("constitutive_law_$(Int(target_vf*100))percent.png", fig_constitutive)
fig_constitutive
#'
#'
#' The following table summarizes the final results.
#'
#'
#'Also create summary table
println("\n  Summary Results for ϕ = $(target_vf*100)% Volume Fraction ")
#'
#'Find the index closest to target porosity
idx = argmin(abs.(vf_values .- target_vf))
mean_C_std_at_target = [means_std.c1111[idx], means_std.c2222[idx], means_std.c1122[idx], means_std.c1212[idx]]
mean_C_contact_at_target = [means_contact.c1111[idx], means_contact.c2222[idx], means_contact.c1122[idx], means_contact.c1212[idx]]
#'
@printf "Component | Standard (GPa) | Contact (GPa) | Difference (%%)\n"
@printf "-|-||-\n"
comps = ["C₁₁₁₁", "C₂₂₂₂", "C₁₁₂₂", "C₁₂₁₂"];
for i in 1:4
    std_val = mean_C_std_at_target[i]
    con_val = mean_C_contact_at_target[i]
    diff = (con_val / std_val - 1) * 100
    @printf "%-9s | %14.4f | %13.4f | %13.2f\n" comps[i] std_val con_val diff
end
#'

