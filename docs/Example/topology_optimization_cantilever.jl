### topology_optimization_cb.jl
# # Topology Optimization of a 2D Cantilever Beam
#
# This example demonstrates how to perform topology optimization using the
# SIMP (Solid Isotropic Material with Penalization) method on a cantilever beam.
# The goal is to minimize compliance under a volume constraint.

using HomoLib
using CairoMakie, GeometryBasics
using LinearAlgebra, SparseArrays, Statistics
using WriteVTK

# Helper Functions -
function get_center_node(nodes, boundary_nodes)
    boundary_coords = nodes[boundary_nodes, :]
    center_coord = mean(boundary_coords, dims=1)
    distances = [norm(coord' - center_coord) for coord in eachrow(boundary_coords)]
    return boundary_nodes[argmin(distances)]
end

function update_density(ρ, dC, areas, ρ_min, V_f, total_volume; move=0.2, η=0.5)
    n_elem = length(ρ)
    ρ_new = similar(ρ)
    l1, l2 = 0.0, 1e9

    while (l2 - l1) > 1e-6
        λ = 0.5*(l1 + l2)
        base = -dC ./ (λ .* areas)
        base_clamped = max.(base, 0.0)
        ρ_candidate = clamp.(ρ .* base_clamped .^ η, ρ_min, 1.0)
        vol_frac = sum(ρ_candidate .* areas) / total_volume

        vol_frac > V_f ? (l1 = λ) : (l2 = λ)
    end

    λ = 0.5*(l1 + l2)
    base = -dC ./ (λ .* areas)
    base_clamped = max.(base, 0.0)
    ρ_candidate = clamp.(ρ .* base_clamped .^ η, ρ_min, 1.0)

    for e in 1:n_elem
        ρ_new[e] = clamp(ρ_candidate[e], ρ[e] - move, ρ[e] + move)
    end
    return ρ_new
end

function plot_density(nodes, connect, ρ, title_str)
    fig = Figure()
    ax = Axis(fig[1, 1], title=title_str, aspect=DataAspect())
    minρ, maxρ = extrema(ρ)
    colorrange = minρ ≈ maxρ ? (minρ, minρ + 1e-5) : (minρ, maxρ)
    Ne = size(connect,1)
    polys = [Point2f.(nodes[connect[e, 1:4], 1], nodes[connect[e, 1:4], 2]) for e in 1:Ne]
    poly!(ax, polys, color=ρ, strokecolor=(:black, 0.5), strokewidth=0.5,
        colormap=:viridis, colorrange=colorrange)
    Colorbar(fig[1, 2], limits=colorrange, colormap=:viridis)
    return fig
end

function plot_compliance_history(compliance_history)
    fig = Figure(size=(1000, 800))
    ax = Axis(fig[1, 1],
        title="Compliance History",
        xlabel="Iteration",
        ylabel="Compliance",
        titlesize=24,
        xlabelsize=20,
        ylabelsize=20
    )
    iterations = 1:length(compliance_history)
    lines!(ax, iterations, compliance_history, color=:blue, linewidth=3, label="Compliance")
    scatter!(ax, iterations, compliance_history, color=:red, markersize=12, label="Iterations")
    hlines!(ax, [minimum(compliance_history)], color=:green, linestyle=:dash, linewidth=2, label="Minimum")
    ax.xgridvisible = true; ax.ygridvisible = true
    ax.xgridstyle = :dash; ax.ygridstyle = :dash
    ax.xgridwidth = 0.5; ax.ygridwidth = 0.5
    comp_range = maximum(compliance_history) - minimum(compliance_history)
    if comp_range > 1000
        ax.yscale = Makie.log10
        ax.ylabel = "Log Compliance"
    end
    Legend(fig[1, 2], ax, "Legend"; tellwidth=false, tellheight=false, framevisible=true, backgroundcolor=:white, labelsize=18)
    display(fig)
    save("compliance_history.png", fig)
    return fig
end

function write_vtk_results(iter, nodes, connect, ρ_vec)
    output_dir = "vtk_results"
    !isdir(output_dir) && mkdir(output_dir)
    n_nodes = size(connect, 2)

    cell_type = if n_nodes == 4
        VTKCellTypes.VTK_QUAD
    elseif n_nodes == 9
        VTKCellTypes.VTK_BIQUADRATIC_QUAD
    elseif n_nodes in (16, 25, 36)
        VTKCellTypes.VTK_LAGRANGE_QUADRILATERAL
    elseif n_nodes == 3
        VTKCellTypes.VTK_TRIANGLE
    elseif n_nodes == 6
        VTKCellTypes.VTK_QUADRATIC_TRIANGLE
    elseif n_nodes in (10, 15, 21)
        VTKCellTypes.VTK_LAGRANGE_TRIANGLE
    else
        error("Unsupported element type with $n_nodes nodes")
    end

    Ne = size(connect,1)
    cells = [MeshCell(cell_type, connect[i, :]) for i in 1:Ne]
    vtk_file = vtk_grid("$output_dir/iter_$iter", nodes', cells)
    vtk_file["Density", VTKCellData()] = ρ_vec
    vtk_save(vtk_file)
end

function precompute_element_data(nodes, elems, element_type, int_order, dim)
    gauss_data = compute_shape_function_data(element_type, int_order, 2, size(elems, 2))
    jacobian_cache = compute_element_jacobian_data(elems, nodes, gauss_data, Val(dim))
    n_elem = size(elems,1)
    areas = zeros(n_elem)
    elem_centers = Vector{Vector{Float64}}(undef, n_elem)

    for e in 1:n_elem
        total_vol = 0.0
        center = zeros(2)
        for qp in 1:length(gauss_data.weights)
            detJ, _ = jacobian_cache[qp,e]
            total_vol += abs(detJ) * gauss_data.weights[qp]
            N = gauss_data.N[qp]
            Nn = length(N)
            phys_coord = sum(N[i]*nodes[elems[e,i],:] for i in 1:Nn)
            center .+= phys_coord .* abs(detJ) * gauss_data.weights[qp]
        end
        areas[e] = total_vol
        elem_centers[e] = center ./ total_vol
    end
    return areas, elem_centers, gauss_data, jacobian_cache, sum(areas)
end

function fill_mech_dofs!(out::AbstractVector{<:Integer},
                         conn::AbstractVector{<:Integer},
                         dim::Integer)
    NN = length(conn)
    d  = Int(dim)
    @inbounds for a in 1:NN
        node = Int(conn[a])         # 1-based node id
        base = (node - 1) * d       # global node-wise base dof
        i0   = (a - 1) * d          # element-local offset
        out[i0 + 1] = base + 1
        if d >= 2; out[i0 + 2] = base + 2; end
        if d >= 3; out[i0 + 3] = base + 3; end
    end
    return nothing
end
# Main Optimization Function 
function topology_optimization_CB()
    # Problem Parameters
    b = 6.0    # Width
    h = 1.0    # Height
    lc = 0.1   # Element size
    ν = 0.3    # Poisson's ratio
    E = 1.0    # Young's modulus
    dim = 2
    order = 2
    element_type = :Quad9
    int_order = order * dim

    # Optimization Parameters
    ρ_min = 1e-3
    V_f = 0.4
    max_iter = 200
    tol = 1e-4
    niternp = 20
    pmax = 4.0
    p_current = 1.0
    exponent_update_frequency = 4
    gray_level_threshold = 0.01

    # Generate Mesh
    nodes, elems, border_tag = plaque(b, h, lc, 70, 30, "patch2D", order, :quadrilateral; show_gui=false)
    n_elem = size(elems,1)
    top_middle_node_idx = argmin((nodes[:, 1].- b/2).^2 + (nodes[:, 2].- h).^2)

    # Precompute Element Data (areas/centers for filter)
    areas, elem_centers, gauss_data, jacobian_cache, total_volume =
        precompute_element_data(nodes, elems, element_type, int_order, dim)
    avg_elem_size = sqrt(mean(areas))
    filter_radius = 1.5 * avg_elem_size

    # Initialize Optimization
    ρ_vec = fill(V_f, n_elem)
    compliance_history = Float64[]
    exponent_counter = 0
    old_compliance = Inf

    #  FEM Precompute (new convention) 
    # NOTE: use the new geometric precompute once; re-use every iteration.
    #       For SIMP, only C changes via material interpolation.
    Geometric_Data = precompute_geometric_data(
        element_type, int_order, dim,
        elems, nodes, ElasticMaterial(2, E=E, ν=ν)
    )

    for iter in 1:max_iter
        # Material Interpolation
        E_vec = E * (ρ_min .+ (1 - ρ_min) .* ρ_vec.^p_current)
        materials = [ElasticMaterial(2, E=E_vec[e], ν=ν) for e in 1:n_elem]

        # FEM Analysis (new assemble_KMF API & types)
        PointForces = Dict(top_middle_node_idx => [0.0, -1.0])  # NEW: Dict node=>vector(dofs_per_node)
        K,_,F = assemble_KMF(
            elems, nodes, materials, dim, collect(1:n_elem), Geometric_Data;
            compute_mass=false, NodalForces=nothing, PointForces=PointForces, BoundaryFace=nothing
        )
       

        # Boundary Conditions
        u = zeros(size(K,1))
        bottom_left_node_idx  = argmin(nodes[:, 1].^2 + nodes[:, 2].^2)
        bottom_right_node_idx = argmin((nodes[:, 1].- b).^2 + nodes[:, 2].^2)
        bc1 = DirichletBC([bottom_left_node_idx], [false,true], [0.0,0.0], 2)
        bc2 = DirichletBC([bottom_right_node_idx], [true,true], [0.0,0.0], 2)
        U = solve!(K, F, u, [bc1, bc2])

        # Compliance
        compliance = dot(F, U)
        push!(compliance_history, compliance)

        # Adaptive Exponent Update
        if iter > 1
            gray_level = 4 * sum((ρ_vec .- ρ_min) .* (1 .- ρ_vec) .* areas) / total_volume
            if iter > niternp &&
               abs(compliance - old_compliance) < gray_level_threshold * compliance_history[1] &&
               exponent_counter >= exponent_update_frequency
                p_increment = 0.2 * (1 - gray_level)
                p_current = min(p_current * (1 + p_increment), pmax)
                exponent_counter = 0
            else
                exponent_counter += 1
            end
        end

        #  Sensitivity Analysis (no assemble_global_dofs) 
        dC = zeros(n_elem)

        # Reuse the already built geometric data
        B_all  = Geometric_Data.differential_operator_matrices
        J_all  = Geometric_Data.jacobian_transformation_data
        shp    = Geometric_Data.shape_function_data

        # Preallocate a large enough DOF buffer and slice per element
        maxNN = size(elems, 2)        # nodes per element for this mesh
        gUbuf = Vector{Int}(undef, dim*maxNN)

        for e in 1:n_elem
            conn = @view elems[e, :]
            NN = length(conn)

            # Fill node-wise mechanical DOFs for this element (slice to actual length)
            gU = @view gUbuf[1:dim*NN]
            fill_mech_dofs!(gU, conn, dim)

            uₑ = @view U[gU]
            element_energy = 0.0

            @inbounds for qp in eachindex(shp.weights)
                B    = B_all[e].voigt_gradient_operator[:, :, qp]    
                detJ = (J_all[qp, e])[1]
                w    = shp.weights[qp]

                ε = B * uₑ
                σ = materials[e].C * ε      # direct access to element C
                element_energy += dot(σ, ε) * abs(detJ) * w
            end

            dC[e] = -p_current * (1 - ρ_min) * ρ_vec[e]^(p_current - 1) * element_energy
        end

        # Density Filtering
        ρ_filtered = similar(ρ_vec)
        for i in 1:n_elem
            total_weight = 0.0
            weighted_sum = 0.0
            for j in 1:n_elem
                dist = norm(elem_centers[i] - elem_centers[j])
                weight = max(0, filter_radius - dist)
                total_weight += weight
                weighted_sum += weight * ρ_vec[j]
            end
            ρ_filtered[i] = weighted_sum / max(total_weight, eps())
        end

        # Density Update
        ρ_old = copy(ρ_vec)
        ρ_vec = update_density(ρ_filtered, dC, areas, ρ_min, V_f, total_volume)
        change = norm(ρ_vec - ρ_old, Inf)
        old_compliance = compliance

        # Visualization / Output
        if iter % 20 == 0 || iter == 1 || change < tol
            title_str = "Iter $iter: p=$p_current, C=$(round(compliance, digits=4)), Δρ=$change"
            fig = plot_density(nodes, elems, ρ_vec, title_str)
            # display(fig)
            write_vtk_results(iter, nodes, elems, ρ_vec)
        end

        # Convergence
        if change < tol
            println("Converged at iteration $iter with p=$p_current")
            break
        end
    end

    fig = plot_density(nodes, elems, ρ_vec, "Final Density")
    # display(fig)
    return ρ_vec, compliance_history
end

ρ_opt, compliance_hist = topology_optimization_CB();



