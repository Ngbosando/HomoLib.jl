using Revise
using HomoLib: material_def,assemble_global_matrix,plaque,
            extract_border_nodes_from_elements,force_computation,
            shape_data,jacobian_data,build_B_matrices,assemble_global_dofs,
            plaque,getBoundaryNodes,BoundaryCondition, solve!,compute_element_stiffness,
            plot_champ_scalaire
using CairoMakie,GeometryBasics
using LinearAlgebra, SparseArrays, Statistics
using WriteVTK

# Topology optimization of cantilever beam

function topology_optimization_CB()
    # Problem parameters
    b = 2.0    # Width
    h = 1.0    # Height
    lc = 0.1   # Element size
    ν = 0.3    # Poisson's ratio
    E = 1.0    # Young's modulus
    f = -1     # Load
    int_order = 1
    dim = 2

    # Optimization parameters
    ρ_min = 1e-3  # Minimum density
    V_f = 0.4     # Volume fraction
    max_iter = 100  # Increased for better convergence
    tol = 1e-3

    # Enhanced SIMP parameters
    niternp = 20       # Number of non-penalized iterations
    pmax = 4.0         # Maximum SIMP exponent
    p_current = 1.0    # Start with non-penalized SIMP
    exponent_update_frequency = 4
    gray_level_threshold = 0.01
    exponent_counter = 0
    old_compliance = Inf

    # Generate mesh 
    nodes2D, elems2D, border_tag = plaque(b, h, lc, 60, 40, "patch2D", 1, :quadrilateral; show_gui=false)
    n_elem = size(elems2D,1)
    
    # Initialize density
    ρ_vec = fill(V_f, n_elem)
    compliance_history = Float64[]
    
    # Precompute element centers and areas
    areas = zeros(n_elem)
    elem_centers = Vector{Vector{Float64}}(undef, n_elem)
    gauss_data = shape_data(:Quad4, int_order, dim)
    jacobian_cache = jacobian_data(elems2D, nodes2D, gauss_data, dim)
    
    for e in 1:n_elem
        J_data = jacobian_cache[e]
        total_vol = 0.0
        center = zeros(2)
        
        for qp in 1:length(gauss_data.weights)
            detJ, _ = J_data[qp]
            total_vol += abs(detJ) * gauss_data.weights[qp]
            
            # Get shape functions at quadrature point
            N = gauss_data.N[qp]
         
            
            # Compute physical coordinate of quadrature point
            phys_coord = sum(N[i] * nodes2D[elems2D[e, i], :] for i in 1:4)
            center .+= phys_coord .* abs(detJ) * gauss_data.weights[qp]
        end
        
        areas[e] = total_vol
        elem_centers[e] = center ./ total_vol
    end
    
    total_volume = sum(areas)
    avg_elem_size = sqrt(mean(areas))
    filter_radius = 1.5 * avg_elem_size
    
    # Main optimization loop
    for iter in 1:max_iter
        # ===== MATERIAL INTERPOLATION =====
        # Enhanced SIMP interpolation (more stable)
        E_vec = E * (ρ_min .+ (1 - ρ_min) .* ρ_vec.^p_current)
        materials = [material_def([:elastic], 2, :isotropic; E=E_vec[e], ν=ν) for e in 1:n_elem]

        # ===== FEM ANALYSIS =====
        # Reuse existing precomputation
        B_dicts = build_B_matrices(nodes2D, elems2D, materials[1], gauss_data, jacobian_cache)
        Geometric_Data = (gauss_data = gauss_data, jacobian_cache = jacobian_cache, B_dicts = B_dicts)
        
        # Assemble stiffness matrix (using your existing function)
        K, F = assemble_global_matrix(elems2D, nodes2D, :Quad4, int_order, materials, dim, collect(1:n_elem), Geometric_Data;
                                    PointForces = [(2, [0.0, -1.0])])
                
        # Apply boundary conditions
        u = zeros(size(K,1))
        dof_mask = [true, true]
        values = [0.0, 0.0]
        dirichlet_bc = BoundaryCondition(
            :dirichlet,
            unique(border_tag[:left]),
            dof_mask,
            values,
            2
        )
        
        # Solve system (using your existing solver)
        U = solve!(K, F, u, [dirichlet_bc])
        
        # Compute compliance
        compliance = dot(F, U)
        push!(compliance_history, compliance)
        
        # ===== ADAPTIVE EXPONENT UPDATE =====
        if iter > 1
            # Compute gray level indicator (0 = binary, 1 = gray)
            gray_level = 4 * sum((ρ_vec .- ρ_min) .* (1 .- ρ_vec) .* areas) / total_volume
            
            # Update exponent only after non-penalized phase
            if iter > niternp
                # Check if we should update exponent
                if (abs(compliance - old_compliance) < gray_level_threshold * compliance_history[1] &&
                    exponent_counter >= exponent_update_frequency)
                    
                    # Adaptive exponent increase
                    p_increment = 0.3^(1 + gray_level/2)
                    p_current = min(p_current * (1 + p_increment), pmax)
                    println("Iter $iter: Adaptive p-increase to $p_current (gray_level=$gray_level)")
                    exponent_counter = 0
                else
                    exponent_counter += 1
                end
            end
        end

        # ===== SENSITIVITY ANALYSIS =====
        dC = zeros(n_elem)
        ref_mat = materials[1]
        global_dofs = assemble_global_dofs(elems2D, ref_mat)
        
        for e in 1:n_elem
            mat = materials[e]
            B_dict = B_dicts[e]
            J_data = jacobian_cache[e]
            elem_dof = global_dofs[e, :]
            uₑ = U[elem_dof]
            
            # Compute element energy (more accurate than stiffness matrix approach)
            element_energy = 0.0
            for qp in 1:length(gauss_data.weights)
                B = B_dict[:strain][qp]
                detJ, _ = J_data[qp]
                
                # Compute strain at quadrature point
                ε = B * uₑ
                σ = mat.tensors[1] * ε
                
                # Energy density (corrected with 0.5 factor)
                element_energy += 0.5 * dot(σ, ε) * abs(detJ) * gauss_data.weights[qp]
            end
            
            # Corrected sensitivity for enhanced SIMP
            dC[e] = -p_current * (1 - ρ_min) * ρ_vec[e]^(p_current - 1) * element_energy / (ρ_min + (1 - ρ_min) * ρ_vec[e]^p_current)
        end
        
        # ===== DENSITY FILTERING =====
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
            
            ρ_filtered[i] = weighted_sum / total_weight
        end
        
        # ===== DENSITY UPDATE =====
        ρ_old = copy(ρ_vec)
        ρ_vec = update_density(ρ_filtered, dC, areas, ρ_min, V_f, total_volume)
        change = norm(ρ_vec - ρ_old, Inf)
        
        # Store compliance for next iteration
        old_compliance = compliance
        
        # ===== VISUALIZATION =====
        if iter % 1 == 0 || iter == 1 || change < tol
            # title_str = "Iter $iter: p=$p_current, C=$(round(compliance, digits=4)), Δρ=$change"
            # fig = plot_density(nodes2D, elems2D, ρ_vec, title_str)
            # display(fig)
            
            # VTK export
            write_vtk_results(iter, nodes2D, elems2D, ρ_vec)
        end
        
        # ===== CONVERGENCE CHECK =====
        if change < tol
            println("Converged at iteration $iter with p=$p_current")
            break
        end
    end
    
    # Final visualization
    fig = plot_density(nodes2D, elems2D, ρ_vec, "Final Density")
    display(fig)
    
    return ρ_vec, compliance_history
end

# Sensitivity analysis
function compute_sensitivity(U, materials, connectivity, n_elem, p, ρ_min, ρ_vec, Geometric_Data)
    dC = zeros(n_elem)
    
    gauss_data = Geometric_Data.gauss_data
    jacobian_cache = Geometric_Data.jacobian_cache
    B_dicts = Geometric_Data.B_dicts
    ref_mat = materials[1]
    global_dofs = assemble_global_dofs(connectivity, ref_mat)
    for e in 1:n_elem
        # Element data
            mat = materials[e]
            B_dict = B_dicts[e]
            J_data = jacobian_cache[e]
            elem_dof = global_dofs[e, :]
        
        # Stiffness computation
            kₑ = compute_element_stiffness(mat, B_dict, J_data, gauss_data)
            uₑ = U[elem_dof]
        
        # # Sensitivity = -p(ρ_min + ρ)⁽ᵖ⁻¹⁾ * uₑᵀ * kₑ * uₑ
        #  dC[e] = -p * (ρ_min + ρ_vec[e])^(p-1) * dot(uₑ, kₑ, ue)
        # Sensitivity = -p * (uₑ' * kₑ * uₑ) / (ρ_min + ρ_vec[e])
        dC[e] = -p * dot(uₑ, kₑ, uₑ) / (ρ_min + ρ_vec[e])
    end
    
    return dC
end

# Density update using Optimality Criteria
function update_density(ρ, dC, areas, ρ_min, V_f, total_volume)
    l1 = 0.0
    l2 = 1e5
    move = 0.2
    η = 0.5

    # Initialize ρ_new to the current densities
    ρ_new = copy(ρ)

    # Bisection method for Lagrange multiplier
    while (l2 - l1) > 1e-4
        lmid = 0.5*(l1+l2)
        ρ_new = max.(ρ_min, min.(1.0, ρ .* (-dC./(lmid*areas)).^η))

        if sum(ρ_new .* areas) > V_f * total_volume
            l1 = lmid
        else
            l2 = lmid
        end
    end

    return ρ_new
end

function plot_density(nodes, connect, ρ, title_str)
    fig = Figure()
    ax = Axis(fig[1, 1], title=title_str, aspect=DataAspect())
    
    # Handle constant density case
    minρ = minimum(ρ)
    maxρ = maximum(ρ)
    colorrange = minρ ≈ maxρ ? (minρ, minρ + 1e-5) : (minρ, maxρ)
    Ne = size(connect,1)
    # Create polygons
    polys = [Point2f.(nodes[connect[e, :], 1], nodes[connect[e, :], 2]) 
             for e in 1:Ne]
    
    # Plot with explicit color range
    poly!(ax, polys, color=ρ, strokecolor=(:black, 0.5), strokewidth=0.5,
          colormap=:viridis, colorrange=colorrange)
    
    Colorbar(fig[1, 2], limits=colorrange, colormap=:viridis)
    
    display(fig)
    return fig
end
function write_vtk_results(iter, nodes, connect, ρ_vec)
    output_dir = "vtk_results"
    !isdir(output_dir) && mkdir(output_dir)
    Ne = size(connect, 1)
    # Convert to MeshCell array with correct quadrilateral type
    cells = [MeshCell(VTKCellTypes.VTK_QUAD, connect[i, :]) for i in 1:Ne]
    
    # Create VTK file
    vtk_file = vtk_grid("$output_dir/iter_$iter", nodes', cells)  # Note the transpose for node matrix
    
    # Add density field
    vtk_file["Density", VTKCellData()] = ρ_vec
    
    # Save file
    outfiles = vtk_save(vtk_file)
    # println("Saved: $(outfiles[1])")
end
# Run optimization
ρ_opt, compliance_hist = topology_optimization_CB();