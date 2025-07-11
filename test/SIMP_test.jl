using Revise
using HomoLib: material_def,assemble_global_matrix,plaque,
            extract_border_nodes_from_elements,force_computation,
            shape_data,jacobian_data,build_B_matrices,assemble_global_dofs,
            plaque,getBoundaryNodes,BoundaryCondition, solve!,compute_element_stiffness
using CairoMakie
using LinearAlgebra, SparseArrays, Statistics

# Topology optimization of cantilever beam

function topology_optimization_CB()
    # Problem parameters
    b = 6.0    # Width
    h = 1.0    # Height
    lc = 0.1   # Element size
    ν = 0.3    # Poisson's ratio
    E = 1.0    # Young's modulus
    f = -1     # Load
    int_order = 1
    dim = 2

    # Optimization parameters
    p = 3.0       # SIMP exponent
    ρ_min = 1e-3  # Minimum density
    V_f = 0.4     # Volume fraction
    max_iter = 80
    tol = 1e-3

    # Generate mesh 
    nodes2D, elems2D, border_tag = plaque(b, h, lc, 40, 30, "patch2D", 1, :quadrilateral; show_gui=false)
    n_elem = size(elems2D,1)
    
    # Initialize density
    ρ_vec = fill(V_f, n_elem)
    C_prev = Inf
    compliance_history = Float64[]
    
    # Main optimization loop
    for iter in 1:max_iter
       

        # Compute Young's modulus for each element
        E_vec = E * (ρ_min .+ ρ_vec).^p

        # Material definition (linear isotropic elasticity) 

        # need to change most of my function now since all consider unique material
        materials = [material_def([:elastic], 2, :isotropic; E=E_vec[e], ν=ν) for e in 1:n_elem];
        connect_elem_phase = collect(1:n_elem) 
        # force initialisation
        
        # NodalForces = Dict(
        # :right => (fᵥ = nothing, fₜ = force_computation((x, y)  -> f)));
        # BoundaryFace = Dict(
        # :right => (
        #     element_border = [2],
        #     element_type = :Lin2,
        #     dim = 1,
        #     nodes = nodes2D,
        #     int_order = 1
        # ));
        PointForces = [
            (2, [0.0, -1.0]),  # Node 10, Y-direction force
        ]

        # Precompute data
        gauss_data = shape_data(:Quad4, int_order, dim)
        jacobian_cache = jacobian_data(elems2D, nodes2D, gauss_data, dim)
        B_dicts = build_B_matrices(nodes2D, elems2D, materials[1], gauss_data, jacobian_cache)
        Geometric_Data = (gauss_data = gauss_data,jacobian_cache = jacobian_cache,B_dicts = B_dicts )
        # Precompute element volumes
        areas = zeros(n_elem)
        for e in 1:n_elem
            J_data = jacobian_cache[e]
            total_vol = 0.0
            for qp in 1:length(gauss_data.weights)
                detJ, _ = J_data[qp]
                total_vol += abs(detJ) * gauss_data.weights[qp]
            end
            areas[e] = total_vol
        end
        total_volume = sum(areas)

        # Assemble stiffness matrix
        K, F = assemble_global_matrix(elems2D, nodes2D, :Quad4, int_order, materials, dim, connect_elem_phase,Geometric_Data;
                                    PointForces)
                
        # Apply boundary conditions (MBB specific)
        u = zeros(size(K,1))
        
        dof_mask = [true, true]    # Constrain both ux (1st DOF) and uy (2nd DOF)
        values =  [0.0, 0.0]
       
        dirichlet_bc = BoundaryCondition(
        :dirichlet,
        unique(border_tag[:left]),
        dof_mask,
        values,
        2  # dofs_per_node
        );
        # Solve system
        U = solve!(K, F, u, [dirichlet_bc])
        
        # Compute compliance
        C = dot(F, U)
        push!(compliance_history, C)
        
        # Compute sensitivity
        dC = compute_sensitivity(U, materials, elems2D, n_elem, p, ρ_min, ρ_vec, Geometric_Data)
        
        # Update design variables
        ρ_old = copy(ρ_vec)
        ρ_vec = update_density(ρ_old, dC, areas, ρ_min, V_f, total_volume)
        
        # Check convergence
        change = norm(ρ_vec - ρ_old, Inf)
       
        
        # Visualize every 5 iterations
        if iter % 5 == 0
            # println("Iteration $iter")
            # println("Compliance: $C, Change: $change")
            plot_density(nodes2D, elems2D, ρ_vec, "Density at iteration $iter")
        end
        
        if change < tol
            println("Converged at iteration $iter")
            break
        end
        C_prev = C
    end
    
    # Final visualization
    plot_density(nodes2D, elems2D, ρ_vec, "Final Density")
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

# Density visualization

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

# Run optimization
ρ_opt, compliance_hist = topology_optimization_CB()