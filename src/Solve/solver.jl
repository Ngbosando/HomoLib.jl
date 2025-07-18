# =============================================
# Data Structures
# =============================================
    """
        BoundaryCondition

        Standard boundary condition for finite element analysis.

        # Fields
        - "type::Symbol": :dirichlet, :neumann, or :periodic
        - "nodes_dof::Union{Vector{Int}, Tuple{Vector{Int}, Vector{Int}}}":
        - For Dirichlet/Neumann: Node IDs
        - For Periodic: (master_nodes, slave_nodes)
        - "dof_mask::Union{Vector{Bool}, Nothing}": Boolean mask for constrained DOFs
        - "macro_matrix::Union{Vector{Float64}, Nothing}": Prescribed values
        - "dofs_per_node::Int": Degrees of freedom per node

        # Examples
        ```julia
        # Fix x-direction for nodes 1-5
        bc = BoundaryCondition(:dirichlet, 1:5, [true, false], [0.0], 2)

        # Periodic constraints between left and right faces
        bc = BoundaryCondition(:periodic, (left_nodes, right_nodes), nothing, nothing, 3)
    """
    struct BoundaryCondition
        type::Symbol           # :dirichlet, :neumann, :periodic
        nodes_dof::Union{Vector{Int}, Tuple{Vector{Int}, Vector{Int}}}
        dof_mask::Union{Vector{Bool}, Nothing}        # True/false mask for DOF constraints
        macro_matrix::Union{Vector{Float64}, Nothing} # Imposed values (Dirichlet/Neumann)
        dofs_per_node::Int
    end

    function Base.show(io::IO, bc::BoundaryCondition)
        println(io, "BoundaryCondition(type = $(bc.type),")
        println(io, "  nodes_dof = $(bc.nodes_dof),")
        println(io, "  dofs_per_node = $(bc.dofs_per_node),")
        println(io, "  dof_mask = ", bc.dof_mask === nothing ? "nothing" : bc.dof_mask)
        println(io, "  macro_matrix = ", bc.macro_matrix === nothing ? "nothing" : bc.macro_matrix, ")")
    end
    """
        HomogenizationBC

        Specialized boundary condition for computational homogenization.

        # Fields
        - "type::Symbol": :dirichlet, :neumann, or :periodic
        - "macro_matrix": Macroscopic strain/stress tensor
        - "node_pairs::Union{Tuple{Vector{Int}, Vector{Int}}, Nothing}": (master_nodes, slave_nodes)
        - "nodes_dof::Union{Vector{Int}, Nothing}": Optional direct DOF specification
        - "dof_mask::Union{Vector{Bool}, Nothing}": Boolean mask for constrained DOFs
        - "coords::Matrix{Float64}": Nodal coordinates (N×dim matrix)
        - "dofs_per_node::Int": Degrees of freedom per node

        # Notes
        - "macro_matrix" interpretation:
        - Dirichlet: Strain tensor (voigt or matrix form)
        - Neumann: Stress tensor
        - "node_pairs" should be matched using "match_opposite_periodic_nodes"
        - "coords" used to compute position-dependent constraints

        # Examples
        ```julia
        # 2% strain in xx-direction
        bc = HomogenizationBC(:dirichlet, [0.02 0; 0 0], nothing, [true, false], nodes, 2)

        # Periodic BC with 500 node pairs
        bc = HomogenizationBC(:periodic, G_macro, node_pairs, nothing, trues(3), nodes, 3)
    """
    struct HomogenizationBC
        type::Symbol                    # :dirichlet, :neumann, :periodic
        macro_matrix::Union{Matrix{Float64}, Vector{Float64}, Nothing} # Macro-scale matrix (Gmacro or Emacro)
        node_pairs::Union{Tuple{Vector{Int}, Vector{Int}}, Nothing}
        nodes_dof::Union{Vector{Int}, Nothing}
        dof_mask::Union{Vector{Bool}, Nothing}
        coords::Union{Matrix{Float64}, Nothing}       # N×dim matrix of nodal coordinates
        dofs_per_node::Int
    end

    function Base.show(io::IO, bc::HomogenizationBC)
        println(io, "HomogenizationBC(type = $(bc.type),")
        println(io, "  nodes_dof = $(bc.nodes_dof),")
        println(io, "  node_pairs = $(bc.node_pairs),")
        println(io, "  dofs_per_node = $(bc.dofs_per_node),")
        println(io, "  dof_mask = ", bc.dof_mask === nothing ? "nothing" : bc.dof_mask)
        println(io, "  macro_matrix = ", bc.macro_matrix === nothing ? "nothing" : "$((bc.macro_matrix))")
        println(io, "  coords = ", bc.coords === nothing ? "nothing" : "$(size(bc.coords))", ")")
    end

# =============================================
# Main Solver Function
# =============================================

    function solve!(
        K::Union{SparseMatrixCSC{Float64,Int}, Matrix{Float64}},
        f::Vector{Float64},
        u::Vector{Float64},
        bcs::Union{BoundaryCondition, HomogenizationBC, AbstractVector};
        problem_type::Symbol=:linear)
        # Normalize input to vector of BCs
        bcs = bcs isa Union{BoundaryCondition, HomogenizationBC} ? [bcs] : bcs
        
        # Validate all boundary conditions
        for bc in bcs
            bc isa Union{BoundaryCondition, HomogenizationBC} || 
            error("Invalid BC type: $(typeof(bc)). Must be BoundaryCondition or HomogenizationBC")
        end
        
        # Sort BCs: Dirichlet first, then Neumann, then Periodic
        sorted_bcs = sort(bcs, by=bc -> bc.type == :dirichlet ? 1 : bc.type == :neumann ? 2 : 3)
        
        # Initialize system matrices
        K_current = copy(K)
        F_current = copy(f)
        U_current = copy(u)
        
        periodic_applied = false
        homogenization_applied = false
        
        # Apply boundary conditions in sorted order
        for bc in sorted_bcs
            if bc isa BoundaryCondition
                
                 if bc.type == :dirichlet
                    K_current, F_current = apply_dirichlet_problem!(K_current, F_current, bc)
                elseif bc.type == :neumann
                    K_current = apply_neumann_problem!(K_current, bc)
                elseif bc.type == :periodic
                    K_current, F_current, U_current = apply_periodic_problem!(K_current, F_current, U_current, bc)
                    periodic_applied = true
                end
               
            else  # Standard BC
                homogenization_applied = true
                if bc.type == :dirichlet
                    F_current, U_current, _ = apply_dirichlet_homogenization!(bc, K_current, F_current, U_current)
                elseif bc.type == :neumann
                    F_current, U_current, _ = apply_neumann_homogenization!(bc, K_current, F_current, U_current)
                elseif bc.type == :periodic
                    K_current, F_current, U_current = apply_periodic_homogenization!(bc, K_current, F_current, U_current)
                    periodic_applied = true
                end
            end
        end
        
        # Solve system based on applied BC types
        if homogenization_applied && !periodic_applied && all(iszero, u)
            u = U_current  # Use precomputed solution
        elseif periodic_applied
            # Extract solution for original DOFs
            u = U_current[1:length(u)]
        else
            # Solve linear system
            u = K_current \ F_current
        end

        return u
    end

# =============================================
# Utility Functions
# =============================================

    function _get_global_dofs(nodes::AbstractVector{Int}, dofs_per_node::Int)
        # Get global DOFs for a set of nodes
        return vcat([(node-1)*dofs_per_node .+ (1:dofs_per_node) for node in nodes]...)
    end

    function _dimension(coords::Matrix)
        # Return spatial dimension (1, 2, or 3)
        return size(coords, 2)
    end

# =============================================
# Standard Boundary Condition Application
# =============================================

    function apply_dirichlet_problem!(K, f, bc::BoundaryCondition)
        # Get global DOFs and values
        global_dofs = _get_global_dofs(bc.nodes_dof, bc.dofs_per_node)
        values = repeat(bc.macro_matrix, length(bc.nodes_dof))
        
        # Apply DOF mask if provided
        if bc.dof_mask !== nothing
            full_mask = repeat(bc.dof_mask, length(bc.nodes_dof))
            global_dofs = global_dofs[full_mask]
            values = values[full_mask]
        end
      
        # Apply constraints
        for gdof in global_dofs
            # Zero out row and set diagonal to 1
            for j in axes(K, 2)
                K[gdof, j] = 0.0
            end
            K[gdof, gdof] = 1.0
            
            # Set RHS value
            f[gdof] = values[findfirst(==(gdof), global_dofs)]
        end
        
        return K, f
    end

    function apply_neumann_problem!(K, bc::BoundaryCondition)
        # Get global DOFs
        global_dofs = _get_global_dofs(bc.nodes_dof, bc.dofs_per_node)
        
        # Apply DOF mask if provided
        if bc.dof_mask !== nothing
            full_mask = repeat(bc.dof_mask, length(bc.nodes_dof))
            global_dofs = global_dofs[full_mask]
        end
        
        # Apply constraints
        for gdof in global_dofs
            # Zero out row and set diagonal to 1
            for j in axes(K, 2)
                K[gdof, j] = 0.0
            end
            K[gdof, gdof] = 1.0
        end
        
        return K
    end

    function apply_periodic_problem!(K, f, u, bc::BoundaryCondition)
        # Unpack master and slave nodes
        master_nodes, slave_nodes = bc.nodes_dof
        m = length(master_nodes) * bc.dofs_per_node
        
        # Build constraint matrix
        C = spzeros(m, size(K, 2))
        row = 1
        for (m_node, s_node) in zip(master_nodes, slave_nodes)
            for dof in 1:bc.dofs_per_node
                m_dof = (m_node-1)*bc.dofs_per_node + dof
                s_dof = (s_node-1)*bc.dofs_per_node + dof
                
                C[row, m_dof] = -1.0
                C[row, s_dof] = 1.0
                row += 1
            end
        end
        
        # Build augmented system
        K_aug = [K  C'
                C  spzeros(m, m)]
        f_aug = [f; zeros(m)]
        u_aug = [u; zeros(m)]
        
        return K_aug, f_aug, u_aug
    end

# =============================================
# Homogenization Boundary Condition Application
# =============================================

    function apply_dirichlet_homogenization!(bc::HomogenizationBC, K, F_glob, U_result)
        # Get active DOFs
        total_dofs = size(K, 1)
        active_dofs = collect(1:total_dofs)
        constrained_dofs = _get_global_dofs(bc.nodes_dof, bc.dofs_per_node)
         
        dim = _dimension(bc.coords)
        #  # Apply DOF mask if provided
        if bc.dof_mask !== nothing
            full_mask = repeat(bc.dof_mask, length(bc.nodes_dof))
          
            constrained_dofs = constrained_dofs[full_mask]
          
        end
        # Apply constraints
        for cdof in constrained_dofs
            # Compute constrained value from macro matrix
            node = div(cdof - 1, bc.dofs_per_node) + 1
            local_dof = (cdof - 1) % bc.dofs_per_node + 1
            val = dot(bc.macro_matrix[local_dof, 1:dim], bc.coords[node, 1:dim])
            
            # Apply constraint
            U_result[cdof] = val
            F_glob .-= K[:, cdof] .* val
        end
        
        # Solve for free DOFs
        free_dofs = setdiff(active_dofs, constrained_dofs)
        U_result[free_dofs] = K[free_dofs, free_dofs] \ F_glob[free_dofs]
        
        return F_glob, U_result, free_dofs
    end

    function apply_neumann_homogenization!(bc::HomogenizationBC, K, F_glob, U_result)
        # Get global DOFs
        global_dofs = _get_global_dofs(bc.nodes_dof, bc.dofs_per_node)
        dim = _dimension(bc.coords)
        macro_vals = zeros(length(global_dofs))
        
        # Compute macro values at nodes
        for (i, gdof) in enumerate(global_dofs)
            node = div(gdof - 1, bc.dofs_per_node) + 1
            local_dof = (gdof - 1) % bc.dofs_per_node + 1
            macro_vals[i] = dot(bc.macro_matrix[local_dof, 1:dim], bc.coords[node, 1:dim])
        end
        
        # Apply to RHS
        F_glob[global_dofs] .+= macro_vals
        
        return F_glob, U_result, collect(1:size(K, 1))
    end

    function apply_periodic_homogenization!(bc::HomogenizationBC, 
        K::Union{SparseMatrixCSC{Float64,Int}, Matrix{Float64}}, 
        F_glob::Vector{Float64}, 
        U_result::Vector{Float64})

        master_nodes, slave_nodes = bc.node_pairs
        dofs_per_node = bc.dofs_per_node
        Gmacro = bc.macro_matrix
        n_dofs = length(U_result)
        coords = bc.coords
        dim = size(coords, 2)
        Ne = size(coords, 1)

        # Identify corner nodes
        min_coords = minimum(coords, dims=1)
        max_coords = maximum(coords, dims=1)
        is_corner = [all(coords[i,d] ≈ min_coords[d] || coords[i,d] ≈ max_coords[d] for d in 1:dim) for i in 1:Ne]
        corners = Set(findall(is_corner))

        # Classify node pairs
        regular_pairs = Tuple{Int,Int}[]
        corner_pairs  = Tuple{Int,Int}[]
        for (m, s) in zip(master_nodes, slave_nodes)
            ((m in corners) || (s in corners)) ? push!(corner_pairs, (m, s)) : push!(regular_pairs, (m, s))
        end

        dof_mask = bc.dof_mask === nothing ? trues(dofs_per_node) : bc.dof_mask
        active_dofs = sum(dof_mask)

        function build_constraints(pairs::Vector{Tuple{Int,Int}})
            num_rows = length(pairs) * active_dofs
            C = spzeros(num_rows, n_dofs)
            R = zeros(num_rows)

            row = 1
            for (m_node, s_node) in pairs
                Δ = coords[m_node, :] .- coords[s_node, :]
                for d in 1:dofs_per_node
                    if !dof_mask[d]; continue; end
                    mgdof = (m_node-1)*dofs_per_node + d
                    sgdof = (s_node-1)*dofs_per_node + d
                    C[row, sgdof] = -1.0
                    C[row, mgdof] = 1.0
                    R[row] = dot(Gmacro[d, 1:dim], Δ)
                    row += 1
                end
            end
            return C[1:row-1, :], R[1:row-1]
        end

        # Regular constraints
        C_reg, R_reg = build_constraints(regular_pairs)

        # Dynamically generate corner constraints based on mesh geometry
        function generate_corner_constraints()
            if dim == 2
                corner_coords = [(x,y) for x in (min_coords[1], max_coords[1]), y in (min_coords[2], max_coords[2])]
                idx_map = Dict{Tuple{Float64,Float64}, Int}()
                for cc in corner_coords
                    idx = findmin([norm(coords[i,:] .- collect(cc)) for i in corners])[2]
                    idx_map[cc] = collect(corners)[idx]
                end
                return [
                    (idx_map[(min_coords[1], min_coords[2])], idx_map[(max_coords[1], max_coords[2])]),
                    (idx_map[(min_coords[1], max_coords[2])], idx_map[(max_coords[1], min_coords[2])]),
                    (idx_map[(min_coords[1], min_coords[2])], idx_map[(max_coords[1], min_coords[2])])
                ]
            elseif dim == 3
                return [] # Extend for 3D if needed
            else
                return []
            end
        end

        corner_constraints = generate_corner_constraints()
        C_corner, R_corner = build_constraints(corner_constraints)

        # Combine constraints
        C = [C_reg; C_corner]
        R = [R_reg; R_corner]

        # Penalty augmentation
        α = maximum(abs.(K))
        K_aug = [K α*C'; α*C spzeros(size(C,1), size(C,1))]
        F_aug = [F_glob; α*R]
        U_aug = K_aug \ F_aug

        return K_aug, F_aug, U_aug
    end
    
    """
        match_opposite_periodic_nodes(coords, tol=1e-8)

        Match nodes on opposite faces for periodic boundary conditions.

        # Arguments
        - "coords": N×dim matrix of nodal coordinates
        - "tol": Geometric tolerance for matching

        # Returns
        Vector of (master, slave) node pairs

        # Algorithm
        1. For each dimension:
        - Find nodes on min/max faces
        - Match nodes by closest distance in perpendicular directions
        2. Returns all matched pairs

        # Notes
        - Essential for homogenization with periodic BCs
        - Handles 2D and 3D geometries
        - Ensures compatible mesh for periodic constraints
    """
    function match_opposite_periodic_nodes(coords::Matrix{Float64}, tol::Float64=1e-8)
        dim = size(coords, 2)
        pairs = Tuple{Int, Int}[]

        for d in 1:dim
            min_val = minimum(coords[:, d])
            max_val = maximum(coords[:, d])

            side_min = findall(x -> abs(x[d] - min_val) < tol, eachrow(coords))
            side_max = findall(x -> abs(x[d] - max_val) < tol, eachrow(coords))

            for i in side_min
                ci = coords[i, setdiff(1:dim, [d])]
                closest_j = argmin([norm(ci - coords[j, setdiff(1:dim, [d])]) for j in side_max])
                j = side_max[closest_j]
                push!(pairs, (i, j))
            end
        end

        return pairs
    end
