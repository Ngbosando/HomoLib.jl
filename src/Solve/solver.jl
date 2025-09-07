#  ABSTRACT TYPE DEFINITIONS

abstract type AbstractBoundaryCondition end
abstract type StandardBC <: AbstractBoundaryCondition end
abstract type HomogenizationBC{BC_TYPE, M} <: AbstractBoundaryCondition end

# Standard Boundary Conditions
struct DirichletBC{T<:Real, VI<:AbstractVector{Int}, VV<:AbstractVector{T}, VB<:AbstractVector{Bool}} <: StandardBC
    nodes::VI
    dof_mask::VB
    values::VV
    dofs_per_node::Int
end

struct PeriodicBC{VI<:AbstractVector{Int}} <: StandardBC
    master_nodes::VI
    slave_nodes::VI
    dofs_per_node::Int
end

# Homogenization Boundary Conditions
struct DirichletHomogenizationBC{M<:AbstractMatrix{<:Real}, VI<:AbstractVector{Int}, VB<:AbstractVector{Bool}, MC<:AbstractMatrix{<:Real}} <: HomogenizationBC{:dirichlet, M}
    macro_matrix::M
    boundary_nodes::VI
    coords::MC
    dofs_per_node::Int
    dof_mask::VB
end

struct PeriodicHomogenizationBC{M<:AbstractMatrix{<:Real}, VI<:AbstractVector{Int}, MC<:AbstractMatrix{<:Real}, VB<:AbstractVector{Bool}} <: HomogenizationBC{:periodic, M}
    macro_matrix::M
    node_pairs::Tuple{VI, VI}
    coords::MC
    dofs_per_node::Int
    dof_mask::VB
end

#  UTILITY FUNCTIONS 

function get_global_dofs(nodes::AbstractVector{Int}, dofs_per_node::Int)
    num_nodes = length(nodes)
    global_dofs = Vector{Int}(undef, num_nodes * dofs_per_node)
    @inbounds for (i, node) in enumerate(nodes)
        offset = (i - 1) * dofs_per_node
        base_dof = (node - 1) * dofs_per_node
        for j in 1:dofs_per_node
            global_dofs[offset + j] = base_dof + j
        end
    end
    return global_dofs
end

get_dimension(coords::Matrix) = size(coords, 2)

#  BOUNDARY CONDITION APPLICATION 

function apply_bc_modifier!(K_in::SparseMatrixCSC, F_in::Vector{Float64}, U_in::Vector{Float64}, bc::DirichletBC)
    K_out, F_out = apply_dirichlet_problem!(K_in, F_in, bc)
    return K_out, F_out, U_in
end

function apply_bc_modifier!(K_in::SparseMatrixCSC, F_in::Vector{Float64}, U_in::Vector{Float64}, bc::PeriodicBC)
    K_aug, F_aug, U_aug = apply_periodic_problem!(K_in, F_in, U_in, bc)
    return K_aug, F_aug, U_aug
end

function apply_bc_modifier!(K_in::SparseMatrixCSC, F_in::Vector{Float64}, U_in::Vector{Float64}, bc::DirichletHomogenizationBC)
    F_out, U_out, free_dofs = apply_dirichlet_homogenization!(bc, K_in, F_in, U_in)
    return F_out, U_out, free_dofs
end

function apply_bc_modifier!(K_in::SparseMatrixCSC, F_in::Vector{Float64}, U_in::Vector{Float64}, bc::PeriodicHomogenizationBC)
    K_aug, F_aug, U_aug = apply_periodic_homogenization!(bc, K_in, F_in, U_in)
    return K_aug, F_aug, U_aug
end

# Periodic constraint application
function apply_periodic_problem!(K::SparseMatrixCSC, f::Vector{Float64}, u::Vector{Float64}, bc::PeriodicBC)
    master, slave = bc.master_nodes, bc.slave_nodes
    n_pairs = length(master)
    m = n_pairs * bc.dofs_per_node
    n_dofs = size(K, 2)
    
    # Preallocate COO arrays
    I_C = Vector{Int}(undef, 2*m)
    J_C = Vector{Int}(undef, 2*m)
    V_C = Vector{Float64}(undef, 2*m)
    
    idx = 1
    @inbounds for i in 1:n_pairs
        m_node = master[i]
        s_node = slave[i]
        base_m = (m_node - 1) * bc.dofs_per_node
        base_s = (s_node - 1) * bc.dofs_per_node
        
        for j in 1:bc.dofs_per_node
            # Constraint: u_m - u_s = 0
            row = (i-1)*bc.dofs_per_node + j
            I_C[idx] = row; J_C[idx] = base_m + j; V_C[idx] = 1.0; idx += 1
            I_C[idx] = row; J_C[idx] = base_s + j; V_C[idx] = -1.0; idx += 1
        end
    end
    
    C = sparse(view(I_C,1:idx-1), view(J_C,1:idx-1), view(V_C,1:idx-1), m, n_dofs)
    K_aug = [K C'; C spzeros(m, m)]
    return K_aug, [f; zeros(m)], [u; zeros(m)]
end

function apply_dirichlet_problem!(K::SparseMatrixCSC, f::Vector{Float64}, bc::DirichletBC)
    # Get global DOFs and values
    global_dofs = get_global_dofs(bc.nodes, bc.dofs_per_node)
    values = repeat(bc.values, length(bc.nodes))
    
    # Apply mask if it exists
    if !isempty(bc.dof_mask)
        full_mask = repeat(bc.dof_mask, length(bc.nodes))
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

function apply_dirichlet_homogenization!(bc::DirichletHomogenizationBC, K::SparseMatrixCSC, F_glob::Vector{Float64}, U_result::Vector{Float64})
    nodes = bc.boundary_nodes
    n_nodes = length(nodes)
    dofs_per_node = bc.dofs_per_node
    dim = get_dimension(bc.coords)
    
    # Precompute constrained DOFs
    constrained_dofs = get_global_dofs(nodes, dofs_per_node)
    if !isempty(bc.dof_mask)
        mask = repeat(bc.dof_mask, n_nodes)
        constrained_dofs = constrained_dofs[mask]
    end
    
    # In-place updates
    F_out = copy(F_glob)
    U_out = copy(U_result)
    
    @inbounds for cdof in constrained_dofs
        node = div(cdof - 1, dofs_per_node) + 1
        local_dof = (cdof - 1) % dofs_per_node + 1
        val = dot(@view(bc.macro_matrix[local_dof, 1:dim]), @view(bc.coords[node, 1:dim]))
        
        U_out[cdof] = val
        for i in eachindex(F_out)
            F_out[i] -= K[i, cdof] * val
        end
    end
    
    # Return free DOFs as range for efficiency
    free_dofs = setdiff(1:length(F_glob), constrained_dofs)
    return F_out, U_out, free_dofs
end

function apply_periodic_homogenization!(bc::PeriodicHomogenizationBC, 
                                      K::SparseMatrixCSC{Float64,Int}, 
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

        return K_aug, F_aug, U_result
end

#  MAIN SOLVER FUNCTION 

function solve!(
    K::SparseMatrixCSC{Float64,Int},
    f::Vector{Float64},
    u::Vector{Float64},
    bcs::AbstractVector{<:AbstractBoundaryCondition})
    
    # Process boundary conditions
    K_current = copy(K)
    F_current = copy(f)
    U_current = copy(u)
    free_dofs = Colon()  # Default to all DOFs
    condensed = false
    
    @inbounds for bc in bcs
        if bc isa DirichletHomogenizationBC
            F_current, U_current, free_dofs = apply_dirichlet_homogenization!(bc, K_current, F_current, U_current)
            condensed = true
        else
            K_new, F_new, U_new = apply_bc_modifier!(K_current, F_current, U_current, bc)
            if size(K_new) != size(K_current)
                # System was augmented
                K_current, F_current, U_current = K_new, F_new, U_new
                condensed = false
            else
                K_current, F_current, U_current = K_new, F_new, U_new
            end
        end
    end

    # Solve system
    if condensed
        # Solve condensed system
        K_ff = K_current[free_dofs, free_dofs]
        F_f = F_current[free_dofs]
        U_current[free_dofs] = K_ff \ F_f
        return U_current
    else
        U_current = K_current \ F_current
        return U_current[1:length(u)]
    end
end