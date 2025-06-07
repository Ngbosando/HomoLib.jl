

###############################################################################
# Structures de données
###############################################################################
struct BoundaryCondition
    type::Symbol           # :dirichlet, :neumann, :periodic
    nodes_dof::Union{Vector{Int}, Tuple{Vector{Int}, Vector{Int}}}
    dof_mask::Union{Vector{Bool}, Nothing}        # Vrai ou faux selon qu'on impose le dof
    macro_matrix::Union{Vector{Float64}, Nothing} # Valeurs imposées (Dirichlet/Neumann)
    dofs_per_node::Int
end
function Base.show(io::IO, bc::BoundaryCondition)
    println(io, "BoundaryCondition(type = $(bc.type),")
    println(io, "  nodes_dof = $(bc.nodes_dof),")
    println(io, "  dofs_per_node = $(bc.dofs_per_node),")
    println(io, "  dof_mask = ", bc.dof_mask === nothing ? "nothing" : bc.dof_mask)
    println(io, "  macro_matrix = ", bc.macro_matrix === nothing ? "nothing" : bc.macro_matrix, ")")
end

struct HomogenizationBC
    type::Symbol                    # :dirichlet, :neumann, :periodic
    macro_matrix::Union{Matrix{Float64}, Nothing} # Gmacro ou Emacro
    node_pairs::Union{Tuple{Vector{Int}, Vector{Int}}, Nothing}
    nodes_dof::Union{Vector{Int}, Nothing}
    dof_mask::Union{Vector{Bool}, Nothing}
    coords::Union{Matrix{Float64}, Nothing}       # Matrice N×dim contenant x,y,(z) par nœud
    dofs_per_node::Int
end
function Base.show(io::IO, bc::HomogenizationBC)
    println(io, "HomogenizationBC(type = $(bc.type),")
    println(io, "  nodes_dof = $(bc.nodes_dof),")
    println(io, "  node_pairs = $(bc.node_pairs),")
    println(io, "  dofs_per_node = $(bc.dofs_per_node),")
    println(io, "  dof_mask = ", bc.dof_mask === nothing ? "nothing" : bc.dof_mask)
    println(io, "  macro_matrix = ", bc.macro_matrix === nothing ? "nothing" : bc.macro_matrix)
    println(io, "  coords = ", bc.coords === nothing ? "nothing" : "$(size(bc.coords))", ")")
end


###############################################################################
# Fonction de résolution principale
###############################################################################
function solve!(K::Union{SparseMatrixCSC{Float64,Int},Matrix{Float64}},
    f::Vector{Float64},
    u::Vector{Float64},
    bcs::Vector{Any}; problem_type=:linear)
    # Valider tous les éléments
    for bc in bcs
        bc isa Union{BoundaryCondition, HomogenizationBC} || 
        error("All boundary conditions must be BoundaryCondition or HomogenizationBC, got: $(typeof(bc))")
    end
    # Convertir vers bon type statique
    typed_bcs = Union{BoundaryCondition, HomogenizationBC}[bc for bc in bcs]
    return solve!(K, f, u, typed_bcs; problem_type)
end

function solve!(
    K::Union{SparseMatrixCSC{Float64,Int}, Matrix{Float64}},
    f::Vector{Float64},
    u::Vector{Float64},
    bcs::Union{BoundaryCondition, HomogenizationBC, Vector{<:Union{BoundaryCondition, HomogenizationBC}}};
    problem_type=:linear)

    bcs = bcs isa Union{BoundaryCondition, HomogenizationBC} ? [bcs] : bcs

    sorted_bcs = sort(bcs, by = bc -> if bc.type == :dirichlet 1 elseif bc.type == :neumann 2 else 3 end)

    K_current = copy(K)
    F_current = copy(f)
    U_current = copy(u)

    periodic_applied = false
    homogenization_applied = false
   
        for bc in sorted_bcs
            if bc isa HomogenizationBC
                if bc.type == :dirichlet
                    F_current, U_current, _ = apply_dirichlet_homogenization!(bc, K_current, F_current, U_current)
                    homogenization_applied = true
                elseif bc.type == :neumann
                    F_current, U_current, _ = apply_neumann_homogenization!(bc, K_current, F_current, U_current)
                    homogenization_applied = true
                elseif bc.type == :periodic
                    K_current, F_current, U_current = apply_periodic_homogenization!(bc, K_current, F_current, U_current)
                    periodic_applied = true
                    homogenization_applied = true
             
                end
            else
                if bc.type == :dirichlet
                    K_current, F_current = apply_dirichlet_problem!(K_current, F_current, bc)
                    
                elseif bc.type == :neumann
                    apply_neumann_problem!(F_current, bc)
                elseif bc.type == :periodic
                    K_current, F_current, U_current = apply_periodic_problem!(K_current, F_current, U_current, bc)
                    periodic_applied = true
                end
            end
        end
    
        
    if homogenization_applied && !periodic_applied && all(iszero, u)
        u = U_current
        elseif periodic_applied
        # If periodic constraints were applied, extract the solution for the original DOFs
        ndofs = length(u)
        u = U_current[1:ndofs]
    else
        u_result = K_current \ F_current
        u = u_result
        
    end

    return u
end



###############################################################################
# Fonctions utilitaires
###############################################################################
function _get_global_dofs(nodes::AbstractVector{Int}, dofs_per_node::Int)
    # Récupère la liste des dofs globaux d’un ensemble de nœuds
    return vcat([(node-1)*dofs_per_node .+ (1:dofs_per_node) for node in nodes]...)
end

function _dimension(coords)
    # Retourne 1, 2 ou 3 si coords est N×1, N×2, N×3, etc.
    return size(coords, 2)
end

###############################################################################
# BC standard : Dirichlet, Neumann, Périodique
###############################################################################
function apply_dirichlet_problem!(K, f, bc::BoundaryCondition)
    global_dofs = _get_global_dofs(bc.nodes_dof, bc.dofs_per_node)
    values = repeat(bc.macro_matrix, length(bc.nodes_dof))
    N = size(K, 2)
 
    for gdof in global_dofs
      
        for j in 1: N
            K[gdof, j] = 0.0
        end
      

        K[gdof, gdof] = 1.0
       
    end
    
    f[global_dofs] .= values
   
    return K, f
end


function apply_neumann_problem!(f, bc::BoundaryCondition)
    # Ajout direct de la valeur imposée à f[gdof]
    global_dofs = _get_global_dofs(bc.nodes_dof, bc.dofs_per_node)
    imposed_vals = repeat(bc.macro_matrix, length(bc.nodes_dof))

    for (i, gdof) in enumerate(global_dofs)
        f[gdof] += imposed_vals[i]
    end
    return f
end

function apply_periodic_problem!(K, f, u, bc::BoundaryCondition)
    # Suppose bc.nodes_dof = (master_nodes, slave_nodes)
    master_nodes_dof, slave_nodes_dof = bc.nodes_dof
    m = length(master_nodes_dof) * bc.dofs_per_node

    # Matrice de contraintes
    C = spzeros(m, size(K, 2))
    row = 1
    for (m_node, s_node) in zip(master_nodes_dof, slave_nodes_dof)
        for dof in 1:bc.dofs_per_node
            C[row, (s_node-1)*bc.dofs_per_node + dof] =  1.0
            C[row, (m_node-1)*bc.dofs_per_node + dof] = -1.0
            row += 1
        end
    end

    # Système augmenté
    K_aug = [K  C'
             C  spzeros(m,m)]
    f_aug = [f; zeros(m)]
    u_aug = [u; zeros(m)]

    return K_aug, f_aug, u_aug
end

###############################################################################
# BC « Homogenization »
# (cas où on impose un champ macro, typiquement dans la résolution RVE)
###############################################################################
function apply_dirichlet_homogenization!(bc::HomogenizationBC, K, F_glob, U_result)
    total_dofs = size(K, 1)
    ListDOF = collect(1:total_dofs)
    constrained_dofs = _get_global_dofs(bc.nodes_dof, bc.dofs_per_node)
   

    # On ne contraint que ceux encore "actifs"
    dofs_to_constrain = intersect(constrained_dofs, ListDOF)
  
    dim = _dimension(bc.coords)
 
    for cdof in dofs_to_constrain
        node = div(cdof - 1, bc.dofs_per_node) + 1
        local_dof = (cdof - 1) % bc.dofs_per_node + 1
       
        # Produit scalaire macro_matrix[local_dof, :] avec coords[node, :]
        # On suppose que bc.macro_matrix est dimensionné en (dofs_per_node, dim)
        val = dot(bc.macro_matrix[local_dof, 1:dim], bc.coords[node, 1:dim])

        # Imposer la valeur
        U_result[cdof] = val
        # Mettre à jour le second membre
        for j in ListDOF
            if j != cdof
                F_glob[j] -= K[j, cdof] * val
            end
        end
    end

    # Mise à jour de la matrice K pour supprimer les cdof
    new_ListDOF = setdiff(ListDOF, dofs_to_constrain)
    K_free = K[new_ListDOF, new_ListDOF]
    U_result[new_ListDOF] = K_free \ F_glob[new_ListDOF]

    return F_glob, U_result, new_ListDOF
end

function apply_neumann_homogenization!(bc::HomogenizationBC, K, F_glob, U_result)
    total_dofs = size(K, 1)
    ListDOF = collect(1:total_dofs)
    global_dofs = _get_global_dofs(bc.nodes_dof, bc.dofs_per_node)

    dim = _dimension(bc.coords)
    macro_matrix_loc = zeros(length(global_dofs))

    # On calcule la « traction » en fonction de la macro_matrix et des coords
    # Ex : bc.macro_matrix[d, :] ~ un vecteur de dimension "dim"
    for (i, gdof) in enumerate(global_dofs)
        node = div(gdof - 1, bc.dofs_per_node) + 1
        local_dof = (gdof - 1) % bc.dofs_per_node + 1
        macro_matrix_loc[i] = dot(bc.macro_matrix[local_dof, 1:dim], bc.coords[node, 1:dim])
    end

    # Ajout au second membre
    F_glob[global_dofs] .+= macro_matrix_loc

    return F_glob, U_result, ListDOF
end

# function apply_periodic_homogenization!(
#     bc::HomogenizationBC, 
#     K::Union{SparseMatrixCSC{Float64,Int}, Matrix{Float64}}, 
#     F_glob::Vector{Float64}, 
#     U_result::Vector{Float64})
#     # Récupération
#     master_nodes, slave_nodes = bc.node_pairs
#     dofs_per_node = bc.dofs_per_node
#     Gmacro = bc.macro_matrix
#     n_dofs = length(U_result)
#     coords = bc.coords
#     dim = size(coords, 2)  # 1D, 2D, ou 3D
#     Ne = size(coords, 1)   # nombre de nœuds

#     # -------------------------------------------------------------------------
#     # Détection des nœuds corners via min_coords / max_coords
#     # -------------------------------------------------------------------------
#     min_coords = minimum(coords, dims=1)  # vector(1×dim)
#     max_coords = maximum(coords, dims=1)  # idem

#     # Pour chaque nœud i, on vérifie s'il est min ou max dans TOUTES les directions
#     is_corner = [all(coords[i,d] ≈ min_coords[d] || coords[i,d] ≈ max_coords[d] 
#                      for d in 1:dim) 
#                  for i in 1:Ne]
#     corners = findall(is_corner)   # indices des nœuds "coin"

#     # -------------------------------------------------------------------------
#     # On sépare master/slave pairs en corner_pairs et regular_pairs
#     # -------------------------------------------------------------------------
#     regular_pairs = Tuple{Int,Int}[]
#     corner_pairs  = Tuple{Int,Int}[]

#     for (m, s) in zip(master_nodes, slave_nodes)
#         # si un des deux est un corner, on le met dans corner_pairs
#         if (m in corners) || (s in corners)
#             push!(corner_pairs, (m, s))
#         else
#             push!(regular_pairs, (m, s))
#         end
#     end

#     # -------------------------------------------------------------------------
#     # On gère un dof_mask pour n'imposer que certains DDL (ex. dim=2 => d=1,2)
#     # -------------------------------------------------------------------------
#     dof_mask = bc.dof_mask === nothing ? trues(dofs_per_node) : bc.dof_mask
#     active_dofs = sum(dof_mask)

#     # -------------------------------------------------------------------------
#     # Contrainte sur paires "regular"
#     # -------------------------------------------------------------------------
#     num_reg = length(regular_pairs)*active_dofs
#     C_reg = spzeros(num_reg, n_dofs)
#     R_reg = zeros(num_reg)

#     row = 1
#     for (m_node, s_node) in regular_pairs
#         Δ = coords[m_node, :] .- coords[s_node, :]

#         for d in 1:dofs_per_node
#             if !dof_mask[d]; continue; end
#             mgdof = (m_node-1)*dofs_per_node + d
#             sgdof = (s_node-1)*dofs_per_node + d

#             C_reg[row, sgdof] = -1.0
#             C_reg[row, mgdof] =  1.0
#             R_reg[row] = dot(Gmacro[d, 1:dim], Δ)  # impose Gmacro * Δ

#             row += 1
#         end
#     end
#     C_reg = C_reg[1:row-1, :]
#     R_reg = R_reg[1:row-1]

#     # -------------------------------------------------------------------------
#     # Contrainte sur paires "corner"
#     # => On utilise le même "pattern" que le code d'origine :
#     #    corner_constraints = if dim==2 ...
#     # -------------------------------------------------------------------------
#     # Dans le code initial, on construisait un tableau corner_constraints
#     #    si dim == 2 => [ (4,3), (1,4), (2,3) ]
#     #    si dim == 3 => un exemple "placeholder"
#     # etc.
#     # Ensuite on applique un 2e "C_corner" distinct.
#     # -------------------------------------------------------------------------
#     corner_constraints = if dim == 2
#         # Sélection EXACTE du code d'origine
#         [(4,3), (1,4), (2,3)]
#     elseif dim == 3
#         # Placeholder : tu mettras ici les paires corner
#         # Ex: [(8,7), (1,8), (5,8)] ou tout autre
#         # en fonction de ton maillage
#         []
#     else
#         # Pour dim==1, en général on a 2 "coins", à adapter éventuellement
#         []
#     end

#     # Construit la matrice de contraintes "corner" en se basant sur corner_constraints
#     num_corner = length(corner_constraints)*active_dofs
#     C_corner = spzeros(num_corner, n_dofs)
#     R_corner = zeros(num_corner)

#     c_row = 1
#     for (m_node, s_node) in corner_constraints
#         Δ = coords[m_node, :] .- coords[s_node, :]

#         for d in 1:dofs_per_node
#             if !dof_mask[d]; continue; end
#             mgdof = (m_node-1)*dofs_per_node + d
#             sgdof = (s_node-1)*dofs_per_node + d

#             C_corner[c_row, sgdof] = -1.0
#             C_corner[c_row, mgdof] =  1.0
#             R_corner[c_row] = dot(Gmacro[d, 1:dim], Δ)

#             c_row += 1
#         end
#     end
#     C_corner = C_corner[1:c_row-1, :]
#     R_corner = R_corner[1:c_row-1]

#     # -------------------------------------------------------------------------
#     # Assemblage total des contraintes
#     # -------------------------------------------------------------------------
#     C = [C_reg; C_corner]
#     R = [R_reg; R_corner]

#     # Système augmenté (pénalisation)
#     α = maximum(abs.(K))  # facteur de pénalisation
#     K_aug = [
#         K     α*C'
#         α*C   spzeros(size(C,1), size(C,1))
#     ]
#     F_aug = [F_glob; α*R]

#     # On résout
#     U_aug = K_aug \ F_aug

#     return K_aug, F_aug, U_aug
# end
# function apply_periodic_homogenization!(bc::HomogenizationBC, 
#     K::Union{SparseMatrixCSC{Float64,Int}, Matrix{Float64}}, 
#     F_glob::Vector{Float64}, 
#     U_result::Vector{Float64})

#     master_nodes, slave_nodes = bc.node_pairs
#     dofs_per_node = bc.dofs_per_node
#     Gmacro = bc.macro_matrix
#     n_dofs = length(U_result)
#     coords = bc.coords
#     dim = size(coords, 2)
#     Ne = size(coords, 1)

#     # Identify corner nodes
#     min_coords = minimum(coords, dims=1)
#     max_coords = maximum(coords, dims=1)
#     is_corner = [all(coords[i,d] ≈ min_coords[d] || coords[i,d] ≈ max_coords[d] for d in 1:dim) for i in 1:Ne]
#     corners = Set(findall(is_corner))

#     # Classify node pairs
#     regular_pairs = Tuple{Int,Int}[]
#     corner_pairs  = Tuple{Int,Int}[]
#     for (m, s) in zip(master_nodes, slave_nodes)
#         ((m in corners) || (s in corners)) ? push!(corner_pairs, (m, s)) : push!(regular_pairs, (m, s))
#     end

#     dof_mask = bc.dof_mask === nothing ? trues(dofs_per_node) : bc.dof_mask
#     active_dofs = sum(dof_mask)

#     function build_constraints(pairs::Vector{Tuple{Int,Int}})
#         num_rows = length(pairs) * active_dofs
#         C = spzeros(num_rows, n_dofs)
#         R = zeros(num_rows)

#         row = 1
#         for (m_node, s_node) in pairs
#             Δ = coords[m_node, :] .- coords[s_node, :]
#             for d in 1:dofs_per_node
#                 if !dof_mask[d]; continue; end
#                 mgdof = (m_node-1)*dofs_per_node + d
#                 sgdof = (s_node-1)*dofs_per_node + d
#                 C[row, sgdof] = -1.0
#                 C[row, mgdof] = 1.0
#                 R[row] = dot(Gmacro[d, 1:dim], Δ)
#                 row += 1
#             end
#         end
#         return C[1:row-1, :], R[1:row-1]
#     end

#     # Regular constraints
#     C_reg, R_reg = build_constraints(regular_pairs)

#     # Dynamically generate corner constraints based on mesh geometry
#     function generate_corner_constraints()
#         if dim == 2
#             corner_coords = [(x,y) for x in (min_coords[1], max_coords[1]), y in (min_coords[2], max_coords[2])]
#             idx_map = Dict{Tuple{Float64,Float64}, Int}()
#             for cc in corner_coords
#                 idx = findmin([norm(coords[i,:] .- collect(cc)) for i in corners])[2]
#                 idx_map[cc] = collect(corners)[idx]
#             end
#             return [
#                 (idx_map[(min_coords[1], min_coords[2])], idx_map[(max_coords[1], max_coords[2])]),
#                 (idx_map[(min_coords[1], max_coords[2])], idx_map[(max_coords[1], min_coords[2])]),
#                 (idx_map[(min_coords[1], min_coords[2])], idx_map[(max_coords[1], min_coords[2])])
#             ]
#         elseif dim == 3
#             return [] # Extend for 3D if needed
#         else
#             return []
#         end
#     end

#     corner_constraints = generate_corner_constraints()
#     C_corner, R_corner = build_constraints(corner_constraints)

#     # Combine constraints
#     C = [C_reg; C_corner]
#     R = [R_reg; R_corner]

#     # Penalty augmentation
#     α = maximum(abs.(K))
#     K_aug = [K α*C'; α*C spzeros(size(C,1), size(C,1))]
#     F_aug = [F_glob; α*R]
#     U_aug = K_aug \ F_aug

#     return K_aug, F_aug, U_aug
# end

# function match_opposite_periodic_nodes(coords::Matrix{Float64}, tol::Float64=1e-8)
#     dim = size(coords, 2)
#     pairs = Tuple{Int, Int}[]

#     for d in 1:dim
#         min_val = minimum(coords[:, d])
#         max_val = maximum(coords[:, d])

#         side_min = findall(x -> abs(x[d] - min_val) < tol, eachrow(coords))
#         side_max = findall(x -> abs(x[d] - max_val) < tol, eachrow(coords))

#         for i in side_min
#             ci = coords[i, setdiff(1:dim, [d])]
#             closest_j = argmin([norm(ci - coords[j, setdiff(1:dim, [d])]) for j in side_max])
#             j = side_max[closest_j]
#             push!(pairs, (i, j))
#         end
#     end

#     return pairs
# end
function apply_periodic_homogenization!(bc::HomogenizationBC, K, F_glob, U_result)

    # Retrieve problem parameters
    master_nodes, slave_nodes = bc.node_pairs
    dofs_per_node = bc.dofs_per_node
    Gmacro = bc.macro_matrix
    n_dofs = length(U_result)
    coords = bc.coords
    dim = size(coords, 2)
    n_nodes = size(coords, 1)
    tol = 1e-10  # Tolerance for floating point comparisons

    # Find min/max coordinates for corner detection
    min_coords = dropdims(minimum(coords, dims=1), dims=1)
    max_coords = dropdims(maximum(coords, dims=1), dims=1)

    # Detect corner nodes (min or max in all dimensions)
    is_corner = fill(false, n_nodes)
    for i in 1:n_nodes
        is_corner[i] = all(
            (abs(coords[i,d] - min_coords[d]) < tol) || 
            (abs(coords[i,d] - max_coords[d]) < tol) 
            for d in 1:dim
        )
    end
    corners = findall(is_corner)

    # Separate pairs into regular and corner
    regular_pairs = Tuple{Int,Int}[]
    corner_pairs  = Tuple{Int,Int}[]
    for (m, s) in zip(master_nodes, slave_nodes)
        if m in corners || s in corners
            push!(corner_pairs, (m, s))
        else
            push!(regular_pairs, (m, s))
        end
    end

    # Determine active DOFs
    dof_mask = bc.dof_mask === nothing ? trues(dofs_per_node) : bc.dof_mask
    active_dofs = sum(dof_mask)

    # ========================================================================
    # Regular constraints (all pairs)
    # ========================================================================
    n_reg_constraints = length(regular_pairs) * active_dofs
    C_reg = spzeros(n_reg_constraints, n_dofs)
    R_reg = zeros(n_reg_constraints)
    row = 1

    for (m_node, s_node) in regular_pairs
        Δ = coords[m_node, :] .- coords[s_node, :]
        
        for d in 1:dofs_per_node
            dof_mask[d] || continue  # Skip inactive DOFs
            
            mgdof = (m_node-1)*dofs_per_node + d
            sgdof = (s_node-1)*dofs_per_node + d

            C_reg[row, mgdof] =  1.0
            C_reg[row, sgdof] = -1.0
            R_reg[row] = dot(Gmacro[d, 1:dim], Δ)
            
            row += 1
        end
    end
    C_reg = C_reg[1:row-1, :]
    R_reg = R_reg[1:row-1]

    # ========================================================================
    # Corner constraints (spanning tree approach)
    # ========================================================================
    if !isempty(corner_pairs)
        # Find reference corner (min_coords in all dimensions)
        ref_corner = findfirst(i -> all(abs.(coords[i, :] .- min_coords) .< tol), 1:n_nodes)
        ref_corner === nothing && error("Reference corner node not found")

        # Build graph of corner nodes
        graph = Dict{Int, Vector{Int}}()
        for node in corners
            graph[node] = Int[]
        end

        # Add edges from corner pairs (undirected)
        for (m, s) in corner_pairs
            if m in keys(graph) && s in keys(graph)
                push!(graph[m], s)
                push!(graph[s], m)
            end
        end

        # BFS to find spanning tree edges
        spanning_tree_edges = Tuple{Int,Int}[]
        visited = Set{Int}([ref_corner])
        queue = [ref_corner]

        while !isempty(queue)
            u = popfirst!(queue)
            for v in graph[u]
                if !(v in visited)
                    push!(spanning_tree_edges, (u, v))
                    push!(visited, v)
                    push!(queue, v)
                end
            end
        end

        # Apply constraints for spanning tree edges
        n_corner_constraints = length(spanning_tree_edges) * active_dofs
        C_corner = spzeros(n_corner_constraints, n_dofs)
        R_corner = zeros(n_corner_constraints)
        c_row = 1

        for (u, v) in spanning_tree_edges
            Δ = coords[u, :] .- coords[v, :]
            for d in 1:dofs_per_node
                dof_mask[d] || continue
                
                u_dof = (u-1)*dofs_per_node + d
                v_dof = (v-1)*dofs_per_node + d
                
                C_corner[c_row, u_dof] =  1.0
                C_corner[c_row, v_dof] = -1.0
                R_corner[c_row] = dot(Gmacro[d, 1:dim], Δ)
                
                c_row += 1
            end
        end
        C_corner = C_corner[1:c_row-1, :]
        R_corner = R_corner[1:c_row-1]
    else
        C_corner = spzeros(0, n_dofs)
        R_corner = zeros(0)
    end

    # ========================================================================
    # Combine constraints and solve
    # ========================================================================
    C = [C_reg; C_corner]
    R = [R_reg; R_corner]

    # Penalty method
    α = maximum(abs, K) * 1e6  # Penalty factor
    K_aug = [
        K           α * C';
        α * C   spzeros(size(C, 1), size(C, 1))
    ]
    F_aug = [F_glob; α * R]
    U_aug = K_aug \ F_aug

    return K_aug, F_aug, U_aug
end