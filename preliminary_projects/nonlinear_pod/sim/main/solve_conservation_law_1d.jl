using LinearAlgebra, SparseArrays

function solve_conservation_law_1d(μ, Ω, Nh, BC, T, δt)
        # ==============================================================================
        # 1. Problem Parameters
        # ==============================================================================

        # Boundaries and boundary conditions 
        a = Ω[1]
        b = Ω[2]
        L = b - a
        a_bc = BC[1]
        b_bc = BC[2]

        # Spatial discretization
        n_nodes = Nh + 1

        # Time discretization
        n_time_steps = Int(T/δt)

        # Parametrised advection field
        ɸ(μ, t) = μ*cos(pi * t)

        # Spatial discretization parameters
        nodes = LinRange(-1.0, 1.0, n_nodes)
        internal_indices = 2:(n_nodes - 1)
        n_int = length(internal_indices)
        h = nodes[2] - nodes[1] # Element length

        # Initial Condition
        u0_func(x) = sech(20.0 * x)
        U0_vector = u0_func.(nodes)
        U0_int = U0_vector[internal_indices]

        # ==============================================================================
        # 2. Local and Global Matrix Assembly
        # ==============================================================================

        # Local element mass and advection matrices
        # For piecewise linear basis functions on a local element [0, h]
        # Local mass matrix M_local = integral of Phi_i * Phi_j
        M_local = h / 6 * [2 1; 1 2]

        # Local advection matrix A_local = integral of Phi_i * Phi_j'
        # Phi_1'(x) = -1/h, Phi_2'(x) = 1/h
        # A_local = [ -1/2  1/2 ]
        #           [ -1/2  1/2 ]
        A_local = [ -0.5 0.5; -0.5 0.5 ]

        # Global sparse mass and advection matrices
        M_global = spzeros(n_nodes, n_nodes)
        A_global = spzeros(n_nodes, n_nodes)

        # Loop over elements and assemble the global matrices
        for i in 1:n_nodes-1
                global_indices = [i, i+1]
                M_global[global_indices, global_indices] .+= M_local
                A_global[global_indices, global_indices] .+= A_local
        end

        # Partition the global matrices for the internal nodes
        M_int = M_global[internal_indices, internal_indices]
        A_int_int = A_global[internal_indices, internal_indices]

        # ==============================================================================
        # 3. Offline Stage: Snapshot Generation and SVD
        # ==============================================================================

        # Snapshot matrix stores the internal nodal values at each time step
        snapshot_matrix = Matrix{Float64}(undef, n_int, n_time_steps + 1)
        snapshot_matrix[:, 1] = U0_int

        # Full-order solution vector for the offline stage
        global U_current_full = copy(U0_int)

        for n in 1:n_time_steps
                t_next = n*δt
                # LHS matrix for the full-order linear system
                LHS_full = M_int - δt*ɸ(μ, t_next)*A_int_int
    
                # RHS vector
                RHS_full = M_int*U_current_full

                # Solve the full-order linear system
                U_next_full = LHS_full\RHS_full
    
                snapshot_matrix[:, n+1] = U_next_full
                U_current_full = U_next_full
        end

        # ==============================================================================
        # 6. Data-handling 
        # ==============================================================================
        
        output = (snapshot_matrix, M_int, A_int_int, U0_int)

        return output
end
