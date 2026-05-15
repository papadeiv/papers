include("./newton_raphson.jl")
using LinearAlgebra, SparseArrays

function solve_nonlinear_burgers_1d(μ, Ω, Nh, BC, T, δt)
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

       # K_local_ij = integral of mu * (dphi_i/dx) * (dphi_j/dx) over [0,h]
        # (dphi_i/dx) is -1/h for phi_1 and 1/h for phi_2
        K_local_base = 1 / h * [1 -1; -1 1] # This is the base for stiffness, will be multiplied by mu

        # Global sparse mass and advection matrices
        M_global = spzeros(n_nodes, n_nodes)
        K_global_base = spzeros(n_nodes, n_nodes)

        # Loop over elements and assemble the global matrices
        for i in 1:n_nodes-1
                global_indices = [i, i+1]
                M_global[global_indices, global_indices] .+= M_local
                K_global_base[global_indices, global_indices] .+= K_local_base
        end

        # Partition the global matrices for the internal nodes
        M_int = M_global[internal_indices, internal_indices]
        K_int_int_base = K_global_base[internal_indices, internal_indices]

        # ==============================================================================
        # 3. Offline Stage: Snapshot Generation and SVD
        # ==============================================================================

        # Snapshot matrix stores the internal nodal values at each time step
        snapshot_matrix = Matrix{Float64}(undef, n_int, n_time_steps + 1)
        snapshot_matrix[:, 1] = U0_int

        # Full-order solution vector for the offline stage
        global U_current_full = copy(U0_int)

        for n in 1:n_time_steps
                # Call the nonlinear solver for the full system.
                # Phi is implicitly Identity (default).
                U_next_full = newton_raphson(U_current_full, U_current_full, M_int, K_int_int_base, μ, δt, h, n_nodes, internal_indices)
    
                snapshot_matrix[:, n + 1] = U_next_full
                global U_current_full = U_next_full
        end

        # ==============================================================================
        # 4. Data-handling 
        # ==============================================================================
        
        output = (snapshot_matrix, M_int, K_int_int_base, U0_int, h)

        return output
end
