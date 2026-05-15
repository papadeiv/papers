# This script solves the 1D heat equation with Dirichlet boundary conditions
# using a Finite Element Method for spatial discretization and Implicit Euler
# for time-stepping.

using LinearAlgebra, SparseArrays

function solve_heat_equation_1d(μ, Ω, Nh, BC, T, δt)
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

        # ==============================================================================
        # 2. Mesh and Basis Functions
        # ==============================================================================

        # Domain partition
        nodes = LinRange(a, b, n_nodes)
        h = nodes[2] - nodes[1]

        # Basis function properties:
        # The local basis functions (hat functions) are defined on a local element [0, h]
        # Phi_1(x) = (h-x)/h, Phi_2(x) = x/h
        # Local derivatives: Phi_1'(x) = -1/h, Phi_2'(x) = 1/h

        # Local element mass and stiffness matrices
        # Local mass matrix M_local = integral of Phi_i * Phi_j over [0,h]
        # M_local = h/6 * [2  1]
        #                 [1  2]
        M_local = (h/6)*[2 1; 1 2]

        # Local stiffness matrix K_local = integral of Phi_i' * Phi_j' over [0,h]
        # K_local = 1/h * [ 1 -1]
        #                 [-1  1]
        K_local = (1/h)*[1 -1; -1 1]

        # ==============================================================================
        # 3. Global Matrix Assembly
        # ==============================================================================

        # Global sparse mass and stiffness matrices
        M_global = spzeros(n_nodes, n_nodes)
        K_global = spzeros(n_nodes, n_nodes)

        # Loop over elements and assemble the global matrices
        for i in 1:Nh
                # The global indices for this element are i and i+1
                global_indices = [i, i+1]
    
                # Add local matrices to global matrices
                M_global[global_indices, global_indices] .+= M_local
                K_global[global_indices, global_indices] .+= μ.*K_local
        end

        # ==============================================================================
        # 4. Initial Condition and Boundary Condition Handling
        # ==============================================================================

        # Define a non-trivial initial condition function
        # u0_func(x) = exp(-100*(x^2)) + (a_bc/2)*(1-x) + (b_bc/2)*(1+x) 
        u0_func(x) = ((a_bc/2)*(1-x) + (b_bc/2)*(1+x)) + sin(3*pi*x) 

        # The linear part (x+1)*(b-a)/2 + a is added to satisfy BCs
        # The sin part is the non-trivial IC part

        # Project the initial condition onto the FEM basis to get the initial solution vector U0
        # This is done by solving M_global * U0 = b0, where b0_i = integral of u0 * phi_i
        # We approximate the integral using a simple quadrature rule (summing u0 at nodes)
        # This is a simplification; for a more accurate projection, one would solve the system
        U0_vector = u0_func.(nodes)

        # The solution will be for the internal nodes only (from index 2 to n_nodes-1)
        internal_indices = 2:(n_nodes - 1)
        global U_current = U0_vector[internal_indices]

        # Partition the global matrices for the internal nodes
        M_int = M_global[internal_indices, internal_indices]
        K_int_int = K_global[internal_indices, internal_indices]
        K_int_bc = K_global[internal_indices, [1, n_nodes]]

        # The matrix for the linear system (LHS) is constant, so we compute it once.
        A = M_int + δt*K_int_int

        # ==============================================================================
        # 5. Time-stepping Loop (Implicit Euler)
        # ==============================================================================

        # Storage for the full solution at each time step
        all_solutions = [U0_vector[internal_indices]]

        # The time loop
        for n in 1:n_time_steps
                # Right-hand side (RHS) of the linear system
                b = M_int*U_current - δt*K_int_bc*BC
    
                # Solve for the new solution vector U_next using backslash operator
                U_next_int = A\b
    
                # Store the solution and update for the next time step
                push!(all_solutions, U_next_int)
                U_current = U_next_int
        end

        # ==============================================================================
        # 6. Data-handling 
        # ==============================================================================
        
        output = (all_solutions, M_int, K_int_int, K_int_bc, U0_vector[internal_indices])

        return output 
end
