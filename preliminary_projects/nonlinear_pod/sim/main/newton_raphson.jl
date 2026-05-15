# ==============================================================================
#  Nonlinear Solver for Time-Stepping (Newton-Raphson)
# ==============================================================================

# This function solves a nonlinear system G(U) = 0 using Newton's method.
# It is designed to be generic enough to handle both the Full-Order Model (FOM)
# and the Reduced-Order Model (ROM).
# Arguments:
#   U_initial: Initial guess for the solution at the current time step (U_k in Newton's loop)
#   U_previous: Solution from the previous time step (U^n)
#   M: Mass matrix (M_int for FOM, Mr for ROM)
#   K_base: Base stiffness matrix (K_int_int_base for FOM, Kr_base for ROM)
#   mu_val: Current value of the viscosity parameter mu
#   dt: Time step size
#   n_nodes_full: Total number of nodes in the full system (needed for compute_nonlinear_term_and_jacobian)
#   internal_indices_map: The range or array mapping internal indices to global indices
#   Phi: POD basis matrix (optional, default is Identity for FOM). For ROM,
#        this is used to reconstruct the full-order solution U_int from U_r.

include("./compute_nonlinear_term_and_jacobian.jl")

const NEWTON_MAX_ITER = 15 # Increased iterations for potentially tougher nonlinearity
const NEWTON_TOL = 1e-7    # Convergence tolerance

function newton_raphson(U_initial, U_previous, M, K_base, mu_val, dt, h, n_nodes_full, internal_indices_map; Phi=I)
    U_k = copy(U_initial) # Current guess for U^{n+1} (U_r for ROM, U_int for FOM)
    
    # Linear part of the Jacobian: (M + dt * K)
    # K = mu_val * K_base
    K_linear = mu_val * K_base
    J_linear_part = M + dt * K_linear

    for k in 1:NEWTON_MAX_ITER
        # Reconstruct the full-order solution approximation from U_k.
        # If FOM, Phi is Identity, so U_full_approx = U_k.
        # If ROM, U_full_approx = Phi * U_k.
        U_full_approx = Phi * U_k 

        # Compute the full-order nonlinear term N(U_full_approx) and its Jacobian J_N_full
        # Note: mu_val is not used in compute_nonlinear_term_and_jacobian for N_vector and J_N_full,
        # as the -u*u_x term does not depend on viscosity.
        N_full_approx, J_N_full = compute_nonlinear_term_and_jacobian(U_full_approx, h, n_nodes_full, internal_indices_map)
        
        # Project N_full_approx back to the current space (ROM or FOM)
        # N_k_projected = Phi' * N(U_full_approx) for ROM, or N_k_projected = N(U_full_approx) for FOM
        N_k_projected = Phi' * N_full_approx

        # Evaluate the nonlinear function G(U_k) for the current system (FOM or ROM)
        # G(U) = (M + dt*K)*U - M*U_previous - dt*N(U)
        G_U_k = J_linear_part * U_k - M * U_previous - dt * N_k_projected
        
        # Check for convergence: ||G(U_k)|| < tolerance
        if norm(G_U_k) < NEWTON_TOL
            return U_k # Converged solution
        end

        # Assemble the Jacobian matrix J for the current system (FOM or ROM)
        # J = J_linear_part - dt * (Phi' * J_N_full * Phi)
        J = J_linear_part - dt * (Phi' * J_N_full * Phi)
        
        # Solve the linear system J*delta_U = -G
        delta_U = J \ (-G_U_k)
        
        # Update the solution: U_{k+1} = U_k + delta_U
        U_k .+= delta_U
    end
    println("Newton-Raphson did not converge after $NEWTON_MAX_ITER iterations.")
    return U_k # Return best guess if not converged
end
