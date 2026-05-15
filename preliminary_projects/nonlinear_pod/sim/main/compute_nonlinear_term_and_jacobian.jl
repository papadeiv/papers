# ==============================================================================
# Helper Functions for Nonlinear Term and its Jacobian
# ==============================================================================

# This function computes the nonlinear term N(U) = - integral(phi_i * u * u_x) dx
# and its Jacobian J_N = dN/dU.
# It operates on full-order solution vector `U_full_approx`.
# U_full_approx is assumed to be the vector of internal (non-boundary) nodal values.
# n_nodes_full: Total number of nodes in the full system (e.g., 101)
# internal_indices_map: The range or array mapping internal indices to global indices (e.g., 2:100)
function compute_nonlinear_term_and_jacobian(U_full_approx, h, n_nodes_full, internal_indices_map)
    n_int_current = length(U_full_approx)
    N_vector = zeros(n_int_current)
    J_N_matrix = spzeros(n_int_current, n_int_current)

    # Reconstruct full solution including boundary values (which are 0 for this problem)
    # This U_full_with_bc will have size n_nodes_full (101)
    U_full_with_bc = zeros(n_nodes_full)
    U_full_with_bc[internal_indices_map] = U_full_approx # Explicitly map internal values

    for e in 1:(n_nodes_full - 1) # Loop over elements
        # Global indices of nodes for this element
        global_idx_m = e
        global_idx_mp1 = e + 1

        # Get nodal values on this element (U_m, U_{m+1}) from U_full_with_bc
        U_m = U_full_with_bc[global_idx_m]
        U_mp1 = U_full_with_bc[global_idx_mp1]

        # Local derivative u_x on this element (constant for linear elements)
        u_x_e = (U_mp1 - U_m) / h
        
        # Map global node indices to internal node indices for N_vector and J_N_matrix
        # If a global_idx is a boundary node, its corresponding local_idx will be 0.
        local_idx_m = (global_idx_m == 1 || global_idx_m == n_nodes_full) ? 0 : (global_idx_m - 1)
        local_idx_mp1 = (global_idx_mp1 == 1 || global_idx_mp1 == n_nodes_full) ? 0 : (global_idx_mp1 - 1)

        # Contribution to N_vector (from element e to global internal nodes)
        # N_i = - integral(phi_i * u * u_x) dx
        # For linear elements, integral(phi_1 * u * u_x) = (U_mp1 - U_m)/6 * (2*U_m + U_mp1)
        # integral(phi_2 * u * u_x) = (U_mp1 - U_m)/6 * (U_m + 2*U_mp1)
        
        # Add contributions to N_vector
        if local_idx_m != 0
            N_vector[local_idx_m] += - ((U_mp1 - U_m)/6) * (2*U_m + U_mp1)
        end
        if local_idx_mp1 != 0
            N_vector[local_idx_mp1] += - ((U_mp1 - U_m)/6) * (U_m + 2*U_mp1)
        end

        # Contributions to Jacobian J_N = dN/dU
        # dN_i/dU_j terms for local element
        # dN_m/dU_m = -1/6 * (-4*U_m + U_mp1)
        # dN_m/dU_mp1 = -1/6 * (U_m + 2*U_mp1)
        # dN_{m+1}/dU_m = -1/6 * (-2*U_m - U_mp1)
        # dN_{m+1}/dU_mp1 = -1/6 * (-U_m + 4*U_mp1)

        # Assemble into sparse Jacobian matrix J_N_matrix
        # Only add if both indices are internal nodes
        if local_idx_m != 0 && local_idx_m != 0
            J_N_matrix[local_idx_m, local_idx_m] += -1/6 * (-4*U_m + U_mp1)
        end
        if local_idx_m != 0 && local_idx_mp1 != 0
            J_N_matrix[local_idx_m, local_idx_mp1] += -1/6 * (U_m + 2*U_mp1)
        end
        if local_idx_mp1 != 0 && local_idx_m != 0
            J_N_matrix[local_idx_mp1, local_idx_m] += -1/6 * (-2*U_m - U_mp1)
        end
        if local_idx_mp1 != 0 && local_idx_mp1 != 0
            J_N_matrix[local_idx_mp1, local_idx_mp1] += -1/6 * (-U_m + 4*U_mp1)
        end
    end
    return N_vector, J_N_matrix
end
