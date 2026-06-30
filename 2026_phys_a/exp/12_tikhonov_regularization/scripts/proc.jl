"""
    Postprocessing script

Collection of quantities and functions used to postprocess and analyse the results of a simulation.
"""

# Solve the RLLS problem
function solve_lls(solution; α = 1e-2)
        # Define the observation vectors 
        Xn = solution[1:end-1]
        Y  = (solution[2:end] .- solution[1:end-1])./dt

        # Assemble the model matrix
        A = hcat(ones(length(Xn)), Xn, Xn.^2)

        # Solve the regularized (linear) least-squares problem
        β = (A'*A + α.*I(size(A,2)))\(A'*Y)

        # Compute the coefficients of the potential
        θ = [-β[1], -β[2]/2, -β[3]/3]
        return θ 
end
 
