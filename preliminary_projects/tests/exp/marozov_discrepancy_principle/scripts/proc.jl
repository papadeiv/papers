"""
    Postprocessing script

Collection of quantities and functions used to postprocess and analyse the results of a simulation.
"""

# Scalar potential of the conservative system 
U(x, μ) =  + μ*x + (1.0/3.0)*x^3                # Potential (ground truth)
V(x, c) = c[1]*x + c[2]*(x^2) + c[3]*(x^3)      # Potential

# Solve the LLS problem
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

        # Compute the misfit between the model's prediction and the observations
        r = norm(Y - [V(x, θ) for x in Xn])/norm(Y)

        return (
                solution = θ,
                misfit = r
               ) 
end
