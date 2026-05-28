"""
    Postprocessing script

Collection of quantities and functions used to postprocess and analyse the results of a simulation.
"""

# Scalar potential of the conservative system 
U(x, μ) =  + μ*x + (1.0/3.0)*x^3                # Potential (ground truth)
V(x, c) = c[1]*x + c[2]*(x^2) + c[3]*(x^3)      # Potential

# Stationary solution of the FPE
ρ(x, μ) = exp(-U(x, μ)/D)

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
        return θ 
end
 
# Compute the variance and modified escape EWS
function compute_ews(solution; α = 1e-2)
        # Solve the LLS problem
        θ = solve_lls(solution, α=α)

        # Compute estimated stable and unstable equilibria of the cubic approximation
        xs = +(1/(3*θ[3]))*(sqrt((θ[2])^2 - 3*θ[1]*θ[3]) - θ[2])
        xu = -(1/(3*θ[3]))*(sqrt((θ[2])^2 - 3*θ[1]*θ[3]) + θ[2])

        # Compute the modified escape EWS 
        ΔV = abs(V(xu, θ) - V(xs, θ))
        escape = exp(-ΔV)

        # Compute the sample variance
        variance = var(solution)

        # Return the EWS 
        return (
                escape = escape,
                variance = variance
               )
end
