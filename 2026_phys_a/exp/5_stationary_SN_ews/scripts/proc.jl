"""
    Postprocessing script

Collection of quantities and functions used to postprocess and analyse the results of a simulation.
"""

# Scalar potential of the conservative system 
U(x, μ) =  + μ*x + (1.0/3.0)*x^3                # Potential (ground truth)
V(x, c) = c[1]*x + c[2]*(x^2) + c[3]*(x^3)      # Potential

# Stationary solution of the FPE
ρ(x, μ) = exp(-U(x, μ)/D)

# Build the bifurcation diagram
function compute_bif_diag()
        # Redefine the vector field for the API 
        F(x, μ) = @. f(x, μ[1])

        # Define the 0-problem
        zero_problem = BifurcationProblem(F, [sqrt(-(-1.0))], -1.0, 1; record_from_solution = (x, p; k...) -> x = x[1])

        # Compute the branches of the diagram
        bifurcation_diagram = continuation(zero_problem, PALC(), ContinuationPar(p_min = -1.0, p_max = 0.1, dsmax = 0.01, max_steps=1000))

        # Store the bifurcation diagram in the appropriate data structure
        branch = Matrix{Float64}(undef, length(bifurcation_diagram.sol), 3)
        for (index, (point, eigenvalue)) in enumerate(zip(bifurcation_diagram.sol, bifurcation_diagram.eig))
                branch[index, 1] = point.p
                branch[index, 2] = (point.x)[1]
                branch[index, 3] = real((eigenvalue.eigenvals)[1])
        end

        # Store the bifurcation points
        bifurcations = Matrix{Float64}(undef, 2, 2) 
        for (index, point) in enumerate(bifurcation_diagram.specialpoint)
                if point.type == :bp
                        bifurcations[index, 1] = point.param
                        bifurcations[index, 2] = point.norm
                end
        end

        return (
                branch = branch,
                bifurcations = bifurcations
               )
end

# Solve the LLS problem
function solve_lls(solution)
        # Define the observation vectors 
        Xn = solution[1:end-1]
        Y  = (solution[2:end] .- solution[1:end-1])./dt

        # Assemble the model matrix
        A = hcat(ones(length(Xn)), Xn, Xn.^2)

        # Solve the (linear) least-squares problem
        β = A\Y

        # Compute the coefficients of the potential
        θ = [-β[1], -β[2]/2, -β[3]/3]
        return θ 
end
 

function compute_ews(solution)
        # Solve the LLS problem
        θ = solve_lls(solution)

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
