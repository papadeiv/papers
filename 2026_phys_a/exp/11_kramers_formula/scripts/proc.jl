"""
    Postprocessing script

Collection of quantities and functions used to postprocess and analyse the results of a simulation.
"""

# Scalar potential of the conservative system and its second derivative 
U(x, μ) = -(r/2.0)*x^2 + (r/(3.0*k))*x^3 + μ*(x - sqrt(h)*atan(x/sqrt(h)))
Uxx(x, μ) = -1 + x/5 + (2*μ)*(x/(x^2 + 1)^2)

# Build the bifurcation diagram
function compute_bif_diag(initial_condition)
        # Redefine the vector field for the API 
        F(x, μ) = @. f(x, μ[1])

        # Define the 0-problem
        zero_problem = BifurcationProblem(F, initial_condition, c0, 1; record_from_solution = (x, p; k...) -> x = x[1])

        # Compute the branches of the diagram
        bifurcation_diagram = continuation(zero_problem, PALC(), ContinuationPar(p_min = c0, p_max = cf, dsmax = 0.01, max_steps=1000))

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
