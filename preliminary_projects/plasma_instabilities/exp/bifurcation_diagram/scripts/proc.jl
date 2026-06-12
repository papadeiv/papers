"""
    Postprocessing script

Collection of quantities and functions used to postprocess and analyse the results of a simulation.
"""

# Build the bifurcation diagram
function compute_bif_diag(initial_condition)
        # Redefine the vector field for the API 
        F(x, μ) = @. f(x, μ[1])

        # Define the 0-problem
        zero_problem = BifurcationProblem(F, [1.0,-1.0], μ_set[1], 1; record_from_solution = (x, p; k...) -> x = x[1])

        # Compute the branches of the diagram
        bifurcation_diagram = continuation(zero_problem, PALC(), ContinuationPar(p_min = -μ_set[end], p_max = μ_set[end], dsmax = 0.01, max_steps=1000))

        # Store the bifurcation diagram in the appropriate data structure
        branch = Matrix{Float64}(undef, length(bifurcation_diagram.sol), 3)
        for (index, (point, eigenvalue)) in enumerate(zip(bifurcation_diagram.sol, bifurcation_diagram.eig))
                branch[index, 1] = point.p
                branch[index, 2] = (point.x)[1]
                branch[index, 3] = real((eigenvalue.eigenvals)[1])
        end

        return (
                branch = branch,
                bifurcations = Float64[] 
               )
end
