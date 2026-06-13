"""
    Postprocessing script

Collection of quantities and functions used to postprocess and analyse the results of a simulation.
"""

# Data structure to store bifurcation points
bifurcations = Matrix{Float64}(undef, 2, 3)

# Build the bifurcation diagram
function compute_bif_diag(initial_condition, μ_range)

        zero_problem = BifurcationProblem(F, initial_condition, μ, (@optic _.η2); record_from_solution = (x, p; k...)->(T = x[1], S = x[2]),) 

        bifurcation_diagram = continuation(zero_problem, PALC(), ContinuationPar(p_min = μ_range[1], p_max = μ_range[2], dsmin = 1e-6, dsmax = 1e-2, max_steps=10000); bothside = true)

        branch = Matrix{Float64}(undef, length(bifurcation_diagram.sol), 4)
        for (index, (point, eigenvalue)) in enumerate(zip(bifurcation_diagram.sol, bifurcation_diagram.eig))
                branch[index, 1] = point.p
                branch[index, 2] = (point.x)[1]
                branch[index, 3] = (point.x)[2]
                branch[index, 4] = maximum(real(eigenvalue.eigenvals))
        end

        # Store the bifurcation point
        bifurcations = Vector{Float64}(undef, 3)
        for (index, point) in enumerate(bifurcation_diagram.specialpoint)
                if point.type == :bp
                        bifurcations[1] = point.param
                        bifurcations[2] = (point.x)[1]
                        bifurcations[3] = (point.x)[2]
                end
        end

        return (
                branch = branch,
                bifurcations = bifurcations
               )
end
