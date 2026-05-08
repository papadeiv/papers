"""
    Postprocessing script

Collection of quantities and functions used to postprocess and analyse the results of a simulation.
"""

# Define the (relative) sliding window size
window_size = 0.5::Float64

# Build the bifurcation diagram
function compute_bif_diag(initial_condition)
        # Redefine the vector field for the API 
        F(x, μ) = @. f(x, μ[1])

        # Define the 0-problem
        zero_problem = BifurcationProblem(F, initial_condition, c0, 1; record_from_solution = (x, p; k...) -> x = x[1])

        # Compute the branches of the diagram
        bifurcation_diagram = continuation(zero_problem, PALC(), ContinuationPar(p_min = c0, p_max = cf, dsmax = 0.01, max_steps=1000))

        # Store the bifurcation diagram in the appropriate data structure
        branch_1 = Matrix{Float64}(undef, length(bifurcation_diagram.sol), 3)
        for (index, (point, eigenvalue)) in enumerate(zip(bifurcation_diagram.sol, bifurcation_diagram.eig))
                branch_1[index, 1] = point.p
                branch_1[index, 2] = (point.x)[1]
                branch_1[index, 3] = real((eigenvalue.eigenvals)[1])
        end

        return (
                branch_1 = branch_1 
               )
end

# Converts the sliding window problem into an ensemble one 
function preprocess_solution(timestamps, timeseries, width)
        # Assemble the sliding window
        window = build_window(length(timeseries), width)
        Nw = window.size 
        Ns = window.strides

        # Convert the sliding window subseries into an ensemble of timeseries
        printstyled("Converting the truncated sample path to an ensemble of ", Ns," trajectories of ", Nw, " steps\n"; bold=true, underline=true, color=:light_blue)
        timesteps = [@view timestamps[n:(n+Nw-1)] for n in 1:Ns] 
        ensemble = [@view timeseries[n:(n+Nw-1)] for n in 1:Ns] 

        # Export the ensemble problem 
        return (
                timesteps = timesteps,
                trajectories = ensemble 
               ) 
end
