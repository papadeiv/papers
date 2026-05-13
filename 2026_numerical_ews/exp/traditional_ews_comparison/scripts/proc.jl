"""
    Postprocessing script

Collection of quantities and functions used to postprocess and analyse the results of a simulation.
"""

# Scalar potential of the conservative system 
U(x, μ) = -(r/2.0)*x^2 + (r/(3.0*k))*x^3 + μ*(x - sqrt(h)*atan(x/sqrt(h)))    # Ground truth
V(x, c) = c[1]*x + c[2]*(x^2) + c[3]*(x^3)                                    # Arbitrary cubic potential

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
function preprocess_solution(timestamps, timeseries, width; verbose = true)
        # Assemble the sliding window
        window = build_window(length(timeseries), width)
        Nw = window.size 
        Ns = window.strides

        # Convert the sliding window subseries into an ensemble of timeseries
        if verbose
                printstyled("Converting the truncated sample path to an ensemble of ", Ns," trajectories of ", Nw, " steps\n"; bold=true, underline=true, color=:light_blue)
        end
        timesteps = [@view timestamps[n:(n+Nw-1)] for n in 1:Ns] 
        ensemble = [@view timeseries[n:(n+Nw-1)] for n in 1:Ns] 

        # Export the ensemble problem 
        return (
                timesteps = timesteps,
                trajectories = ensemble 
               ) 
end

# Solve the MLE on the EM approximation
function solve_EM_MLE(solution)
        Xn = solution[1:end-1]
        Y  = (solution[2:end] .- solution[1:end-1])./δt
        Φ = hcat(ones(length(Xn)), Xn, Xn.^2)
        β = Φ\Y
        c = [-β[1], -β[2]/2, -β[3]/3]
        return c
end

# Compute the early-warning signal
function compute_ews(solutions)
        # Compute estimated stable and unstable equilibria of the cubic
        xs = +(1/(3*solutions[3]))*(sqrt((solutions[2])^2 - 3*solutions[1]*solutions[3]) - solutions[2])
        xu = -(1/(3*solutions[3]))*(sqrt((solutions[2])^2 - 3*solutions[1]*solutions[3]) + solutions[2])

        # Compute the ews
        ΔV = abs(V(xu, solutions) - V(xs, solutions))
        return exp(-ΔV)
end
