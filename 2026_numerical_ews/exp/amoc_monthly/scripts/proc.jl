"""
    Postprocessing script

Collection of quantities and functions used to postprocess and analyse the results of a simulation.
"""

# Parameters of the scalar potential method
window_size = 0.1                               # Relative width of the slinding window
idx = 21000                                     # Time index of the tipping point 
Na = convert(Int64, 1e4)                        # Number of attempts per guess 
β = 1e-2                                        # Std of the guess perturbation 

# Arbitrary cubic potential
V(x, c) = c[1]*x + c[2]*(x^2) + c[3]*(x^3)

# Data structures
solutions = Vector{Vector{Float64}}()           # Solutions of the inference method 
time_idx = Vector{Float64}()                    # Index of timesteps 
ews = Vector{Float64}()                         # Statistical distribution over the ensemble

# Converts the non-stationary timeseries into an ensemble of subseries associated to the strides of a sliding window 
function preprocess_solution(timestamps, timeseries, width)
        # Extract the subseries up to the tipping
        t = timestamps[1:idx] 
        u = timeseries[1:idx]
        Nt = length(u)

        # Assemble the sliding window
        window = build_window(Nt, width)
        Nw = window.size 
        Ns = window.strides

        # Convert the sliding window subseries into an ensemble of timeseries
        printstyled("Converting the truncated sample path to an ensemble of ", Ns," trajectories of ", Nw, " steps\n"; bold=true, underline=true, color=:light_blue)
        timesteps = [@view t[n:(n+Nw-1)] for n in 1:Ns] 
        ensemble = [@view u[n:(n+Nw-1)] for n in 1:Ns] 

        # Export the parameters of the ensemble problem 
        return (
                tipping_point = idx,
                timesteps = timesteps,
                trajectories = ensemble 
               ) 
end

function analyse(solutions)
        # Compute estimated stable and unstable equilibria of the cubic
        xs = +(1/(3*solutions[3]))*(sqrt((solutions[2])^2 - 3*solutions[1]*solutions[3]) - solutions[2])
        xu = -(1/(3*solutions[3]))*(sqrt((solutions[2])^2 - 3*solutions[1]*solutions[3]) + solutions[2])

        # Compute the ews
        ΔV = abs(V(xu, solutions) - V(xs, solutions))
        push!(ews, exp(-ΔV))
end
