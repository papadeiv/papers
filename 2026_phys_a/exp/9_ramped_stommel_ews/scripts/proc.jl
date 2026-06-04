"""
    Postprocessing script

Collection of quantities and functions used to postprocess and analyse the results of a simulation.
"""

# Scalar potential of the conservative system 
V(x, c) = c[1]*x + c[2]*(x^2) + c[3]*(x^3)

# Define the (relative) sliding window size
window_size = 0.10::Float64

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
 
# Compute the modified escape EWS
function compute_ews(θ)
        # Compute estimated stable and unstable equilibria of the cubic approximation
        xs = +(1/(3*θ[3]))*(sqrt((θ[2])^2 - 3*θ[1]*θ[3]) - θ[2])
        xu = -(1/(3*θ[3]))*(sqrt((θ[2])^2 - 3*θ[1]*θ[3]) + θ[2])

        # Compute the modified escape EWS 
        ΔV = abs(V(xu, θ) - V(xs, θ))
        escape = exp(-ΔV)

        # Return the EWS 
        return escape 
end
