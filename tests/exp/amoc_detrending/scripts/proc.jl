"""
    Postprocessing script

Collection of quantities and functions used to postprocess and analyse the results of a simulation.
"""

# Parameters of the sliding window 
window_size = 0.4                               # Relative width of the slinding window
idx = 21000                                     # Time index of the tipping point 

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
