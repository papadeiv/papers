# Computes the temporal variance of a timeseries over a sliding window of specified width

function skew(time, data, width::Float64)
        # Get the total number of steps in the timeseries 
        Nt = length(time)

        # Assemble the sliding window
        window = get_window_parameters(Nt, width)
        Nw = window[1]
        Ns = window[2]

        # Get the truncated timesteps
        time_ews = time[Nw:end] 

        # Compute the variance ews across the sliding window
        printstyled("Computing the skewness EWS across the sliding window\n"; bold=true, underline=true, color=:light_blue)
        ews = [skewness(data[n:(n + Nw - 1)]) for n in 1:Ns]

        # Return the ews
        return time_ews, ews
end
