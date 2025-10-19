"""
Utilities to analyse the sample paths solutions of stochastic processes whose determinist term is non-autonomus (i.e. time-dependent).

Author: Davide Papapicco
Affil: U. of Auckland
Date: 03-09-2025
"""

function detrend(timeseries; alg = "exact", timestamps = Float64[], qse = Float64[])
        # Initialise arrays for the trend and the residuals
        trend = Float64[]
        residuals = Float64[]

        # Detrend the timeseries
        if alg == "mean"
                # Compute the mean of the timeseries to define the trend
                trend = mean(timeseries).*ones(length(timeseries))
                # Remove the trend to find the residuals
                residuals = timeseries - trend 

        elseif alg == "linear"
                # Assemble the model matrix
                X = hcat(ones(length(timestamps)), timestamps)
                # Solve the least-squares problem
                c = X\timeseries
                # Compute the linear trend and the residuals
                trend = X*c
                residuals = timeseries - trend

        else alg == "exact"
                trend = qse 
                residuals = timeseries - trend
        end

        return (
                trend = trend, 
                residuals = residuals
               )
end

function build_window(Nt::Int64, width::Float64)
        # Compute the size Nw of the window
        Nw = convert(Int64, floor(width*Nt))

        # Compute the number Ns of subseries (the number of strides is Ns - 1)
        Ns = (Nt - Nw) + 1::Int64

        # Return the window parameters
        return (
                size = Nw, 
                strides = Ns 
               )
end

function find_tipping(ut::AbstractVector{Float64}; check = 0.100::Float64, criterion = 0.010::Float64, verbose=true)
        # Get the length of the timeseries
        Nt = length(ut)
        
        # Get the number of realisations in a window to check for tipping
        window = build_window(Nt, check)
        Nw = window.size
        Ns = window.strides

        # Compute the split index between the control and the target subseries
        Np = convert(Int64, floor(criterion*Nw))

        # Initialise the tipping point flags
        tip_chk = false
        tip_idx = Nt

        # Loop over the number of strides of the sliding window
        for n in 1:Ns
                # Extract the indices of the control and target subseries in the window
                n_end = Nw + n - 1::Int64
                n_split = n_end - Np

                # Split the subseries into the control and the target subseries
                u_ctl = ut[n:n_split]
                u_trg = ut[(n_split+1):end]

                # Compute the mean and std in the control subseries
                u_mean = mean(u_ctl)
                u_std = std(u_ctl)

                # Define the boundaries for tipping
                u_inf = u_mean - 3*u_std
                u_sup = u_mean + 3*u_std

                # Determine whether all the realisations in the target subseries have tipped
                if !any(u -> u_inf <= u <= u_sup, u_trg)
                        # Store the tip flags and stop the loop
                        tip_chk = true
                        tip_idx = n_split 
                        break
                end
        end

        if tip_chk && verbose
                # Print this message if tipping has been bound
                printstyled("Tipping point found at timestep ", tip_idx, " of ", Nt,"\n"; bold=true, underline=true, color=:green)
          
        elseif verbose
                # Print this message if not tipping has been bound
                printstyled("No tipping point found\n"; bold=true, underline=true, color=:red)
        end

        # Return the boolean check and the tip index
        return (
                check = tip_chk,
                index = tip_idx
               )
end
