"""
Utilities to analyse the sample paths solutions of stochastic processes whose determinist term is non-autonomus (i.e. time-dependent).

Author: Davide Papapicco
Affil: U. of Auckland
Date: 03-09-2025
"""

#----------------------------#
#                            # 
#   transient_processes.jl   #                     
#                            #
#----------------------------#

"""
    detrend(timestamps, timeseries; kwargs...)

Removes the trend from a `timeseries` and computes the stationary residuals.
"""
function detrend(timeseries; alg = "exact", timestamps = Float64[], qse = Float64[])
        # Initialise arrays for the trend and the residuals
        trend = Vector{Float64}(undef, length(timeseries))
        residuals = Vector{Float64}(undef, length(timeseries))

        # Detrend the timeseries
        if alg == "mean"
                # Compute the mean of the timeseries to define the trend
                trend = mean(timeseries).*ones(length(timeseries))
                # Remove the trend to find the residuals
                residuals = timeseries - trend 

        elseif alg == "linear"
                # Solve the linear least-squares problem 
                coefficients = approximate_potential(timestamps, timeseries, degree=1).coeffs 
                # Define the linear trend
                trend = [coefficients[1] + coefficients[2]*t for t in timestamps]
                # Remove the trend to find the residuals
                residuals = timeseries - trend

        elseif alg == "nonlinear"
                trend = mean(timeseries).*ones(length(timeseries))
                residuals = timeseries - trend

        elseif alg == "emd"
                trend = mean(timeseries).*ones(length(timeseries))
                residuals = timeseries - trend

        elseif alg == "exact"
                trend = qse 
                residuals = timeseries - trend

        else
                display("The algorithm you specified does not exist or it is not implemented for this method")
        end

        return (
                trend = trend, 
                residuals = residuals
               )

        #=
        py"""
        import sys
        sys.path.append('/home/dpap666/Libraries/STEWS/src/python')

        from Process import Process
        from Estimator import Estimator 
        from TimeSeries import TimeSeries

        import numpy as np

        def analyse(data, savepath):
                ts = Process()
                ts.realizations = TimeSeries(realizations=data)
                #ts.detrend(mode='EMD', order=2)
                #imfs = ts.imfs
                #np.savetxt(savepath, imfs[:,(imfs.shape)[1]-1], delimiter=',')
                ts.detrend(mode='diff')
                np.savetxt(savepath, ts.detrended.ts, delimiter=',')

        """

        pyanalyse = py"analyse"
        pyanalyse(timeseries, savepath)
        =#
end

"""
    get_window_parameters(Nt, width)

Compute the number of timesteps in a window of `width` relative size w.r.t. a total of `Nt` timesteps and the number of strides.
"""
function get_window_parameters(Nt::Int64, width::Float64)
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

"""
    find_tipping(ut; check=0.1, criterion=0.1)

Identify the tipping point in a timeseries `u` by checking the deviation of the last `criterion` realisations away from the mean of the first `check` realisations.
"""
function find_tipping(ut::Vector{Float64}; check = 0.100::Float64, criterion = 0.010::Float64)
        # Get the length of the timeseries
        Nt = length(ut)
        
        # Get the number of realisations in a window to check for tipping
        window = get_window_parameters(Nt, check)
        Nw = window.size
        Ns = window.strides

        # Compute the split index between the control and the target subseries
        Np = convert(Int64, floor(criterion*Nw))

        # Initialise the tipping point flags
        tip_chk = false
        tip_idx = Nt

        # Loop over the number of strides of the sliding window
        printstyled("Searching for tipping points\n"; bold=true, underline=true, color=:light_blue)
        @showprogress for n in 1:Ns
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

        if tip_chk
                # Print this message if tipping has been bound
                printstyled("Tipping point found at timestep ", tip_idx, " of ", Nt,"\n"; bold=true, underline=true, color=:green)
          
        else
                # Print this message if not tipping has been bound
                printstyled("No tipping point found\n"; bold=true, underline=true, color=:red)
        end

        # Return the boolean check and the tip index
        return (
                check = tip_chk,
                index = tip_idx
               )
end
