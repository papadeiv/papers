"""
Utilities to analyse the sample paths solutions of stochastic processes whose determinist term is non-autonomus (i.e. time-dependent).

Author: Davide Papapicco
Affil: U. of Auckland
Date: 03-09-2025
"""

function detrend(timeseries; alg = "exact", timestamps = Float64[], qse = Float64[], n_modes=0::Integer)
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

        elseif alg == "exact"
                trend = qse 
                residuals = timeseries - trend
        else alg == "emd"
                emd = PyEMD.EMD()
                imfs = Array(emd(timeseries))
                trend = sum(imfs[(end-n_modes):end,:], dims=1)[:]
                residuals = timeseries .- trend
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

function find_tipping(ut::AbstractVector{Float64}; width = 200::Integer, threshold = 2.0::Float64, verbose=true)
        # Extract the length of the timeseries
        Nt = length(ut)

        # Loop over the window strides
        scores = zeros(Nt)
        for n in (width+1):(Nt-width)
                # Compute the mean and standard deviation of the left and right subseries in the window
                μL = mean(@view ut[n-width:n-1])
                μR = mean(@view ut[n+1:n+width])
                σL = std(@view ut[n-width:n-1])
                σR = std(@view ut[n+1:n+width])

                # Compute the score matching
                scores[n] = abs(μR - μL)/(0.5*(σL + σR) + 1e-12)
        end

        # Find the maximum score
        tipping_idx = argmax(scores)

        # Impose the tipping criterion
        if scores[tipping_idx] > threshold
                # Print message
                if verbose
                        printstyled("Tipping point found at timestep ", tipping_idx, " of ", Nt,"\n"; bold=true, underline=true, color=:green)
                end

                # Return the tipping point index 
                return tipping_idx 
        else
                # Print message
                if verbose
                        printstyled("No tipping point found\n"; bold=true, underline=true, color=:red)
                end

                # Return the tipping point index 
                return Nt 
        end
end
