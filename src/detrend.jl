using Statistics, PyCall

# Detrend the timeseries using different algorithms (according user's choice) 
function detrend(timestamps, timeseries; alg = "exact", qse = Float64[])
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

        return [trend, residuals]

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
