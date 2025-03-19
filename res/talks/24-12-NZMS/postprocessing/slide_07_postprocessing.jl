include("../../../../inc/IO.jl")
include("../../../../inc/TimeseriesAnalysis.jl")
include("../../../../inc/EarlyWarningSignals.jl")

import .EarlyWarningSignals as ews

# Read in the timeseries of the non-stationary process
t = readin("../data/slide_07/time.csv")
u = readin("../data/slide_07/u.csv")

# Compute the variance of the non-stationary process
t_var, u_var = ews.variance(t, u, 0.1)

writeout(hcat(t_var, u_var), "../data/slide_07/ews_variance.csv")

# Detrend the mean to make the process weakly stationary
detrend(u, "../data/slide_07/ut.csv")

# Compute the variance of the residuals 
