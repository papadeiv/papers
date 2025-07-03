include("../../../../inc/IO.jl")
include("../../../../inc/TimeseriesAnalysis.jl")
include("../../../../inc/EarlyWarningSignals.jl")

import .EarlyWarningSignals as ews

# Read in the timeseries of the non-stationary process
solution = readin("../data/slide_05/solution.csv")
t = solution[:,1]
u = solution[:,3]

# Compute the variance of the non-stationary process
t_var, u_var = ews.variance(t, u, 0.2)

writeout(hcat(t_var, u_var), "../data/slide_05/ews_variance.csv")
