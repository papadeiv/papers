include("../../../../inc/IO.jl")
include("../../../../inc/TimeseriesAnalysis.jl")
include("../../../../inc/EarlyWarningSignals.jl")

import .EarlyWarningSignals as ews

# Read in the timeseries of the non-stationary process
solution = readin("../data/slide_03/solution.csv")
t = solution[:,1]
u = solution[:,3]

# Compute the variance and skewness of the non-stationary process
width = 0.1::Float64
t_var, u_var = ews.variance(t, u, width)

writeout(hcat(t_var, u_var), "../data/slide_03/ews.csv")
