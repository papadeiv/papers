include("../../../../inc/IO.jl")
include("../../../../inc/TimeseriesAnalysis.jl")

# Read in the timeseries of the non-stationary process
u = readin("../data/slide_06/u.csv")

# Detrend the mean to make the process weakly stationary
detrend(u, "../data/slide_06/ut.csv")
