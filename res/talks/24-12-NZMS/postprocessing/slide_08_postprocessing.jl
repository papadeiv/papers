include("../../../../inc/IO.jl")
include("../../../../inc/TimeseriesAnalysis.jl")
include("../../../../inc/EarlyWarningSignals.jl")

import .EarlyWarningSignals as ews

# Read in the timeseries of lattice simulation 
t = readin("../data/slide_08/time.csv")
μ = readin("../data/slide_08/K.csv")
u = readin("../data/slide_08/u.csv")

# Get number of lattice sites
Ns = size(u)[1]

# Get number of timesteps
Nt = size(u)[2]

# Compute the mean field timeseries (spatial aggregate)
u_mean = [(1.0::Float64/Ns)*sum(u[:,t]) for t=1:Nt] 

# Detrend the mean to make the process weakly stationary
detrend(u_mean, "../data/slide_08/detrended_mean_field.csv")
ut_mean = readin("../data/slide_08/detrended_mean_field.csv")

# Compute the temporal variance of the detrended mean field timeseries
μ_var, u_var = ews.variance(μ, ut_mean, 0.1)

# Export the data
writeout(u_mean, "../data/slide_08/mean_field.csv")
writeout(hcat(μ_var, u_var), "../data/slide_08/ews_variance.csv")
