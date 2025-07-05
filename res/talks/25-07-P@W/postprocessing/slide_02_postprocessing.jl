include("../../../../inc/IO.jl")
include("../../../../inc/PotentialLearning.jl")

# Import data from csv 
solution = readin("../data/slide_02/solution.csv")
t = solution[:,1]
u = solution[:,2]
Nt = length(t)

# Define number of bins for the histogram (0.25% of the total number of timesteps)
Nb = convert(Int64, floor(0.00005*Nt))

# Define the domain for the histogram
domain = [minimum(u), maximum(u)]

# Compute the histogram of the process upto the current timestep
bins, pdf = fit_distribution(u, interval = domain, n_bins=Nb+1)

# Export the data 
writeout(hcat(bins, pdf), "../data/slide_02/histogram.csv")
