include("../../../../inc/IO.jl")
include("../../../../inc/TimeseriesAnalysis.jl")
include("../../../../inc/EarlyWarningSignals.jl")

# Read in the timeseries of the stationary process
solution = readin("../data/slide_04/N-tipping.csv")
t = solution[:,1]
u = solution[:,2]
Nt = length(t)

# Define number of bins for the histogram (0.25% of the total number of timesteps)
Nb = convert(Int64, floor(0.0025*Nt))

# Define the domain for the histogram
domain = [minimum(u), maximum(u)]

# Define an empty vector to store the bins of the histograms
global bins = Vector{Float64}(undef, Nb)
# Define an empty matrix to store the histogram at different timesteps
global pdf = Matrix{Float64}(undef, Nb, Nt)

# Loop over the realizations
printstyled("Postprocessing the N-tipping solution using $(Threads.nthreads()) threads\n"; bold=true, underline=true, color=:darkgreen)
@showprogress Threads.@threads for n in 1:length(t)
        # Extract the solution upto the current timestep
        ut = u[1:n]

        # Compute the histogram of the process upto the current timestep
        bins[:], pdf[:,n] = fit_distribution(ut, interval = domain, n_bins=Nb+1)
end

# Export the bins and histograms
writeout(bins, "../data/slide_04/bins.csv")
writeout(pdf, "../data/slide_04/pdf.csv")
