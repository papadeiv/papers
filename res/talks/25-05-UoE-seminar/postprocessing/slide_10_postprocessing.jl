include("../../../../inc/IO.jl")
include("../../../../inc/EscapeProblem.jl")
include("../../../../inc/SystemAnalysis.jl")
include("../../../../inc/PotentialLearning.jl")
include("../../../../inc/TimeseriesAnalysis.jl")

# Import the data from csv 
OUP = readin("../data/slides_10_to_13/OUP.csv")
solution = readin("../data/slides_10_to_13/solution.csv")
t = solution[:,1]
μ = solution[:,2]
u = solution[:,3]
Nt = length(t)

# Define the width of the sliding window (as a relative length of the timeseries)
width = 0.150::Float64

# Assemble the sliding window
window = get_window_parameters(Nt, width)
Nw = window[1]
Ns = window[2]

# Number of bins for the histograms across the sliding window (2% of the number of timesteps in the window)
Nb = convert(Int64, floor(0.02*Nw))
# Number of bins for the histograms of the OUP equivalent (2% of the number of timesteps of the OUP)
Nb_OUP = 20 

# Loop over the window's strides
printstyled("Detrending and binning the timeseries\n"; bold=true, underline=true, color=:light_blue)
@showprogress for n in 1:Ns
        # Define the index of the timestep at the end of the sliding window
        n_end = Nw + n - 1::Int64

        # Detrend the timeseries by removing the mean
        u_det = detrend(t[n:n_end], u[n:n_end], alg = "mean")
        mean_trend = u_det[1]
        mean_residuals = u_det[2]

        # Detrend the timeseries by removing the mean
        u_det = detrend(t[n:n_end], u[n:n_end], alg = "linear")
        lin_trend = u_det[1]
        lin_residuals = u_det[2]

        # Fit an empirical distribution to the detrended timeseries
        mean_bins, mean_pdf = fit_distribution(mean_residuals, n_bins=Nb+1)
        lin_bins, lin_pdf = fit_distribution(lin_residuals, n_bins=Nb+1)

        # Fit a coarser histogram to the OUP equivalent
        OUP_bins, OUP_pdf = fit_distribution(OUP[:,n_end], n_bins=Nb_OUP+1)

        # Export the data 
        writeout(hcat(t[n:n_end], μ[n:n_end], u[n:n_end], mean_trend, mean_residuals, lin_trend, lin_residuals), "../data/slide_10/residuals/$n.csv")
        writeout(hcat(OUP_bins, OUP_pdf), "../data/slide_10/histograms/OUP/$n.csv")
        writeout(hcat(mean_bins, mean_pdf), "../data/slide_10/histograms/SN_mean/$n.csv")
        writeout(hcat(lin_bins, lin_pdf), "../data/slide_10/histograms/SN_lin/$n.csv")
end
