"""
    Main script

Run this file to execute the simulation, analyse and plot the results.
"""

# Import the necessary packages and local modules
include("inc.jl")

# Import the simulation's scripts
include("./scripts/proc.jl")

# Define the main algorithm
function main()
        # Data structures
        solutions = Vector{Vector{Float64}}()           # Solutions of the inference method 
        time_idx = Vector{Float64}()                    # Index of timesteps 
        ews = Vector{Float64}()                         # Statistical distribution over the ensemble

        # Import the AMOC timeseries
        data = NCDataset("../../res/data/amoc/AMOC_transport_depth_0-1000m_monthly.nc")
        time = data["time"][:]
        observable = data["Transport"][:]

        # Convert the analysis of the non-autonomous drift into an ensemble problem
        ensemble = preprocess_solution(time, observable, window_size)
        tipping = ensemble.tipping_point
        Ne = length(ensemble.trajectories)
        Nt = length(ensemble.trajectories[1])

        # Extract windowed subseries
        tw = copy(ensemble.timesteps[1])
        uw = copy(ensemble.trajectories[1]) 

        # Detrend it 
        PyEMD = pyimport("PyEMD")
        emd = PyEMD.EMD()
        imfs = Array(emd(uw))
        trend = sum(imfs[(end-3):end,:], dims=1)[:]
        residual = uw .- trend

        # IMFS figure
        fig = Figure(;size = (2000, 2000))
        ax1 = Axis(fig[1,1], title="IMFs")
        colors = cgrad(:Reds, size(imfs)[1])
        for n in 1:size(imfs)[1]
                lines!(ax1, tw, imfs[n,:], color = colors[n], linewidth = 2.0)
        end

        ax2 = Axis(fig[2,1], title="IMF[11] + IMF[10]")
        lines!(ax2, tw, imfs[end,:], color = :red, linewidth = 2.0)
        lines!(ax2, tw, sum(imfs[(end-1):end,:], dims=1)[:], color = :black, linewidth = 2.0)

        ax3 = Axis(fig[3,1], title="IMF[11] + IMF[10] + IMF[9]")
        lines!(ax3, tw, imfs[end,:], color = :red, linewidth = 2.0)
        lines!(ax3, tw, sum(imfs[(end-2):end,:], dims=1)[:], color = :black, linewidth = 2.0)

        ax4 = Axis(fig[4,1], title="IMF[11] + IMF[10] + IMF[9] + IMF[8]")
        lines!(ax4, tw, imfs[end,:], color = :red, linewidth = 2.0)
        lines!(ax4, tw, sum(imfs[(end-3):end,:], dims=1)[:], color = :black, linewidth = 2.0)

        ax5 = Axis(fig[5,1], title="IMF[11] + IMF[10] + IMF[9] + IMF[8] + IMF[7]")
        lines!(ax5, tw, imfs[end,:], color = :red, linewidth = 2.0)
        lines!(ax5, tw, sum(imfs[(end-4):end,:], dims=1)[:], color = :black, linewidth = 2.0)

        savefig("amoc/residuals/imfs.png", fig)

        # Detrending figure
        fig = Figure()

        ax1 = Axis(fig[1,1:2], title="Original signal w/ trend")
        lines!(ax1, tw, uw, color = :black, linewidth = 1.0)
        lines!(ax1, tw, trend, color = :red, linewidth = 3.0)

        ax2 = Axis(fig[1,3], title="Empirical distribution")
        Nb = convert(Int64, ceil(abs(maximum(uw)-minimum(uw))/(3.49*std(uw)*(length(uw))^(-1.0/3.0))))
        bins, pdf = fit_distribution(uw, n_bins=Nb+1)
        barplot!(ax2, bins, pdf, color = pdf, colormap = [:white,CtpGray], strokecolor = :black, strokewidth = 1)

        ax3 = Axis(fig[2,1:2], title="Detrended signal")
        lines!(ax3, tw, residual, color = :black, linewidth = 1.0)

        ax4 = Axis(fig[2,3], title="Empirical distribution")
        Nb = convert(Int64, ceil(abs(maximum(residual)-minimum(residual))/(3.49*std(residual)*(length(residual))^(-1.0/3.0))))
        bins, pdf = fit_distribution(residual, n_bins=Nb+1)
        barplot!(ax4, bins, pdf, color = pdf, colormap = [:white,CtpGray], strokecolor = :black, strokewidth = 1)

        ax5 = Axis(fig[3,1:2], title="Smallest scale")
        lines!(ax5, tw, imfs[1,:], color = :black, linewidth = 1.0)

        ax6 = Axis(fig[3,3], title="Empirical distribution")
        Nb = convert(Int64, ceil(abs(maximum(imfs[1,:])-minimum(imfs[1,:]))/(3.49*std(imfs[1,:])*(length(imfs[1,:]))^(-1.0/3.0))))
        bins, pdf = fit_distribution(imfs[1,:], n_bins=Nb+1)
        barplot!(ax6, bins, pdf, color = pdf, colormap = [:white,CtpGray], strokecolor = :black, strokewidth = 1)

        savefig("amoc/residuals/detrending.png", fig)
end

# Execute the main
main()
