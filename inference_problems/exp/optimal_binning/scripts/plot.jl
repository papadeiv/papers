"""
    Plotting script

Collection of all the functions used to generate and properly format the figures of the simulations.
"""

# Create empty layouts for the figures
include("./figs.jl")

# Plot the histograms
function plot_hist(samples, bins, pdf, idx)
        # Check whether the plot is for the OUP or Gaussian samples
        if idx == 1
                # Plot the OUP samples
                scatter!(ax3, samples, zeros(length(samples)), color = (CtpRed,0.01), markersize = 30)
                # Plot binned histogram
                barplot!(ax1, bins, pdf, color = pdf, colormap = [(:white,1.0),(CtpRed,1.0)], strokecolor = :black, strokewidth = 1)
                # Plot the stationary density
                domain = collect(LinRange(-0.4,0.4,1000))
                local G(x) = sqrt(θ/(pi*σ^2))*exp(-(θ/(σ^2))*(x - μ)^2) 
                lines!(ax1, domain, [G(x) for x in domain], color = (CtpGreen,0.5), linewidth = 5.0)

                # Plot the misfit
                barplot!(ax5, bins, pdf .- [G(x) for x in bins], color = pdf, colormap = [(:white,1.0),(CtpRed,1.0)], strokecolor = :black, strokewidth = 1)
        else
                # Plot the Gaussian samples
                scatter!(ax4, samples, zeros(length(samples)), color = (CtpBlue,0.01), markersize = 35)
                # Plot binned histogram
                barplot!(ax2, bins, pdf, color = pdf, colormap = [(:white,1.0),(CtpBlue,1.0)], strokecolor = :black, strokewidth = 1)
                # Plot the normal distribution
                domain = collect(LinRange(-0.4,0.4,1000))
                σ_g = sqrt((σ^2)/(2*θ))
                local H(x) = sqrt(1.0/(2*pi*σ_g^2))*exp(-(1.0/(2*σ_g^2))*(x - μ)^2)
                lines!(ax2, domain, [H(x) for x in domain], color = (CtpGreen,0.5), linewidth = 5.0)

                # Plot the misfit
                barplot!(ax6, bins, pdf .- [H(x) for x in bins], color = pdf, colormap = [(:white,1.0),(CtpBlue,1.0)], strokecolor = :black, strokewidth = 1)
        end
end
