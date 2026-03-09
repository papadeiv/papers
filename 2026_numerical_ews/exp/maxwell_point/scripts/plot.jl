"""
    Plotting script

Functions used to create the plots in each figure.
"""

# Create empty layouts for the figures
include("./figs.jl")

# Plot and export the empirical distribution and reconstructed potential
function plot(bins, pdf, coeff, idx)
        # Plot the empirical distribution
        barplot!(ax1, bins, pdf, color = pdf, colormap = [CtpWhite,CtpGray], strokecolor = :black, strokewidth = 2)

        # Plot the ground truth potential and equilibrium distribution
        domain = LinRange(-2.5,2.5,1000)
        lines!(ax1, domain, [ρ(x, μ) for x in domain], color = :red, linewidth = 5.0)
        lines!(ax2, domain, [U(x, μ) for x in domain], color = :red, linewidth = 5.0)

        # Plot the arbitrary (reconstructed) potential and its equilibrium distribution
        lines!(ax1, domain, [p(x, coeff) for x in domain], color = CtpBlue, linewidth = 5.0)
        lines!(ax2, domain, [c0 + V(x, coeff) for x in domain], color = CtpBlue, linewidth = 5.0)

        # Export the figure
        savefig("maxwell/$idx.png", fig1)
end
