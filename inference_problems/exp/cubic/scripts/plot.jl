"""
    Plotting script

Collection of all the functions used to generate and properly format the figures of the simulations.
"""

# Create empty layouts for the figures
include("./figs.jl")

# Plot the timeseries and the solution of the nonlinear reconstruction of the potential
function plot_solutions(bins, pdf, linear_coeff, nonlinear_coeff, linear_potential, nonlinear_potential, run_label)
        # Plot the histograms
        scatter!(ax1, bins, pdf, color = (CtpRed,0.5), markersize = 20, strokewidth = 2.0)
        scatter!(ax2, bins, pdf, color = (CtpRed,0.5), markersize = 20, strokewidth = 2.0)

        # Plot the ground truth 
        domain = LinRange(-1.5,1.5,1000)
        lines!(ax3, domain, [U(x,μ) for x in domain], color = (CtpRed,1.0), linewidth = 5.0)
        lines!(ax4, domain, [U(x,μ) for x in domain], color = (CtpRed,1.0), linewidth = 5.0)

        # Plot the reconstructed probability density
        domain = LinRange(-0.35,0.35,5000)
        lines!(ax1, domain, [p(x,linear_coeff) for x in domain], color = CtpBlue, linewidth = 5.0)
        lines!(ax2, domain, [p(x,nonlinear_coeff) for x in domain], color = CtpBlue, linewidth = 5.0)
        #display(maximum([p(x,linear_coeff) for x in domain]))
        #display(maximum([p(x,nonlinear_coeff) for x in domain]))

        # Plot the reconstructed shifted potentials
        domain = LinRange(-1.5,1.5,1000)
        lines!(ax3, domain, [linear_potential(x) for x in domain], color = (CtpBlue,1.0), linewidth = 5.0)
        lines!(ax4, domain, [nonlinear_potential(x) for x in domain], color = (CtpBlue,1.0), linewidth = 5.0)

        # Export the figure and the solution
        savefig("cubic/"*run_label*".png", fig1)
end
