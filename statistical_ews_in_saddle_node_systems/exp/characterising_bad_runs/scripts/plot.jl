"""
    Plotting script

Collection of all the functions used to generate and properly format the figures of the simulations.
"""

# Create empty layouts for the figures
include("./figs.jl")

# Plot the timeseries and the solution of the nonlinear reconstruction of the potential
function plot_solutions(timesteps, timeseries, shifted_potential, idx)
        # Plot the ground truth 
        domain = LinRange(-1.5,1.5,1000)
        lines!(ax2, domain, [U(x,μ) for x in domain], color = (:red,1.0), linewidth = 3.0)

        # Extract the length of the timeseries
        Nt = length(timeseries)

        # Plot the reconstructed shifted potential
        lines!(ax2, domain, [shifted_potential(x) for x in domain], color = (:black,0.5), linewidth = 2.0)

        # Export the figure
        savefig("dt=$dt/μ=$μ/$idx.png", fig1)
end
