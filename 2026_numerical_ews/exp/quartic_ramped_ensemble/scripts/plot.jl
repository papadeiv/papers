"""
    Plotting script

Functions used to create the plots in each figure.
"""

# Create empty layouts for the figures
include("./figs.jl")

# Plot the reconstruction of the shifted potential
function plot_solution(μf)
        # Plot the shifted ground truth
        domain = collect(LinRange(-1.5, 1.5, 1000))
        lines!(ax1, domain, [U(x, μf) for x in domain], color = CtpRed, linewidth = 5.0)

        # Loop over the solutions of the nonlinear least-squares reconstruction
        for solution in solutions
                # Extract the relevant data 
                x0 = solution[1]
                μ = solution[2]
                c = solution[3:5]

                # Compute the shifted potential
                xs, Vs = shift_potential(x0, μ, c)

                # Plot the reconstruction
                lines!(ax1, domain, [Vs(x) for x in domain], color = (:black, 0.25), linewidth = 3.0)
        end

        # Export the figure
        savefig("ramped_ensemble/$glb_idx.png", fig1)
end
