"""
    Plotting script
"""

# Create empty layouts for the figures
include("./figs.jl")

# Array of colors for different trajectories
colors = [CtpRed, CtpBlue, CtpPeach, CtpMauve, CtpGreen]

# Plot the 'time-paths' of the map
function plot_paths(U)
        # Extract the state variables
        gt, At, et, Lt = columns(U)

        # Define the time domain
        t = collect(UnitRange(0, N))

        # Plot the state variables' time-paths
        scatterlines!(ax1, t, gt, color = colors[plt_idx], linewidth = 2.0, markersize = 12)
        scatterlines!(ax2, t, At, color = colors[plt_idx], linewidth = 2.0, markersize = 12)
        scatterlines!(ax3, t, et, color = colors[plt_idx], linewidth = 2.0, markersize = 12)
        scatterlines!(ax4, t, Lt, color = colors[plt_idx], linewidth = 2.0, markersize = 12)

        #=
        # Specify the limits of the plots
        ax1.limits = ((t[1], t[end]), (minimum(gt) - 0.05*(maximum(gt) - minimum(gt)), maximum(gt) + 0.05*(maximum(gt) - minimum(gt))))
        ax3.limits = ((t[1], t[end]), (minimum(et) - 0.05*(maximum(et) - minimum(et)), maximum(et) + 0.05*(maximum(et) - minimum(et))))
        ax4.limits = ((t[1], t[end]), (minimum(Lt) - 0.05*(maximum(Lt) - minimum(Lt)), maximum(Lt) + 0.05*(maximum(Lt) - minimum(Lt))))

        # Specify the ticks of the plots 
        ax1.yticks = [minimum(gt), 0.8, maximum(gt)]
        ax3.yticks = [minimum(et), 0.025, maximum(et)]
        ax4.yticks = [minimum(Lt), 6, maximum(Lt)]
        =#

        # Update global plotting index
        global plt_idx = plt_idx + 1
end
