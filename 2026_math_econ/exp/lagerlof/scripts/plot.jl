"""
    Plotting script
"""

# Create empty layouts for the figures
include("./figs.jl")

# Plot the 'time-paths' of the map
function plot_paths(U)
        # Extract the state variables
        gt, At, et, Lt = columns(U)

        # Define the time domain
        t = collect(UnitRange(0, N))

        # Plot the state variables' time-paths
        lines!(ax1, t, gt, color = CtpRed, linewidth = 5.0)
        lines!(ax2, t, At, color = CtpRed, linewidth = 5.0)
        lines!(ax3, t, et, color = CtpRed, linewidth = 5.0)
        lines!(ax4, t, Lt, color = CtpRed, linewidth = 5.0)

        # Plot boxes for detailed subplots
        poly!(ax1, Point2f[(0, 0), (50, 0), (50, 0.8), (0, 0.8)], color = :transparent, strokecolor = :black, strokewidth = 5.0)
        poly!(ax3, Point2f[(0, minimum(et) - 0.025*(maximum(et) - minimum(et))), (50, minimum(et) - 0.025*(maximum(et) - minimum(et))), (50, 0.025), (0, 0.025)], color = :transparent, strokecolor = :black, strokewidth = 5.0)
        poly!(ax4, Point2f[(0, minimum(Lt) - 0.025*(maximum(Lt) - minimum(Lt))), (50, minimum(Lt) - 0.025*(maximum(Lt) - minimum(Lt))), (50, 6), (0, 6)], color = :transparent, strokecolor = :black, strokewidth = 5.0)

        # Specify the limits of the plots
        ax1.limits = ((t[1], t[end]), (minimum(gt) - 0.05*(maximum(gt) - minimum(gt)), maximum(gt) + 0.05*(maximum(gt) - minimum(gt))))
        ax3.limits = ((t[1], t[end]), (minimum(et) - 0.05*(maximum(et) - minimum(et)), maximum(et) + 0.05*(maximum(et) - minimum(et))))
        ax4.limits = ((t[1], t[end]), (minimum(Lt) - 0.05*(maximum(Lt) - minimum(Lt)), maximum(Lt) + 0.05*(maximum(Lt) - minimum(Lt))))

        # Specify the ticks of the plots 
        ax1.yticks = [minimum(gt), 0.8, maximum(gt)]
        ax3.yticks = [minimum(et), 0.025, maximum(et)]
        ax4.yticks = [minimum(Lt), 6, maximum(Lt)]

        # Plot the detailed subplots
        lines!(ax5, t, gt, color = CtpRed, linewidth = 5.0)
        lines!(ax6, t, et, color = CtpRed, linewidth = 5.0)
        lines!(ax7, t, Lt, color = CtpRed, linewidth = 5.0)
        
        # Specify the limits of the plots
        ax5.limits = ((t[1], 50), (0, 0.8))
        ax6.limits = ((t[1], 50), (minimum(et) - 0.025*(maximum(et) - minimum(et)), 0.025))
        ax7.limits = ((t[1], 50), (minimum(Lt) - 0.025*(maximum(Lt) - minimum(Lt)), 6))

        # Specify the ticks of the plots 
        ax5.yticks = [minimum(gt), 0.8]
        ax6.yticks = [minimum(et), 0.025]
        ax7.yticks = [minimum(Lt), 6]

        # Export the figure
        savefig("timepaths.png", fig)
end
