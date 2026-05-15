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
        scatterlines!(ax1, t, gt, color = CtpRed, linewidth = 2.0, markersize = 12)
        scatterlines!(ax2, t, At, color = CtpRed, linewidth = 2.0, markersize = 12)
        scatterlines!(ax3, t, et, color = CtpRed, linewidth = 2.0, markersize = 12)
        scatterlines!(ax4, t, Lt, color = CtpRed, linewidth = 2.0, markersize = 12)

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
        scatterlines!(ax5, t, gt, color = CtpRed, linewidth = 2.0, markersize = 12)
        scatterlines!(ax6, t, et, color = CtpRed, linewidth = 2.0, markersize = 12)
        scatterlines!(ax7, t, Lt, color = CtpRed, linewidth = 2.0, markersize = 12)
        
        # Specify the limits of the plots
        ax5.limits = ((t[1], 50), (0, 0.8))
        ax6.limits = ((t[1], 50), (minimum(et) - 0.025*(maximum(et) - minimum(et)), 0.025))
        ax7.limits = ((t[1], 50), (minimum(Lt) - 0.025*(maximum(Lt) - minimum(Lt)), 6))

        # Specify the ticks of the plots 
        ax5.yticks = [minimum(gt), 0.8]
        ax6.yticks = [minimum(et), 0.025]
        ax7.yticks = [minimum(Lt), 6]
end

# Plot the observables of the state variables 
function plot_observables(U)
        # Extract the state variables
        gt, At, et, Lt = columns(U)

        # Compute the observables
        ht = [h(e, g) for (e, g) in zip(et, gt)]
        xt = [x(A, L) for (A, L) in zip(At, Lt)]
        zt = [z(g, A, e, L) for (g, A, e, L) in zip(gt, At, et, Lt)]

        # Define the time domain
        t = collect(UnitRange(0, N))

        # Plot the state variables' time-paths
        scatterlines!(ax8, t, ht, color = CtpRed, linewidth = 2.0, markersize = 12)
        scatterlines!(ax9, t, xt, color = CtpRed, linewidth = 2.0, markersize = 12)
        scatterlines!(ax10, t, zt, color = CtpRed, linewidth = 2.0, markersize = 12)

        # Plot the treshold for the switching dynamics of Lt
        lines!(ax10, [t[1], t[end]], (1/(1-γ)).*ones(2), color = (:black,0.5), linestyle = :dash, linewidth = 5.0)
end
