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
        lines!(ax10, [t[1], t[end]], (1/(1-γ)).*ones(2), color = (:black, 0.5), linestyle = :dash, linewidth = 5.0)
end
