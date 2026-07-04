"""
    Plotting script

Functions used to create the plots in each figure.
"""

# Create empty layouts for the figures
include("./figs.jl")

# Plot the full solution, the windowed drift of the stationary distribution and its potential
function plot_solution(parameter, solution, tip_idx, sub_parameter, sub_solution, Vs, idx)
        # Plot the full solution
        ax4.limits = ((parameter[1],parameter[end]),(1.1*minimum(solution),1.1*maximum(solution)))
        ax4.xticks = [parameter[1], sub_parameter[end], parameter[end]]
        ax4.yticks = [1.1*minimum(solution),1.1*maximum(solution)]
        lines!(ax4, parameter[1:tip_idx], solution[1:tip_idx], color = (CtpRed,0.5), linewidth = 3.0)
        lines!(ax4, parameter[tip_idx:end], solution[tip_idx:end], color = (CtpGray,0.25), linewidth = 3.0)

        # Plot the subseries within the window
        lines!(ax4, parameter[idx:(length(sub_parameter) + idx)], solution[idx:(length(sub_solution) + idx)], color = CtpPeach, linewidth = 3.0)

        # Plot the sliding window
        poly!(ax4, Point2f[(parameter[idx], 1.1*minimum(solution)), 
                           (parameter[length(sub_parameter) + idx], 1.1*minimum(solution)),
                           (parameter[length(sub_parameter) + idx], 1.1*maximum(solution)),
                           (parameter[idx], 1.1*maximum(solution)),
                          ], color = (CtpGray, 0.1), strokecolor = :black, strokewidth = 1.0)

        # Plot the (empirical) drifting distribution 
        n_bins = convert(Int64, ceil(abs(maximum(sub_solution)-minimum(sub_solution))/(3.49*std(sub_solution)*(length(sub_solution))^(-1.0/3.0))))
        bins, pdf = fit_distribution(sub_solution, n_bins = n_bins)
        barplot!(ax5, bins, pdf, color = pdf, colormap = [(CtpYellow,0.5),(CtpPeach,1.00)], strokecolor = :black, strokewidth = 1)
        ax5.limits = ((-0.3,0.3),(-0.5,16))
        ax5.xticks = [-0.3, 0, 0.3]

        # Plot the drifting potential
        domain = LinRange(-1.5, 1.5, 1000)
        for μ in sub_parameter
                lines!(ax6, domain, [U(x, μ) for x in domain], color = (CtpPeach,0.01), linewidth = 0.5)
        end

        # Plot the least-square approximation of the drifting potential
        lines!(ax6, domain, [Vs(x) for x in domain], color = CtpBlue, linewidth = 6.0)

        # Export the figure
        savefig("ramped/$idx.png", fig4)
end

# Plot the full solution, escape early-warning and numerical error 
function plot_ews(parameter, solution, timeseries, tip_idx)
        # Extract the three different timeseries
        μ = timeseries[:,1]
        ews = timeseries[:,2]
        error = timeseries[:,3]

        # Plot the three timeseries
        #lines!(ax1, parameter[1:tip_idx], solution[1:tip_idx], color = CtpRed, linewidth = 3.0)
        #lines!(ax1, parameter[tip_idx:end], solution[tip_idx:end], color = CtpGray, linewidth = 3.0)
        lines!(ax2, μ, ews, color = CtpPeach, linewidth = 3.0)
        lines!(ax3, μ, error, color = CtpPeach, linewidth = 3.0)
        ax2.limits = ((μ[1],μ[end]),(0.9*minimum(ews),1.1*maximum(ews)))
        ax3.limits = ((μ[1],μ[end]),(0.9*minimum(error),1.1*maximum(error)))

        #=
        # Customise ticks and limits of the plots
        ax1.limits = ((parameter[1],parameter[end]),(1.1*minimum(solution),1.1*maximum(solution)))
        ax1.xticks = [parameter[1], parameter[end]]
        ax1.yticks = [1.1*minimum(solution),1.1*maximum(solution)]

        ax2.limits = ((parameter[1],parameter[end]),(0.9*minimum(ews),1.1*maximum(ews)))
        ax2.xticks = [parameter[1], parameter[end]]
        ax2.yticks = [0.9*minimum(ews),1.1*maximum(ews)]

        ax3.limits = ((parameter[1],parameter[end]),(0.9*minimum(error),1.1*maximum(error)))
        ax3.xticks = [parameter[1], parameter[end]]
        ax3.yticks = [0.9*minimum(error),1.1*maximum(error)]
        =#

        savefig("ramped/ews.png", fig1)
end
