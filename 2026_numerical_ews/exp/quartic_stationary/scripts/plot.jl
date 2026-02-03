"""
    Plotting script

Functions used to create the plots in each figure.
"""

# Create empty layouts for the figures
include("./figs.jl")

# Plot the timeseries and the solution of the nonlinear reconstruction of the potential
function plot_solutions(timestamps, timeseries, shifted_potential, idx)
        # Plot the ground truth 
        domain = LinRange(-1.5,1.5,1000)
        lines!(ax2, domain, [U(x, μ[idx]) for x in domain], color = (:red,1.0), linewidth = 3.0)

        # Plot the timeseries
        #lines!(ax1, timestamps, timeseries, color = (:black, 0.15))
        
        # Plot the reconstructed shifted potential
        lines!(ax2, domain, [shifted_potential(x) for x in domain], color = (:black,0.25), linewidth = 2.0)

        # Setup ticks and limits for the plot
        ax1.limits = ((timestamps[1], timestamps[end]), (1.1*minimum(timeseries), 1.1*maximum(timeseries)))
end

# Plot the escape rate early-warning signal 
function plot_ews()
        # Plot the ews of the ground truth
        μ_range = LinRange(-0.05, 1.35, 1000)
        Xs = [(get_equilibria(f, μ, domain=[-10,10]).stable)[2] for μ in μ_range]
        Xu = [(get_equilibria(f, μ, domain=[-10,10]).unstable)[1] for μ in μ_range]
        lines!(ax3, μ_range, [exp(-(U(xu, μ) - U(xs, μ))) for (μ, xs, xu) in zip(μ_range, Xs, Xu)], color = CtpRed, linewidth = 5.0)

        # Loop over the parameter values
        for v in μ
                # Import the random variables
                df = CSV.read("../../res/data/weighted/μ=$v.csv", DataFrame)

                # Extract the ones under analysis
                ews = df.LDP

                # Filter out the outliers and plot the ews 
                lower_cutoff = quantile(filter(isfinite, ews), 0.05)
                upper_cutoff = quantile(filter(isfinite, ews), 0.95)
                filtered_ews = filter(x -> lower_cutoff ≤ x ≤ upper_cutoff, ews)
                errorbars!(ax3, [v], [median(filtered_ews)], [median(filtered_ews) - quantile(filtered_ews,0.25)], [quantile(filtered_ews,0.75) - median(filtered_ews)], color = :black, whiskerwidth = 15, linewidth = 3.0)
                scatter!(ax3, v, median(filtered_ews), color = CtpBlue, markersize = 25, strokewidth = 2.0)

                # Print the percentile distribution
                println("μ = ", v)
                println("Lower cutoff = ", lower_cutoff, ". Upper cutoff = ", upper_cutoff)
                println("Q1 = ", quantile(filtered_ews,0.25), ", Q2 = ", median(filtered_ews), ", Q3 = ", quantile(filtered_ews,0.75))
                println("------------------------------------------------------------------")
        end

        # Export the figure
        savefig("ews.png", fig3)
end
