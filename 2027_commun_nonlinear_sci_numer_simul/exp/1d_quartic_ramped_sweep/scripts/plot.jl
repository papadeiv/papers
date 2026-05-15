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
        lines!(ax1, domain, [U(x, μf) for x in domain], color = :red, linewidth = 5.0)

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
end

# Plot the escape rate early-warning signal 
function plot_ews()
        # Plot the ews of the ground truth
        μ_range = LinRange(-0.05, 1.35, 1000)
        Xs = [(get_equilibria(f, μ, domain=[-10,10]).stable)[2] for μ in μ_range]
        Xu = [(get_equilibria(f, μ, domain=[-10,10]).unstable)[1] for μ in μ_range]
        lines!(ax3, μ_range, [exp(-(U(xu, μ) - U(xs, μ))) for (μ, xs, xu) in zip(μ_range, Xs, Xu)], color = CtpRed, linewidth = 5.0)

        # Loop over the parameter values
        for n in 1:length(M0)
                # Import the random variables 
                data = import_data(n)
                coefficients = data.coefficients
                analysis = data.analysis

                # Extract the ones under analysis
                ews = analysis[:,2]

                # Filter out the outliers and plot the ews 
                lower_cutoff = quantile(filter(isfinite, ews), 0.05)
                upper_cutoff = quantile(filter(isfinite, ews), 0.95)
                filtered_ews = filter(x -> lower_cutoff ≤ x ≤ upper_cutoff, ews)
                errorbars!(ax3, [M0[n]], [median(filtered_ews)], [median(filtered_ews) - quantile(filtered_ews,0.25)], [quantile(filtered_ews,0.75) - median(filtered_ews)], color = :black, whiskerwidth = 15, linewidth = 3.0)
                scatter!(ax3, M0[n], median(filtered_ews), color = CtpBlue, markersize = 25, strokewidth = 2.0)

                # Print the percentile distribution
                println("μ = ", M0[n])
                println("Lower cutoff = ", lower_cutoff, ". Upper cutoff = ", upper_cutoff)
                println("Q1 = ", quantile(filtered_ews,0.25), ", Q2 = ", median(filtered_ews), ", Q3 = ", quantile(filtered_ews,0.75))
                println("------------------------------------------------------------------")
        end

        # Export the figure
        savefig("ramped_ensemble/ews.png", fig3)
end
