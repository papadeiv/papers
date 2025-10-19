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

        # Plot the timeseries
        lines!(ax1, timesteps[1:Nt], timeseries[1:Nt], color = (:black,0.15), linewidth = 3.0)
        
        # Plot the reconstructed shifted potential
        lines!(ax2, domain, [shifted_potential(x) for x in domain], color = (:black,0.5), linewidth = 2.0)

        # Setup ticks and limits for the plot
        y_range = maximum(timeseries) - minimum(timeseries)
        ax1.limits = ((timesteps[1], timesteps[end]), (minimum(timeseries) - 0.1*y_range, maximum(timeseries) + 0.1*y_range))
        ax1.xticks = [timesteps[1], timesteps[end]]
        ax1.yticks = [minimum(timeseries), maximum(timeseries)]

        # Export the figure
        savefig("dt=$dt/μ=$μ/$idx.png", fig1)
end

# Plot the escape rate early-warning signal 
function plot_ews()
        # Plot the ews of the ground truth
        lines!(ax3, LinRange(-1.05,-0.05,2000), [exp(-(U(-sqrt(-μ),μ) - U(sqrt(-μ),μ))) for μ in LinRange(-1.05,-0.05,2000)], color = CtpRed, linewidth = 5.0)
        lines!(ax4, LinRange(-1.05,-0.05,2000), [sqrt(abs(Uxx(-sqrt(-μ),μ))*Uxx(sqrt(-μ),μ)) for μ in LinRange(-1.05,-0.05,2000)], color = CtpRed, linewidth = 5.0)
        lines!(ax5, LinRange(-1.05,-0.05,2000), [sqrt(abs(Uxx(-sqrt(-μ),μ))*Uxx(sqrt(-μ),μ))*exp(-(U(-sqrt(-μ),μ) - U(sqrt(-μ),μ))) for μ in LinRange(-1.05,-0.05,2000)], color = CtpRed, linewidth = 5.0)

        # Discrete parameter values
        μ_range = [-1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1]

        # Loop over the parameter values
        for μ in μ_range
                # Import the random variables
                df = CSV.read("../../res/data/dt=$dt/μ=$μ.csv", DataFrame)

                # Extract the ones under analysis
                LDP = df.LDP
                C = df.C
                ews = C.*LDP
                N_escapes = df.Nx

                # Filter out the outliers and plot the LDP
                lower_cutoff = quantile(filter(isfinite, LDP), 0.05)
                upper_cutoff = quantile(filter(isfinite, LDP), 0.95)
                filtered_LDP = filter(x -> lower_cutoff ≤ x ≤ upper_cutoff, LDP)
                errorbars!(ax3, [μ], [median(filtered_LDP)], [median(filtered_LDP) - quantile(filtered_LDP,0.25)], [quantile(filtered_LDP,0.75) - median(filtered_LDP)], color = :black, whiskerwidth = 15, linewidth = 3.0)
                scatter!(ax3, μ, median(filtered_LDP), color = CtpBlue, markersize = 25, strokewidth = 2.0)

                # Plot the empirical escape rate (number of escapes)
                scatter!(ax3, μ, N_escapes[1], color = CtpGreen, markersize = 25, marker = :star5, strokewidth = 2.0)

                # Print the percentile distribution
                println("μ = ", μ)
                println("Lower cutoff = ", lower_cutoff, ". Upper cutoff = ", upper_cutoff)
                println("Q1 = ", quantile(filtered_LDP,0.25), ", Q2 = ", median(filtered_LDP), ", Q3 = ", quantile(filtered_LDP,0.75))
                println("------------------------------------------------------------------")

                # Filter out the outliers and plot the prefactor 
                lower_cutoff = quantile(filter(isfinite, C), 0.05)
                upper_cutoff = quantile(filter(isfinite, C), 0.95)
                filtered_C = filter(x -> lower_cutoff ≤ x ≤ upper_cutoff, C)
                errorbars!(ax4, [μ], [median(filtered_C)], [median(filtered_C) - quantile(filtered_C,0.25)], [quantile(filtered_C,0.75) - median(filtered_C)], color = :black, whiskerwidth = 15, linewidth = 3.0)
                scatter!(ax4, μ, median(filtered_C), color = CtpBlue, markersize = 25, strokewidth = 2.0)

                # Filter out the outliers and plot the ews 
                lower_cutoff = quantile(filter(isfinite, ews), 0.05)
                upper_cutoff = quantile(filter(isfinite, ews), 0.95)
                filtered_ews = filter(x -> lower_cutoff ≤ x ≤ upper_cutoff, ews)
                errorbars!(ax5, [μ], [median(filtered_ews)], [median(filtered_ews) - quantile(filtered_ews,0.25)], [quantile(filtered_ews,0.75) - median(filtered_ews)], color = :black, whiskerwidth = 15, linewidth = 3.0)
                scatter!(ax5, μ, median(filtered_ews), color = CtpBlue, markersize = 25, strokewidth = 2.0)
        end

        # Export the figure
        savefig("dt=$dt/LDP.png", fig3)
        savefig("dt=$dt/prefactor.png", fig4)
        savefig("dt=$dt/ews.png", fig5)
end
