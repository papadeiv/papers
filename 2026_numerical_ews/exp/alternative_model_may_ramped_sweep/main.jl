"""
    Main script

Run this file to execute the simulation, analyse and plot the results.
"""

# Import the necessary packages and local modules
include("inc.jl")

# Import the simulation's scripts
include("./scripts/sim.jl")
include("./scripts/proc.jl")
include("./scripts/plot.jl")

# Define the main algorithm
function main()
        # Loop over the initial parameter values
        for (parameter_index, μ0) ∈ enumerate(M0)
                # Compute the initial condition
                equilibria = get_equilibria(f, μ0, domain=[5,10])
                x0 = [equilibria.stable[1], μ0]

                # Solve the ensemble slow-fast SDEs
                ensemble = evolve(f, η, Λ, x0, steps=Nt, stepsize=δt, particles=Ne)
                t = ensemble.time
                μ = ensemble.parameter 

                # Build the critical manifold of the solution 
                critical_manifold = compute_crit_man(x0, μ)

                # Plot the ground truth potential
                domain = LinRange(-2,12,1000)
                fig = Figure()
                ax1 = Axis(fig[1,1:2])
                ax2 = Axis(fig[2,1:2])
                ax3 = Axis(fig[1,3], limits = (domain[1], domain[end], -6, 2))
                lines!(ax3, domain, [U(x, μ0) for x in domain], color = :blue, linewidth = 3.0)

                # Loop over the ensemble's sample paths
                printstyled("[μ0 = $μ0]: solving the likelihood estimation problems\n"; bold=true, underline=true, color=:light_magenta)
                @showprogress for (solution_index, solution) in enumerate(ensemble.state)
                        # Extract and detrend the trajectory 
                        residuals = detrend(solution, qse = critical_manifold).residuals

                        # Compute the traditional early-warning signals
                        (ews[1])[solution_index, parameter_index] = std(residuals)
                        (ews[2])[solution_index, parameter_index] = cov(residuals[1:end-1], residuals[2:end])/var(residuals[1:end-1])

                        # Compute the parameters estimate and the approximate escape rate
                        c = solve_EM_MLE(residuals)
                        V_EM = shift_potential(x0[1], c, μ[1])
                        (ews[3])[solution_index, parameter_index] = compute_ews(c)
                        
                        # Plot the results
                        lines!(ax1, μ[1:truncation_index], solution[1:truncation_index], color = (:black,0.10), linewidth = 1.0)
                        lines!(ax2, μ[1:truncation_index], residuals[1:truncation_index], color = (:blue,0.10), linewidth = 1.0)
                        lines!(ax3, domain, [V_EM(x) for x in domain], color = (:black,0.10), linewidth = 1.0)
                end
                # Plot the critical manifold
                lines!(ax1, μ[1:truncation_index], critical_manifold[1:truncation_index], color = :blue, linewidth = 3.0)
                savefig("$parameter_index.png", fig)
        end

        # Plot the standard deviation ews
        fig = Figure()
        ax = Axis(fig[1,1], title = L"\textbf{standard-deviation}")
        for (ews_index, dist) in enumerate(eachcol(ews[1]))
                # Filter out the outliers and plot the ews 
                lower_cutoff = quantile(filter(isfinite, dist), 0.05)
                upper_cutoff = quantile(filter(isfinite, dist), 0.95)
                filtered_ews = filter(x -> lower_cutoff ≤ x ≤ upper_cutoff, dist)

                # Plot the interquartile range 
                errorbars!(ax, [M0[ews_index]], [median(filtered_ews)], [median(filtered_ews) - quantile(filtered_ews,0.25)], [quantile(filtered_ews,0.75) - median(filtered_ews)], color = :black, whiskerwidth = 15, linewidth = 3.0)

                # Plot the median
                scatter!(ax, M0[ews_index], median(filtered_ews), color = CtpBlue, markersize = 25, strokewidth = 2.0)
        end
        savefig("ews_std.png", fig)

        # Plot the autocorrelation ews
        fig = Figure()
        ax = Axis(fig[1,1], title = L"\textbf{autocorrelation}")
        for (ews_index, dist) in enumerate(eachcol(ews[2]))
                # Filter out the outliers and plot the ews 
                lower_cutoff = quantile(filter(isfinite, dist), 0.05)
                upper_cutoff = quantile(filter(isfinite, dist), 0.95)
                filtered_ews = filter(x -> lower_cutoff ≤ x ≤ upper_cutoff, dist)

                # Plot the interquartile range 
                errorbars!(ax, [M0[ews_index]], [median(filtered_ews)], [median(filtered_ews) - quantile(filtered_ews,0.25)], [quantile(filtered_ews,0.75) - median(filtered_ews)], color = :black, whiskerwidth = 15, linewidth = 3.0)

                # Plot the median
                scatter!(ax, M0[ews_index], median(filtered_ews), color = CtpBlue, markersize = 25, strokewidth = 2.0)
        end
        savefig("ews_ac1.png", fig)

        # Plot the escape rate ews
        fig = Figure()
        ax = Axis(fig[1,1], title = L"\textbf{escape rate}")
        parameter_range = LinRange(M0[1], M0[end], 1000)
        Xs = [(get_equilibria(f, μ, domain=[0,10]).stable)[end] for μ in parameter_range]
        Xu = [(get_equilibria(f, μ, domain=[0,10]).unstable)[end] for μ in parameter_range]
        lines!(ax, parameter_range, [exp(-(U(xu, μ) - U(xs, μ))) for (μ, xs, xu) in zip(parameter_range, Xs, Xu)], color = CtpRed, linewidth = 5.0)
        for (ews_index, dist) in enumerate(eachcol(ews[3]))
                # Filter out the outliers and plot the ews 
                lower_cutoff = quantile(filter(isfinite, dist), 0.05)
                upper_cutoff = quantile(filter(isfinite, dist), 0.95)
                filtered_ews = filter(x -> lower_cutoff ≤ x ≤ upper_cutoff, dist)

                # Plot the interquartile range 
                errorbars!(ax, [M0[ews_index]], [median(filtered_ews)], [median(filtered_ews) - quantile(filtered_ews,0.25)], [quantile(filtered_ews,0.75) - median(filtered_ews)], color = :black, whiskerwidth = 15, linewidth = 3.0)

                # Plot the median
                scatter!(ax, M0[ews_index], median(filtered_ews), color = CtpBlue, markersize = 25, strokewidth = 2.0)
        end
        savefig("ews_esc.png", fig)
end

# Execute the main
main()
