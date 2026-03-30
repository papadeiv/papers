"""
    Plotting script

Functions used to create the plots in each figure.
"""

# Create empty layouts for the figures
include("./figs.jl")

# Plot the trajectory of the particle inside the double-well potential
function plot_trajectory(x_component, y_component, equilibria, parameter_idx)
        # Define the domain of the potential
        X = collect(LinRange(-2.5, 2.5, 1000))
        Y = collect(LinRange(-2.5, 2.5, 1000))

        # Create figure and axis for the contour
        fig = Figure()
        ax = Axis(fig[1,1], limits = ((-2.5, 2.5), (-2.5, 2.5)))

        # Create figure and axis for the timeseries
        Fig = Figure()
        ax1 = Axis(Fig[1,1])
        ax2 = Axis(Fig[2,1])
        ax3 = Axis(Fig[3,1])

        # Plot the contour of the potential
        U(x, y) = x^4 - 2*x^2 + μ[parameter_idx]*x + y^2 + 1 
        W(x, y) = U(a11*x + a12*y, a21*x + a22*y) 
        contourf!(ax, X, Y, [W(x, y) for x in X, y in Y], levels = 200)

        # Plot the projecting subspaces and the projected timeseries
        α = compute_subspace(μ[parameter_idx])
        lines!(ax, X, [α*x for x in X], color = CtpPink, linewidth = 3.0)
        v = [1, α]
        u = [dot([x_component[n], y_component[n]], v) for n in 1:length(x_component)]
        display(var(u))
        lines!(ax2, LinRange(1, length(u), length(u)), u, color = CtpPink)
        lines!(ax, X, zeros(length(X)), color = CtpTeal, linewidth = 3.0)
        lines!(ax3, LinRange(1, length(u), length(u)), x_component, color = CtpTeal)
        display(var(x_component))

        # Plot the optimal subspace and the projected timeseries 
        α = compute_optimal(hcat(x_component, y_component))
        lines!(ax, X, [α*x for x in X], color = :red, linestyle = :dash, linewidth = 3.0)
        v = [1, α]
        u = [dot([x_component[n], y_component[n]], v) for n in 1:length(x_component)]
        display(var(u))
        println("")
        lines!(ax1, LinRange(1, length(u), length(u)), u, color = :red)

        # Plot the trajectory
        lines!(ax, x_component, y_component, color = (:white, 0.75), linewidth = 0.125)

        # Plot the equilibria
        scatter!(ax, [eq[1] for eq in equilibria.stable], [eq[2] for eq in equilibria.stable], markersize = 15.0, color = :green1, strokewidth = 1.0, strokecolor = :black)
        scatter!(ax, [eq[1] for eq in equilibria.unstable], [eq[2] for eq in equilibria.unstable], markersize = 25.0, marker = :star5, color = :yellow, strokewidth = 1.0, strokecolor = :black)

        # Export the figure
        savefig("2d_potential/μ=$(μ[parameter_idx]).png", fig)
        savefig("projections/μ=$(μ[parameter_idx]).png", Fig)
end

# Plot the timeseries and the solution of the nonlinear reconstruction of the potential
function plot_solutions(timestamps, timeseries, shifted_potential, idx)
        # Plot the ground truth 
        domain = LinRange(-1.5,1.5,1000)
        lines!(ax2, domain, [U0(x, μ[idx]) for x in domain], color = (:red,1.0), linewidth = 3.0)

        # Plot the timeseries
        lines!(ax1, timestamps, timeseries, color = (:black, 0.15))
        
        # Plot the reconstructed shifted potential
        lines!(ax2, domain, [shifted_potential(x) for x in domain], color = (:black,0.25), linewidth = 2.0)

        # Setup ticks and limits for the plot
        ax1.limits = ((timestamps[1], timestamps[end]), (1.1*minimum(timeseries), 1.1*maximum(timeseries)))
end

# Plot the escape rate early-warning signal 
function plot_ews()
        # Plot the ews of the ground truth
        μ_range = LinRange(-0.5, 1.35, 1000)
        f0(x, μ) = f1(x, 0, μ)
        Xs = [(get_equilibria(f0, μ, domain=[-10,10]).stable)[2] for μ in μ_range]
        Xu = [(get_equilibria(f0, μ, domain=[-10,10]).unstable)[1] for μ in μ_range]
        lines!(ax3, μ_range, [exp(-(U0(xu, μ) - U0(xs, μ))) for (μ, xs, xu) in zip(μ_range, Xs, Xu)], color = CtpRed, linewidth = 5.0)

        # Loop over the parameter values
        m = 1::Integer
        for v in μ
                # Import the ews 
                ews = readin("ews/$m.csv")
                m = m + 1

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
