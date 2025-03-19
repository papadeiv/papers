using LaTeXStrings, CairoMakie, Makie.Colors
using CSV, DataFrames
include("ApproximatePotential.jl")
include("EscapeProblem.jl")

# Plot the timeseries of the stochastic process and its histogram
function plot_trajectory(time, states, equilibria, parameter::Float64, noise::Float64)
        # Create the figure 
        CairoMakie.activate!(; px_per_unit = 2)
        fig = Figure(; size = (2000, 800))#, backgroundcolor = :transparent)

        # Create axis for the timeseries 
        ax1 = Axis(fig[1,1:5],
                # Background
                #backgroundcolor = :transparent,
                xgridvisible = false,
                ygridvisible = false,
                limits = ((time[1],time[end]), nothing),
                # Title
                title = L"\mu = %$parameter",
                titlevisible = true,
                titlesize = 22,
                titlealign = :center,
                titlegap = 4.0,
                # Axes labels
                xlabel = "time",
                ylabel = L"x",
                xlabelvisible = false,
                ylabelvisible = false,
                xlabelsize = 20,
                ylabelsize = 20,
                xlabelcolor = :black,
                ylabelcolor = :black,
                xlabelpadding = 0.0,
                ylabelpadding = 0.0,
                # Axes scale, position and direction
                xscale = identity, 
                yscale = identity, #log10,
                xreversed = false,
                yreversed = false,
                xaxisposition = :bottom,
                yaxisposition = :left,
                # Ticks
                xticks = [time[1], time[end]],
                yticks = equilibria[:,1],
                xticksvisible = true,
                yticksvisible = true,
                xticksize = 6,
                yticksize = 6,
                # Ticks labels
                xticklabelsvisible = true,
                yticklabelsvisible = true,
                xticklabelsize = 18,
                yticklabelsize = 18,
                xticklabelalign = (:right, :top),
                yticklabelalign = (:right, :center),
                xticklabelcolor = :black,
                yticklabelcolor = :black,
                xtickformat = "{:.0f}",
                ytickformat = "{:.1f}",
        )
        # Plot the ensemble's individual sample paths
        lines!(ax1, time, states, linewidth = 0.5, color = (:red,0.75))
        # Plot the equilibria of the vector field
        for n in 1:length(equilibria[:,1])
                if equilibria[n,2]==1
                        lines!(ax1, [time[1], time[end]], [equilibria[n,1], equilibria[n,1]], linewidth = 1, color = (:black,0.99))
                end
        end

        # Create axis for the histogram 
        ax2 = Axis(fig[1,6],
                # Background
                #backgroundcolor = :transparent,
                xgridvisible = false,
                ygridvisible = false,
                #limits = ((time[1],time[end]), nothing),
                # Title
                title = L"\mu = %$parameter",
                titlevisible = false,
                titlesize = 22,
                titlealign = :center,
                titlegap = 4.0,
                # Axes labels
                xlabel = "time",
                ylabel = L"x",
                xlabelvisible = false,
                ylabelvisible = false,
                xlabelsize = 20,
                ylabelsize = 20,
                xlabelcolor = :black,
                ylabelcolor = :black,
                xlabelpadding = 0.0,
                ylabelpadding = 0.0,
                # Axes scale, position and direction
                xscale = identity, 
                yscale = identity, #log10,
                xreversed = false,
                yreversed = false,
                xaxisposition = :bottom,
                yaxisposition = :left,
                # Ticks
                #xticks = [time[1], time[end]],
                #yticks = equilibria[:,1],
                xticksvisible = true,
                yticksvisible = true,
                xticksize = 6,
                yticksize = 6,
                # Ticks labels
                xticklabelsvisible = true,
                yticklabelsvisible = true,
                xticklabelsize = 18,
                yticklabelsize = 18,
                xticklabelalign = (:right, :top),
                yticklabelalign = (:right, :center),
                xticklabelcolor = :black,
                yticklabelcolor = :black,
                xtickformat = "{:.0f}",
                ytickformat = "{:.2f}",
        )
        # Plot the histogram of the stochastic process
        hist!(ax2, states, bins = 1000, normalization = :pdf, color = :red, strokecolor = :black, strokewidth = 1, direction = :x)
        # Plot the stationary FP distribution
        x_range = LinRange(-2.5*noise,2.5*noise,1000)
        #lines!(ax2, [gauss(x, parameter, noise) for x in x_range], x_range, color = :blue, linewidth = 2)
 
        # Export the figure
        save("./timeseries.png", fig)
end

# Plot the histogram of the distribution of the realizations and the stationary FP solution
using Statistics
function plot_distribution(bins, distribution, parameter::Float64, noise::Float64; index=1::Int64)
        # Create the figure 
        CairoMakie.activate!(; px_per_unit = 2)
        fig = Figure(; size = (1600, 1200), backgroundcolor = "#eeeeeeff")

        # Customize the appearence
        ax = Axis(fig[1,1],
                # Background
                backgroundcolor = :white,
                spinewidth = 5.0,
                xgridvisible = false,
                ygridvisible = false,
                #limits = ((bins[1],bins[end]), nothing),
                limits = ((bins[1],bins[end]), (-0.5,6.5)),
                # Title
                title = L"\mu = %$parameter",
                titlevisible = false,
                titlesize = 22,
                titlealign = :center,
                titlegap = 4.0,
                # Axes labels
                xlabel = L"\mathbf{x}",
                ylabel = L"\textbf{density}",
                xlabelvisible = true,
                ylabelvisible = true,
                xlabelsize = 50,
                ylabelsize = 50,
                xlabelcolor = :black,
                ylabelcolor = :black,
                xlabelpadding = -50.0,
                ylabelpadding = -50.0,
                # Axes scale, position and direction
                xscale = identity, 
                yscale = identity, #log10,
                xreversed = false,
                yreversed = false,
                xaxisposition = :bottom,
                yaxisposition = :left,
                # Ticks
                xticks = [2.1,3.9],
                yticks = [0,6],
                xticksvisible = true,
                yticksvisible = true,
                xticksize = 20,
                yticksize = 20,
                xtickwidth = 5.0,
                ytickwidth = 5.0,
                # Ticks labels
                xticklabelsvisible = true,
                yticklabelsvisible = true,
                xticklabelsize = 50,
                yticklabelsize = 50,
                xticklabelalign = (:right, :top),
                yticklabelalign = (:right, :center),
                xticklabelcolor = :black,
                yticklabelcolor = :black,
                xtickformat = "{:.0f}",
                ytickformat = "{:.0f}",
        )
        # Plot the binned distribution of the realizations
        scatter!(ax, bins, distribution, color = (:gray69,0.01), strokecolor = :black, strokewidth = 1, markersize = 25)
        # Plot the stationary FP solution
        #lines!(ax, bins, [gauss(x, parameter, noise) for x in bins], linewidth = 3, color = :red)
        
        # Get the mean of the distributions
        μ = Statistics.mean(bins, weights(distribution)) 

        # Extract the datapoints of the distribuzion within 1 standard deviation off the mean
        idx = findall(z -> z > μ-noise && z < μ+noise, bins)
        filtered_bins = [bins[n] for n in idx]
        filtered_dist = [distribution[n] for n in idx]

        # Plot the filtered datapoints
        scatter!(ax, filtered_bins, filtered_dist, color = (:brown2,0.35), strokecolor = :black, strokewidth = 1, markersize = 25)

        # Export the figure
        save("../../../results/precursors/escape_problem/distribution/$index.png", fig)
end

# Plot the potential function and its reconstruction from the stochastic process
function plot_potential(domain, U::Function, Ux::Function, Uxx::Function, bins, distribution, interpolant::Polynomial, parameter::Float64, noise::Float64; index=1::Int64)
        # Compute diffusion coefficient
        D = (noise^2)/2

        # Get the mean of the distribution
        μ = mean(distribution) 

        # Compute the numerical derivatives of the potential distribution
        ddistribution = approx_derivative(bins, distribution) 
        dddistribution = approx_derivative(bins[2:end], ddistribution) 

        # Extract the datapoints of the distribuzion within a certain range
        idx = findall(z -> z > μ-noise && z < μ+noise, bins)
        filtered_bins = [bins[n] for n in idx]
        filtered_dist = [distribution[n] for n in idx]
        filtered_ddist = [ddistribution[n] for n in idx[1:(end-1)]]
        filtered_dddist = [dddistribution[n] for n in idx[1:(end-2)]]

        # Create the figure 
        CairoMakie.activate!(; px_per_unit = 2)
        fig = Figure(; size = (1600, 1200), backgroundcolor = "#eeeeeeff")

        # Create axis for V(x) 
        ax1 = Axis(fig[1,1],
                # Background
                backgroundcolor = :white,
                spinewidth = 5.0,
                xgridvisible = false,
                ygridvisible = false,
                #limits = ((domain[1],domain[end]), (U(domain[end], parameter), U(domain[1], parameter))),
                limits = ((domain[1],domain[end]), (-15,15)),
                # Title
                title = L"\mu = %$parameter",
                titlevisible = false,
                titlesize = 22,
                titlealign = :center,
                titlegap = 4.0,
                # Axes labels
                xlabel = L"\mathbf{x}",
                ylabel = L"\mathbf{V(x)}",
                xlabelvisible = true,
                ylabelvisible = true,
                xlabelsize = 50,
                ylabelsize = 50,
                xlabelcolor = :black,
                ylabelcolor = :black,
                xlabelpadding = -50.0,
                ylabelpadding = -50.0,
                # Axes scale, position and direction
                xscale = identity, 
                yscale = identity, #log10,
                xreversed = false,
                yreversed = false,
                xaxisposition = :bottom,
                yaxisposition = :left,
                # Ticks
                xticks = [-1.9,4.9],
                yticks = [-14,14],
                xticksvisible = true,
                yticksvisible = true,
                xticksize = 20,
                yticksize = 20,
                xtickwidth = 5.0,
                ytickwidth = 5.0,
                # Ticks labels
                xticklabelsvisible = true,
                yticklabelsvisible = true,
                xticklabelsize = 50,
                yticklabelsize = 50,
                xticklabelalign = (:right, :top),
                yticklabelalign = (:right, :center),
                xticklabelcolor = :black,
                yticklabelcolor = :black,
                xtickformat = "{:.0f}",
                ytickformat = "{:.0f}",
        )
        # Plot the analytical potential function
        lines!(ax1, domain, [U(x, parameter) for x in domain], color = (:black, 0.25), linewidth = 4)
        # Plot the potential function from the analytical solution of the stationary FPE
        #scatter!(ax1, bins, [-D*log(gauss(x, parameter, noise)/5.6419) for x in bins], color = (:red,0.5), strokecolor = :black, strokewidth = 1, markersize = 20)
       # Plot the approximate reconstruction of the potential from the interpolation
        lines!(ax1, domain, [interpolant(x) for x in domain], color = (:brown2,0.99), linewidth = 4)
        # Plot the potential distribution of the realizations of the stochastic process
        scatter!(ax1, bins, distribution, marker = :xcross, color = (:brown2,0.35), strokecolor = :black, strokewidth = 1, markersize = 25)
        # Plot the filtered datapoints
        #scatter!(ax1, filtered_bins, filtered_dist, marker = :xcross, color = (:brown2,0.35), strokecolor = :black, strokewidth = 1, markersize = 25)

        #=
        # Create axis for Vx(x) 
        ax2 = Axis(fig[1,2],
                # Background
                #backgroundcolor = :transparent,
                xgridvisible = false,
                ygridvisible = false,
                limits = ((bins[1],bins[end]), nothing),
                # Title
                title = L"\mu = %$parameter",
                titlevisible = true,
                titlesize = 22,
                titlealign = :center,
                titlegap = 4.0,
                # Axes labels
                xlabel = "time",
                ylabel = L"x",
                xlabelvisible = false,
                ylabelvisible = false,
                xlabelsize = 20,
                ylabelsize = 20,
                xlabelcolor = :black,
                ylabelcolor = :black,
                xlabelpadding = 0.0,
                ylabelpadding = 0.0,
                # Axes scale, position and direction
                xscale = identity, 
                yscale = identity, #log10,
                xreversed = false,
                yreversed = false,
                xaxisposition = :bottom,
                yaxisposition = :left,
                # Ticks
                #xticks = equilibria[:,1],
                #yticks = [time[1], time[end]],
                xticksvisible = true,
                yticksvisible = true,
                xticksize = 6,
                yticksize = 6,
                # Ticks labels
                xticklabelsvisible = true,
                yticklabelsvisible = true,
                xticklabelsize = 18,
                yticklabelsize = 18,
                xticklabelalign = (:right, :top),
                yticklabelalign = (:right, :center),
                xticklabelcolor = :black,
                yticklabelcolor = :black,
                xtickformat = "{:.2f}",
                ytickformat = "{:.3f}",
        )
        # Plot the analytical first derivative
        lines!(ax2, bins, [Ux(x, parameter) for x in bins], color = (:black, 0.99), linewidth = 4)
        # Plot the first derivative from the interpolation
        lines!(ax2, bins, [(derivative(interpolant))(x) for x in bins], color = (:blue,0.5), linewidth = 3, linestyle = :dash)
        # Plot the approximate first derivative from the potential distribution
        scatter!(ax2, bins[2:end], ddistribution, marker = :xcross, color = (:gray69,0.5), strokecolor = :black, strokewidth = 1, markersize = 25)
        # Plot the filtered datapoints
        scatter!(ax2, filtered_bins[2:end], filtered_ddist, marker = :xcross, color = (:yellow,0.5), strokecolor = :black, strokewidth = 1, markersize = 25)

         # Create axis for Vxx(x) 
        ax3 = Axis(fig[1,3],
                # Background
                #backgroundcolor = :transparent,
                xgridvisible = false,
                ygridvisible = false,
                limits = ((bins[1],bins[end]), nothing),
                # Title
                title = L"\mu = %$parameter",
                titlevisible = true,
                titlesize = 22,
                titlealign = :center,
                titlegap = 4.0,
                # Axes labels
                xlabel = "time",
                ylabel = L"x",
                xlabelvisible = false,
                ylabelvisible = false,
                xlabelsize = 20,
                ylabelsize = 20,
                xlabelcolor = :black,
                ylabelcolor = :black,
                xlabelpadding = 0.0,
                ylabelpadding = 0.0,
                # Axes scale, position and direction
                xscale = identity, 
                yscale = identity, #log10,
                xreversed = false,
                yreversed = false,
                xaxisposition = :bottom,
                yaxisposition = :left,
                # Ticks
                #xticks = equilibria[:,1],
                #yticks = [time[1], time[end]],
                xticksvisible = true,
                yticksvisible = true,
                xticksize = 6,
                yticksize = 6,
                # Ticks labels
                xticklabelsvisible = true,
                yticklabelsvisible = true,
                xticklabelsize = 18,
                yticklabelsize = 18,
                xticklabelalign = (:right, :top),
                yticklabelalign = (:right, :center),
                xticklabelcolor = :black,
                yticklabelcolor = :black,
                xtickformat = "{:.2f}",
                ytickformat = "{:.3f}",
        )
        # Plot the analytical second derivative 
        lines!(ax3, bins, [Uxx(x, parameter) for x in bins], color = (:black, 0.99), linewidth = 4)
        # Plot the second derivative from the interpolation
        lines!(ax3, bins, [(derivative(derivative(interpolant)))(x) for x in bins], color = (:blue,0.5), linewidth = 3, linestyle = :dash)
        # Plot the approximate second derivative from the potential distribution 
        scatter!(ax3, bins[3:end], dddistribution, marker = :xcross, color = (:gray69,0.5), strokecolor = :black, strokewidth = 1, markersize = 25)
        # Plot the filtered datapoints
        scatter!(ax3, filtered_bins[3:end], filtered_dddist, marker = :xcross, color = (:yellow,0.5), strokecolor = :black, strokewidth = 1, markersize = 25)
        =#
 
        # Export the figure
        save("../../../results/precursors/escape_problem/potential/$index.png", fig)
end

# Plot the timeseries of the ensemble
function plot_ensemble(equilibria, parameter::Float64, idx::Int64)
        # Read the file containing the timestamps
        df = DataFrame(CSV.File("../data/time.csv"; delim=',', header=false))
        time = df[!,1]
        
        # Read the file containing the ensemble
        df = DataFrame(CSV.File("../data/ensemble/$idx.csv"; delim=',', header=false))
        ensemble = df|>CSV.Tables.matrix

        # Get the number of sample paths
        Ne = nrow(df)
        # Get the number of realizations
        Nt = ncol(df)

        # Create the figure 
        CairoMakie.activate!(; px_per_unit = 2)
        fig = Figure(; size = (1200, 800))#, backgroundcolor = :transparent)

        # Customise the apperance
        ax = Axis(fig[1,1],
                # Background
                #backgroundcolor = :transparent,
                xgridvisible = false,
                ygridvisible = false,
                limits = ((time[1],time[end]), nothing),
                # Title
                title = L"\mu = %$parameter",
                titlevisible = true,
                titlesize = 22,
                titlealign = :center,
                titlegap = 4.0,
                # Axes labels
                xlabel = "time",
                ylabel = L"x",
                xlabelvisible = false,
                ylabelvisible = false,
                xlabelsize = 20,
                ylabelsize = 20,
                xlabelcolor = :black,
                ylabelcolor = :black,
                xlabelpadding = 0.0,
                ylabelpadding = 0.0,
                # Axes scale, position and direction
                xscale = identity, 
                yscale = identity, #log10,
                xreversed = false,
                yreversed = false,
                xaxisposition = :bottom,
                yaxisposition = :left,
                # Ticks
                xticks = [time[1], time[end]],
                yticks = equilibria[:,1],
                xticksvisible = true,
                yticksvisible = true,
                xticksize = 6,
                yticksize = 6,
                # Ticks labels
                xticklabelsvisible = true,
                yticklabelsvisible = true,
                xticklabelsize = 18,
                yticklabelsize = 18,
                xticklabelalign = (:right, :top),
                yticklabelalign = (:right, :center),
                xticklabelcolor = :black,
                yticklabelcolor = :black,
                xtickformat = "{:.0f}",
                ytickformat = "{:.2f}",
        )
        # Plot the ensemble's individual sample paths
        for n in 1:Ne
                lines!(ax, time, ensemble[n,:], linewidth = 0.5, color = (:red,0.15))
        end
        # Plot the ensemble's mean timeseries
        lines!(ax, time, [mean(ensemble[:,t]) for t in 1:Nt], linewidth = 2.5, color = (:blue,0.8))
        # Plot the equilibria of the vector field
        for n in 1:length(equilibria[:,1])
                if equilibria[n,2]==1
                        lines!(ax, [time[1], time[end]], [equilibria[n,1], equilibria[n,1]], linewidth = 1, color = (:black,0.99))
                else
                        lines!(ax, [time[1], time[end]], [equilibria[n,1], equilibria[n,1]], linewidth = 1, color = (:black,0.99), linestyle = :dash)
                end
        end

        # Export the figure
        save("../../../results/precursors/escape_problem/timeseries/$idx.png", fig)
end

# Plot the hitting time distribution
function plot_escapes(time, distribution, parameter, plt_idx)
        # Compute the mean of the hitting times distribution
        T = mean(distribution)

        # Create the figure 
        CairoMakie.activate!(; px_per_unit = 2)
        fig = Figure(; size = (1200, 800))#, backgroundcolor = :transparent)

        # Customise the apperance
        ax = Axis(fig[1,1],
                # Background
                #backgroundcolor = :transparent,
                xgridvisible = false,
                ygridvisible = false,
                limits = ((time[1],time[end]), nothing),
                # Title
                title = L"\mu = %$parameter",
                titlevisible = true,
                titlesize = 22,
                titlealign = :center,
                titlegap = 4.0,
                # Axes labels
                xlabel = "time",
                ylabel = L"x",
                xlabelvisible = false,
                ylabelvisible = false,
                xlabelsize = 20,
                ylabelsize = 20,
                xlabelcolor = :black,
                ylabelcolor = :black,
                xlabelpadding = 0.0,
                ylabelpadding = 0.0,
                # Axes scale, position and direction
                xscale = identity, 
                yscale = identity, #log10,
                xreversed = false,
                yreversed = false,
                xaxisposition = :bottom,
                yaxisposition = :left,
                # Ticks
                xticks = [T, time[end]],
                #yticks = equilibria[:,1],
                xticksvisible = true,
                yticksvisible = true,
                xticksize = 6,
                yticksize = 6,
                # Ticks labels
                xticklabelsvisible = true,
                yticklabelsvisible = true,
                xticklabelsize = 18,
                yticklabelsize = 18,
                xticklabelalign = (:right, :top),
                yticklabelalign = (:right, :center),
                xticklabelcolor = :black,
                yticklabelcolor = :black,
                xtickformat = "{:.1f}",
                ytickformat = "{:.2f}",
        )
        # Plot the histogram of the distribution
        hist!(ax, distribution, bins = 200, normalization = :pdf, color = :red, strokecolor = :black, strokewidth = 1)
        # Plot the mean first exit time (MFET)
        scatter!(ax, T, 0, color = :yellow, strokecolor = :black, strokewidth = 1, markersize = 25)

        # Export the figure
        save("../../../results/precursors/escape_problem/hitting_times/$plt_idx.png", fig)
end
