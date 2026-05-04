"""
    Plotting script

Functions used to create the plots in each figure.
"""

# Create empty layouts for the figures
include("./figs.jl")

# Plot the ground-truth domain = LinRange(-1.5,1.5,1000)
lines!(ax2, domain, [U(x,μ) for x in domain], color = (:red,1.0), linewidth = 3.0)

# Plot the timeseries and the solution of the nonlinear reconstruction of the potential
function plot_solutions(timesteps, timeseries, shifted_potential, idx)
        # Extract the length of the timeseries
        Nt = length(timeseries)

        # Plot the timeseries only for the first 10 particles 
        if idx <= 100
                # Plot the timeseries
                lines!(ax1, timesteps[1:1000], timeseries[1:1000], color = (:black,0.15), linewidth = 3.0)
        end
        
        # Plot the reconstructed shifted potential
        lines!(ax2, domain, [shifted_potential(x) for x in domain], color = (:black,0.5), linewidth = 2.0)

        # Setup ticks and limits for the plot
        y_range = maximum(timeseries) - minimum(timeseries)
        ax1.limits = ((timesteps[1], timesteps[1000]), (minimum(timeseries) - 0.1*y_range, maximum(timeseries) + 0.1*y_range))
        ax1.xticks = [timesteps[1], timesteps[1000]]
        ax1.yticks = [minimum(timeseries), maximum(timeseries)]
end

function plot_results(coefficients, transformations)
        # Compute the number of bins for the ensemble
        #n_bins = convert(Int64, floor(0.02*Ne))
        n_bins = 20

        # Extract the single coefficients 
        mat = transpose(reduce(hcat, coefficients))     # Turns a vector of vectors into a matrix
        c1, c2, c3 = eachcol(mat)                       # Turns the columns of the matrix into vectors
 
        # Fit and plot a histogram for c1 
        bins, pdf = fit_distribution(c1, n_bins=n_bins)
        barplot!(ax3, bins, pdf, color = pdf, colormap = [(CtpYellow,0.5),(CtpMauve,0.5)], strokecolor = :black, strokewidth = 1)
        # Plot the c1 value for the ground truth
        #lines!(ax3, [μ,μ], [0,1.1*maximum(pdf)], color = CtpRed, linestyle = :dash, linewidth = 5.0)
        # Setup ticks and limits for the plot
        x_range = bins[end]-bins[1] 
        ax3.limits = ((bins[1] - 0.1*x_range, bins[end] + 0.1*x_range), (0, 1.1*maximum(pdf)))
        ax3.xticks = [bins[1],bins[end]]
        ax3.yticks = [0,maximum(pdf)]

        # Fit and plot a histogram for c2
        bins, pdf = fit_distribution(c2, n_bins=n_bins)
        barplot!(ax4, bins, pdf, color = pdf, colormap = [(CtpYellow,0.5),(CtpMauve,0.5)], strokecolor = :black, strokewidth = 1)
        # Plot the c2 value for the ground truth
        #lines!(ax4, [0,0], [0,1.1*maximum(pdf)], color = CtpRed, linestyle = :dash, linewidth = 5.0)
        # Setup ticks and limits for the plot
        x_range = bins[end]-bins[1] 
        ax4.limits = ((bins[1] - 0.1*x_range, bins[end] + 0.1*x_range), (0, 1.1*maximum(pdf)))
        ax4.xticks = [bins[1],bins[end]]
        ax4.yticks = [0,maximum(pdf)]

        # Fit and plot a histogram for c3
        bins, pdf = fit_distribution(c3, n_bins=n_bins)
        barplot!(ax5, bins, pdf, color = pdf, colormap = [(CtpYellow,0.5),(CtpMauve,0.5)], strokecolor = :black, strokewidth = 1)
        # Plot the c3 value for the ground truth
        #lines!(ax5, [1.0/3.0,1.0/3.0], [0,1.1*maximum(pdf)], color = CtpRed, linestyle = :dash, linewidth = 5.0)
        # Setup ticks and limits for the plot
        x_range = bins[end]-bins[1] 
        ax5.limits = ((bins[1] - 0.1*x_range, bins[end] + 0.1*x_range), (0, 1.1*maximum(pdf)))
        ax5.xticks = [bins[1],bins[end]]
        ax5.yticks = [0,maximum(pdf)]

        # Extract the single transformations 
        mat = transpose(reduce(hcat, transformations))
        xs, Vs, Vxxs, LDPs, xu, Vu, Vxxu, LDPu, ΔV, C, LDP, E = eachcol(mat)

        # Fit and plot a histogram for xs 
        bins, pdf = fit_distribution(xs, n_bins=n_bins)
        barplot!(ax6, bins, pdf, color = pdf, colormap = [(CtpYellow,0.5),(CtpMauve,0.5)], strokecolor = :black, strokewidth = 1)
        # Plot the xs value for the ground truth
        #lines!(ax6, [sqrt(-μ),sqrt(-μ)], [0,1.1*maximum(pdf)], color = CtpRed, linestyle = :dash, linewidth = 5.0)
        # Setup ticks and limits for the plot
        x_range = bins[end]-bins[1] 
        ax6.limits = ((bins[1] - 0.1*x_range, bins[end] + 0.1*x_range), (0, 1.1*maximum(pdf)))
        ax6.xticks = [bins[1],bins[end]]
        ax6.yticks = [0,maximum(pdf)]

        # Fit and plot a histogram for Vs 
        bins, pdf = fit_distribution(Vs, n_bins=n_bins)
        barplot!(ax7, bins, pdf, color = pdf, colormap = [(CtpYellow,0.5),(CtpMauve,0.5)], strokecolor = :black, strokewidth = 1)
        # Plot the Vs value for the ground truth
        #lines!(ax7, [U(sqrt(-μ),μ),U(sqrt(-μ),μ)], [0,1.1*maximum(pdf)], color = CtpRed, linestyle = :dash, linewidth = 5.0)
        # Setup ticks and limits for the plot
        x_range = bins[end]-bins[1] 
        ax7.limits = ((bins[1] - 0.1*x_range, bins[end] + 0.1*x_range), (0, 1.1*maximum(pdf)))
        ax7.xticks = [bins[1],bins[end]]
        ax7.yticks = [0,maximum(pdf)]
        #labels = ["0", "6e5"]
        #ax7.ytickformat = values -> ["$(label)" for label in labels]

        # Fit and plot a histogram for Vxxs 
        bins, pdf = fit_distribution(Vxxs, n_bins=n_bins)
        barplot!(ax8, bins, pdf, color = pdf, colormap = [(CtpYellow,0.5),(CtpMauve,0.5)], strokecolor = :black, strokewidth = 1)
        # Plot the Vxxs value for the ground truth
        #lines!(ax8, [Uxx(sqrt(-μ),μ),Uxx(sqrt(-μ),μ)], [0,1.1*maximum(pdf)], color = CtpRed, linestyle = :dash, linewidth = 5.0)
        # Setup ticks and limits for the plot
        x_range = bins[end]-bins[1] 
        ax8.limits = ((bins[1] - 0.1*x_range, bins[end] + 0.1*x_range), (0, 1.1*maximum(pdf)))
        ax8.xticks = [bins[1],bins[end]]
        ax8.yticks = [0,maximum(pdf)]

        # Fit and plot a histogram for LDPs 
        bins, pdf = fit_distribution(LDPs, n_bins=n_bins)
        barplot!(ax9, bins, pdf, color = pdf, colormap = [(CtpYellow,0.5),(CtpMauve,0.5)], strokecolor = :black, strokewidth = 1)
        # Plot the LDPs value for the ground truth
        #true_LDPs = exp(U(sqrt(-μ),μ))
        #lines!(ax9, [true_LDPs,true_LDPs], [0,1.1*maximum(pdf)], color = CtpRed, linestyle = :dash, linewidth = 5.0)
        # Setup ticks and limits for the plot
        x_range = bins[end]-bins[1] 
        ax9.limits = ((bins[1] - 0.1*x_range, bins[end] + 0.1*x_range), (0, 1.1*maximum(pdf)))
        ax9.xticks = [bins[1],bins[end]]
        ax9.yticks = [0,maximum(pdf)]

        # Fit and plot a histogram for xu 
        bins, pdf = fit_distribution(xu, n_bins=n_bins)
        barplot!(ax10, bins, pdf, color = pdf, colormap = [(CtpYellow,0.5),(CtpMauve,0.5)], strokecolor = :black, strokewidth = 1)
        # Plot the xu value for the ground truth
        #lines!(ax10, [-sqrt(-μ),-sqrt(-μ)], [0,1.1*maximum(pdf)], color = CtpRed, linestyle = :dash, linewidth = 5.0)
        # Setup ticks and limits for the plot
        x_range = bins[end]-bins[1] 
        ax10.limits = ((bins[1] - 0.1*x_range, bins[end] + 0.1*x_range), (0, 1.1*maximum(pdf)))
        ax10.xticks = [bins[1],bins[end]]
        ax10.yticks = [0,maximum(pdf)]

        # Fit and plot a histogram for Vu 
        bins, pdf = fit_distribution(Vu, n_bins=n_bins)
        barplot!(ax11, bins, pdf, color = pdf, colormap = [(CtpYellow,0.5),(CtpMauve,0.5)], strokecolor = :black, strokewidth = 1)
        # Plot the Vu value for the ground truth
        #lines!(ax11, [U(-sqrt(-μ),μ),U(-sqrt(-μ),μ)], [0,1.1*maximum(pdf)], color = CtpRed, linestyle = :dash, linewidth = 5.0)
        # Setup ticks and limits for the plot
        x_range = bins[end]-bins[1] 
        ax11.limits = ((bins[1] - 0.1*x_range, bins[end] + 0.1*x_range), (0, 1.1*maximum(pdf)))
        ax11.xticks = [bins[1],bins[end]]
        ax11.yticks = [0,maximum(pdf)]

        # Fit and plot a histogram for Vxxu 
        bins, pdf = fit_distribution(Vxxu, n_bins=n_bins)
        barplot!(ax12, bins, pdf, color = pdf, colormap = [(CtpYellow,0.5),(CtpMauve,0.5)], strokecolor = :black, strokewidth = 1)
        # Plot the Vxx value for the ground truth
        #lines!(ax12, [Uxx(-sqrt(-μ),μ),Uxx(-sqrt(-μ),μ)], [0,1.1*maximum(pdf)], color = CtpRed, linestyle = :dash, linewidth = 5.0)
        # Setup ticks and limits for the plot
        x_range = bins[end]-bins[1] 
        ax12.limits = ((bins[1] - 0.1*x_range, bins[end] + 0.1*x_range), (0, 1.1*maximum(pdf)))
        ax12.xticks = [bins[1],bins[end]]
        ax12.yticks = [0,maximum(pdf)]

        # Fit and plot a histogram for LDPu
        bins, pdf = fit_distribution(LDPu, n_bins=n_bins)
        barplot!(ax13, bins, pdf, color = pdf, colormap = [(CtpYellow,0.5),(CtpMauve,0.5)], strokecolor = :black, strokewidth = 1)
        # Plot the LDPu value for the ground truth
        #true_LDPu = exp(U(-sqrt(-μ),μ))
        #lines!(ax13, [true_LDPu,true_LDPu], [0,1.1*maximum(pdf)], color = CtpRed, linestyle = :dash, linewidth = 5.0)
        # Setup ticks and limits for the plot
        x_range = bins[end]-bins[1] 
        ax13.limits = ((bins[1] - 0.1*x_range, bins[end] + 0.1*x_range), (0, 1.1*maximum(pdf)))
        ax13.xticks = [bins[1],bins[end]]
        ax13.yticks = [0,maximum(pdf)]

        # Fit and plot a histogram for ΔV 
        bins, pdf = fit_distribution(ΔV, n_bins=n_bins)
        barplot!(ax14, bins, pdf, color = pdf, colormap = [(CtpYellow,0.5),(CtpMauve,0.5)], strokecolor = :black, strokewidth = 1)
        # Plot the ΔV value for the ground truth
        true_ΔV = abs(U(-sqrt(-μ),μ) - U(sqrt(-μ),μ))
        lines!(ax14, [true_ΔV,true_ΔV], [0,1.1*maximum(pdf)], color = CtpRed, linestyle = :dash, linewidth = 5.0)
        # Setup ticks and limits for the plot
        x_range = bins[end]-bins[1] 
        ax14.limits = ((bins[1] - 0.1*x_range, bins[end] + 0.1*x_range), (0, 1.1*maximum(pdf)))
        ax14.xticks = [bins[1],bins[end]]
        ax14.yticks = [0,maximum(pdf)]

        # Fit and plot a histogram for C (prefactor) 
        bins, pdf = fit_distribution(C, n_bins=n_bins)
        barplot!(ax15, bins, pdf, color = pdf, colormap = [(CtpYellow,0.5),(CtpMauve,0.5)], strokecolor = :black, strokewidth = 1)
        # Plot the C value for the ground truth
        #true_prefactor = sqrt(abs(Uxx(-sqrt(-μ),μ))*Uxx(sqrt(-μ),μ))
        #lines!(ax15, [true_prefactor,true_prefactor], [0,1.1*maximum(pdf)], color = CtpRed, linestyle = :dash, linewidth = 5.0)
        # Setup ticks and limits for the plot
        x_range = bins[end]-bins[1] 
        ax15.limits = ((bins[1] - 0.1*x_range, bins[end] + 0.1*x_range), (0, 1.1*maximum(pdf)))
        ax15.xticks = [bins[1],bins[end]]
        ax15.yticks = [0,maximum(pdf)]

        # Fit and plot a histogram for LDP (exponential decay)
        bins, pdf = fit_distribution(LDP, n_bins=n_bins)
        barplot!(ax16, bins, pdf, color = pdf, colormap = [(CtpYellow,0.5),(CtpMauve,0.5)], strokecolor = :black, strokewidth = 1)
        # Plot the LDP value for the ground truth
        true_decay = exp(-true_ΔV)
        lines!(ax16, [true_decay,true_decay], [0,1.1*maximum(pdf)], color = CtpRed, linestyle = :dash, linewidth = 5.0)
        # Setup ticks and limits for the plot
        x_range = bins[end]-bins[1] 
        ax16.limits = ((bins[1] - 0.1*x_range, bins[end] + 0.1*x_range), (0, 1.1*maximum(pdf)))
        ax16.xticks = [bins[1],bins[end]]
        ax16.yticks = [0,maximum(pdf)]

        # Fit and plot a histogram for E (numerical error)
        bins, pdf = fit_distribution(E, n_bins=n_bins)
        barplot!(ax17, bins, pdf, color = pdf, colormap = [(CtpYellow,0.5),(CtpMauve,0.5)], strokecolor = :black, strokewidth = 1)
        # Setup ticks and limits for the plot
        x_range = bins[end]-bins[1] 
        ax17.limits = ((bins[1] - 0.1*x_range, bins[end] + 0.1*x_range), (0, 1.1*maximum(pdf)))
        ax17.xticks = [0,bins[end]]
        ax17.yticks = [0,maximum(pdf)]
end
