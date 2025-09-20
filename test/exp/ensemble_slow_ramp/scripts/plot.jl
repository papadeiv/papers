"""
    Plotting script

Collection of all the functions used to generate the plots of the simulations.
"""

# Create empty layouts for the figures
include("./figs.jl")

# Plot the timeseries 
function plot_solutions(x_data, full, residual)
        # Set up the limits of the plot
        ax1.limits = ((x_data[1],x_data[end]), (minimum(full),maximum(full))) 
        ax2.limits = ((x_data[1],x_data[end]), (minimum(residual),maximum(residual))) 
        # Set up the ticks on the axis of the plot
        ax1.xticks = [x_data[1],x_data[end]]
        ax2.xticks = [x_data[1],x_data[end]]
        ax1.yticks = [minimum(full),maximum(full)]
        ax2.yticks = [minimum(residual),maximum(residual)]
        
        # Plot the solutions
        lines!(ax1, x_data, full, color = (:black,0.15), linewidth = 3.0)
        lines!(ax2, x_data, residual, color = (:black,0.15), linewidth = 3.0)
end

# Plot the drift of the QSE 
function plot_drift(x_data, trend)
        lines!(ax1, x_data, trend, color = (CtpYellow,1.00), linewidth = 5.0)
end

# Define the domain for plotting the potential 
domain = LinRange(0,3.5,1000)

# Plot the true scalar potential 
function plot_ground_truth(parameters)
        # Loop over the paramer values in the windowed subseries
        for n in 1:length(parameters)
                # Plot only every 100 values of the parameter
                if mod(n, 1000) == 0
                        # Determine transparency value
                        alpha = convert(Float64, n/length(parameters))
                        # Plot the ground truth
                        lines!(ax3, domain, [U(x,parameters[n]) for x in domain], color = (:red, alpha), linewidth = 1.0)
                end
        end
end

# Plot the reconstructed potential 
function plot_reconstruction(potential)
        lines!(ax3, domain, [potential(x) for x in domain], color = (:black,0.15), linewidth = 1.0)
end

function plot_results(coefficients, variables)
        # Compute the number of bins for the ensemble
        n_bins = Nb

        println("")
        println("--------- Moments of the NLLS solutions ---------")

        # Fit and plot a histogram to the c1 coefficients
        bins, pdf = fit_distribution(coefficients[:,1], n_bins=n_bins)
        barplot!(ax4, bins, pdf, color = pdf, colormap = [(CtpWhite,1.0),(CtpGray,0.35)], strokecolor = :black, strokewidth = 1)
        # Fit and plot an optimal normal distribution to the c1 coefficients 
        domain = LinRange(bins[1], bins[end], 1000)
        ρ = Distributions.pdf.(fit_mle(Normal, coefficients[:,1]), domain)
        lines!(ax4, domain, ρ, color = ρ, colormap = [CtpYellow, CtpMauve], linewidth = 5)
        # Print the mean and variance of the c1 coefficient
        println("c1:  mean = ", mean(coefficients[:,1]), ", var = ", var(coefficients[:,1]))
        # Set-up axis limits and ticks 
        #ax4.limits = ((minimum(bins),maximum(bins)), (0,1.05*maximum()))
        ax4.xticks = [minimum(bins),maximum(bins)]
        ax4.yticks = [0,1.05*maximum(pdf)]

        # Fit and plot a histogram to the c2 coefficients
        bins, pdf = fit_distribution(coefficients[:,2], n_bins=n_bins)
        barplot!(ax5, bins, pdf, color = pdf, colormap = [(CtpWhite,1.0),(CtpGray,0.35)], strokecolor = :black, strokewidth = 1)
        # Fit and plot an optimal Gamma distribution to the c2 coefficients 
        domain = LinRange(bins[1], bins[end], 1000)
        ρ = Distributions.pdf.(fit_mle(Gamma, abs.(coefficients[:,2])), domain)
        lines!(ax5, domain, ρ, color = ρ, colormap = [CtpYellow, CtpMauve], linewidth = 5)
        # Set-up axis limits and ticks 
        #ax5.limits = ((minimum(bins),maximum(bins)), (0,1.05*maximum()))
        ax5.xticks = [minimum(bins),maximum(bins)]
        ax5.yticks = [0,1.05*maximum(pdf)]
        # Print the mean and variance of the c2 coefficient
        println("c2:  mean = ", mean(coefficients[:,2]), ", var = ", var(coefficients[:,2]))

        # Fit and plot a histogram to the c3 coefficients
        bins, pdf = fit_distribution(coefficients[:,3], n_bins=n_bins)
        barplot!(ax6, bins, pdf, color = pdf, colormap = [(CtpWhite,1.0),(CtpGray,0.35)], strokecolor = :black, strokewidth = 1)
        # Fit and plot an optimal Gamma distribution
        domain = LinRange(bins[1], bins[end], 1000)
        ρ = Distributions.pdf.(fit_mle(Laplace, coefficients[:,3]), domain)
        lines!(ax6, domain, ρ, color = ρ, colormap = [CtpYellow, CtpMauve], linewidth = 5)
        # Set-up axis limits and ticks 
        #ax6.limits = ((minimum(bins),maximum(bins)), (0,1.05*maximum()))
        ax6.xticks = [minimum(bins),maximum(bins)]
        ax6.yticks = [0,1.05*maximum(pdf)]
        # Print the mean and variance of the c3 coefficient
        println("c3:  mean = ", mean(coefficients[:,3]), ", var = ", var(coefficients[:,3]))

        println("")
        println("--------- Transformations of the NLLS solutions ---------")
 
        # Extract and name the variables for clarity
        mat = transpose(reduce(hcat, variables))                             # Turns a vector of vectors into a matrix
        xs, xu, Vs, Vu, Vxxs, Vxxu, ΔV, LDPs, LDPu, LDP = eachcol(mat)       # Turns the columns of the matrix into vectors

        # Fit and plot a histogram to the stable equilibria 
        bins, pdf = fit_distribution(xs, n_bins=n_bins)
        barplot!(ax7, bins, pdf, color = pdf, colormap = [(CtpWhite,1.0),(CtpGray,0.35)], strokecolor = :black, strokewidth = 1)
        # Fit and plot an optimal normal distribution to the stable equilibria
        domain = LinRange(bins[1], bins[end], 1000)
        ρ = Distributions.pdf.(fit_mle(Normal, xs), domain)
        lines!(ax7, domain, ρ, color = ρ, colormap = [CtpYellow, CtpMauve], linewidth = 5)
        # Set-up axis limits and ticks 
        #ax7.limits = ((minimum(bins),maximum(bins)), (0,1.05*maximum()))
        ax7.xticks = [minimum(bins),maximum(bins)]
        ax7.yticks = [0,1.05*maximum(pdf)]
        # Print the mean and variance of the stable equilibrium 
        println("xs:         mean = ", mean(xs), ", var = ", var(xs))

        # Clamp the data in the 5-95 percentile
        low, high = quantile(xu, [0.05, 0.95])
        xu_clipped = clamp.(xu, low, high)
        # Fit and plot a histogram to the unstable equilibria 
        bins, pdf = fit_distribution(xu_clipped, n_bins=n_bins)
        barplot!(ax8, bins, pdf, color = pdf, colormap = [(CtpWhite,1.0),(CtpGray,0.35)], strokecolor = :black, strokewidth = 1)
        # Set-up axis limits and ticks 
        #ax8.limits = ((minimum(bins),maximum(bins)), (0,1.05*maximum(pdf)))
        ax8.xticks = [minimum(bins),maximum(bins)]
        ax8.yticks = [0,1.05*maximum(pdf)]
        # Print the mean and variance of the stable equilibrium 
        println("xu:         mean = ", mean(xu), ", var = ", var(xu))

        # Fit and plot a histogram to the potential value at stable equilibrium
        bins, pdf = fit_distribution(Vs, n_bins=n_bins)
        barplot!(ax9, bins, pdf, color = pdf, colormap = [(CtpWhite,1.0),(CtpGray,0.35)], strokecolor = :black, strokewidth = 1)
        # Set-up axis limits and ticks 
        #ax9.limits = ((minimum(bins),maximum(bins)), (0,1.05*maximum(pdf)))
        ax9.xticks = [minimum(bins),maximum(bins)]
        ax9.yticks = [0,1.05*maximum(pdf)]
        # Print the mean and variance of the potential value 
        println("V(xs):      mean = ", mean(Vs), ", var = ", var(Vs))

        # Clamp the data in the 5-95 percentile
        low, high = quantile(Vu, [0.05, 0.95])
        Vu_clipped = clamp.(Vu, low, high)
        # Fit and plot a histogram to the potential value at unstable equilibriam
        bins, pdf = fit_distribution(Vu_clipped, n_bins=n_bins)
        barplot!(ax10, bins, pdf, color = pdf, colormap = [(CtpWhite,1.0),(CtpGray,0.35)], strokecolor = :black, strokewidth = 1)
        # Set-up axis limits and ticks 
        #ax10.limits = ((minimum(bins),maximum(bins)), (0,1.05*maximum(pdf)))
        ax10.xticks = [minimum(bins),maximum(bins)]
        ax10.yticks = [0,1.05*maximum(pdf)]
        # Print the mean and variance of the potential value 
        println("V(xu):      mean = ", mean(Vu), ", var = ", var(Vu))

        # Fit and plot a histogram to the curvature value at stable equilibrium
        bins, pdf = fit_distribution(Vxxs, n_bins=n_bins)
        barplot!(ax11, bins, pdf, color = pdf, colormap = [(CtpWhite,1.0),(CtpGray,0.35)], strokecolor = :black, strokewidth = 1)
        # Fit and plot an optimal Poisson distribution
        domain = LinRange(bins[1], bins[end], 1000)
        ρ = Distributions.pdf.(fit_mle(Gamma, Vxxs), domain)
        lines!(ax11, domain, ρ, color = ρ, colormap = [CtpYellow, CtpMauve], linewidth = 5)
        # Set-up axis limits and ticks 
        #ax11.limits = ((minimum(bins),maximum(bins)), (0,1.05*maximum()))
        ax11.xticks = [minimum(bins),maximum(bins)]
        ax11.yticks = [0,1.05*maximum(pdf)]
        # Print the mean and variance of the curvature value 
        println("V''(xs):    mean = ", mean(Vxxs), ", var = ", var(Vxxs))

        # Fit and plot a histogram to the curvature value at unstable equilibrium
        bins, pdf = fit_distribution(Vxxu, n_bins=n_bins)
        barplot!(ax12, bins, pdf, color = pdf, colormap = [(CtpWhite,1.0),(CtpGray,0.35)], strokecolor = :black, strokewidth = 1)
        # Set-up axis limits and ticks 
        #ax12.limits = ((minimum(bins),maximum(bins)), (0,1.05*maximum(pdf)))
        ax12.xticks = [minimum(bins),maximum(bins)]
        ax12.yticks = [0,1.05*maximum(pdf)]
        # Print the mean and variance of the curvature value 
        println("V''(xu):    mean = ", mean(Vxxu), ", var = ", var(Vxxu))

        # Fit and plot a histogram to the exponential of the potential value at stable equilibrium
        bins, pdf = fit_distribution(LDPs, n_bins=n_bins)
        barplot!(ax13, bins, pdf, color = pdf, colormap = [(CtpWhite,1.0),(CtpGray,0.35)], strokecolor = :black, strokewidth = 1)
        # Set-up axis limits and ticks 
        #ax13.limits = ((minimum(bins),maximum(bins)), (0,1.05*maximum(pdf)))
        ax13.xticks = [minimum(bins),maximum(bins)]
        ax13.yticks = [0,1.05*maximum(pdf)]
        # Print the mean and variance of the exponential of the potential value 
        println("e^(V(xs)):  mean = ", mean(LDPs), ", var = ", var(LDPs))

        #=
        # Plot the realisations of the potential value at unstable equilibrium
        scatter!(ax14, [0.0 for n in 1:length(LDPu)], LDPu, markersize = 10, color = :red, strokewidth = 1.0)
        # Fit and plot a histogram to the exponential of the potential value at unstable equilibrium
        bins, pdf = fit_distribution(LDPu, n_bins=n_bins)
        barplot!(ax14, bins, pdf, color = pdf, colormap = [(CtpWhite,1.0),(CtpGray,0.35)], strokecolor = :black, strokewidth = 1)
        # Set-up axis limits and ticks 
        ax14.limits = ((-0.2,0.2), nothing)
        ax14.xticks = [0]
        =#
        # Print the mean and variance of the exponential of the potential value 
        println("e^(V(xu)):  mean = ", mean(LDPu), ", var = ", var(LDPu))

        println("")
        println("--------- Large-deviation principles ---------")

        # Clamp the data in the 5-95 percentile
        low, high = quantile(ΔV, [0.05, 0.95])
        ΔV_clipped = clamp.(ΔV, low, high)
        # Fit and plot a histogram to the energy barrier 
        bins, pdf = fit_distribution(ΔV_clipped, n_bins=n_bins)
        barplot!(ax15, bins, pdf, color = pdf, colormap = [(CtpWhite,1.0),(CtpGray,0.35)], strokecolor = :black, strokewidth = 1)
        # Set-up axis limits and ticks 
        #ax15.limits = ((minimum(bins),maximum(bins)), (0,1.05*maximum(pdf)))
        ax15.xticks = [minimum(bins),maximum(bins)]
        ax15.yticks = [0,1.05*maximum(pdf)]
        # Print the mean and variance of the energy barrier 
        println("ΔV:  mean = ", mean(ΔV), ", var = ", var(ΔV))
 
        # Clamp the data in the 5-95 percentile
        low, high = quantile(LDP, [0.05, 0.95])
        LDP_clipped = clamp.(LDP, low, high)
        # Fit and plot a histogram to the large-deviation
        bins, pdf = fit_distribution(LDP_clipped, n_bins=n_bins)
        barplot!(ax16, bins, pdf, color = pdf, colormap = [(CtpWhite,1.0),(CtpGray,0.35)], strokecolor = :black, strokewidth = 1)
        # Set-up axis limits and ticks 
        #ax16.limits = ((minimum(bins),maximum(bins)), (0,1.05*maximum(pdf)))
        ax16.xticks = [minimum(bins),maximum(bins)]
        ax16.yticks = [0,1.05*maximum(pdf)]
        # Print the mean and variance of the energy barrier 
        println("LDP:  mean = ", mean(LDP_clipped), ", var = ", var(LDP_clipped))
        println("")

        writeCSV(hcat(sort(xu),sort(Vu)), "../../res/data/results.csv")
end

# Print the title with info
function print_info(Nt, idx)
        # Compute parameter range
        p0 = μ0[idx]
        pf = μf[idx]

        # Print the info
        Label(fig1[begin-1, 1:2],
              L"\textbf{N_t = %$Nt,\; N_e = %$Ne}",
              fontsize = 50,
              padding = (0,0,0,0),
             )
        Label(fig4[begin-1, 1:3],
              L"\mathbf{\varepsilon = %$ε,\;\mu\in[%$p0,%$pf],\;}\textbf{N_t = %$Nt,\; N_e = %$Ne}",
              fontsize = 50,
              padding = (0,0,0,0),
             )
        Label(fig7[begin-1, 1:4],
              L"\mathbf{\varepsilon = %$ε,\;\mu\in[%$p0,%$pf],\;}\textbf{N_t = %$Nt,\; N_e = %$Ne}",
              fontsize = 50,
              padding = (0,0,0,0),
             )
        Label(fig15[begin-1, 1:2],
              L"\mathbf{\varepsilon = %$ε,\;\mu\in[%$p0,%$pf],\;}\textbf{N_t = %$Nt,\; N_e = %$Ne}",
              fontsize = 50,
              padding = (0,0,0,0),
             )
end

# Plot the escape ews errorbars and analytic ground truth
function plot_ews(escape)
        # Define a "continuous" domain for the parameter range
        μ = collect(LinRange(μ0[1], μf[end], 1000))

        # Loop over the domain
        true_escape = Vector{Float64}(undef, length(μ))
        for n in 1:length(μ)
                # Get the stable and unstable equilibria of the ground truth
                bounds = get_bounds(μ[n], [0,2,1])

                # Define new scalar potential and its derivative
                ground_truth(x) = U(x,μ[n]) 
                second_derivative(x) = Uxx(x,μ[n]) 

                # Compute the escape ews of the ground truth
                true_escape[n] = (estimate_escape(ground_truth, second_derivative, bounds.true_min, bounds.true_max, σ)).LDP
        end

        # Plot the ground truth escape
        lines!(ax11, μ, true_escape, color = :red, linewidth = 3.0)

        #=
        # Define the limits and ticks of the plot
        ax11.limits = (nothing, (true_escape[1],true_escape[end]))
        ax11.yticks = [true_escape[1],true_escape[end]]
        =#

        # Loop over the discretised parameter range
        for n in 1:Nμ
                # Extract the ews rvs
                ews = escape[:,n]

                scatter!(ax11, [μ0[n] for m in 1:Ne], ews, color = :black)

                # Compute mean and variance
                ensemble_mean = mean(ews) 

                scatter!(ax11, μ0[n], ensemble_mean, color = :blue, markersize = 15, strokecolor = :black, strokewidth = 1.0)

                #=
                deviation = std(ews) 
                bottom = ensemble_mean - deviation
                top = ensemble_mean + deviation

                # Extract minimum and maximum of the realisations
                low = minimum(ews)
                high = maximum(ews)

                # Plot the range of the rv
                rangebars!(ax11, [μ0[n]+(δμ/2)], [low], [high], color = :black, linewidth = 5.0, whiskerwidth = 30.0)
                # Plot the mean and deviation errorbar
                crossbar!(ax11, μ0[n]+(δμ/2), ensemble_mean, bottom, top, color = (:gray,0.15), width = 0.05, midlinecolor = :black, midlinewidth = 3.0, strokecolor = :black, strokewidth = 3.0)
                =#
        end

        return true_escape
end
