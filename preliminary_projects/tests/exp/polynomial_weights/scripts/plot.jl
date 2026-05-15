"""
    Plotting script

Collection of all the functions used to generate and properly format the figures of the simulations.
"""

# Create empty layouts for the figures
include("./figs.jl")

# Plot the histograms and the linear solution of the least-squares problem 
function plot_linear_solution(bins, pdf, linear_coeff, linear_potential)
        # Plot the histograms
        scatter!(ax1, bins, pdf, color = (CtpRed,0.5), markersize = 20, strokewidth = 2.0)
        scatter!(ax2, bins, pdf, color = (CtpRed,0.5), markersize = 20, strokewidth = 2.0)

        # Plot the ground truth 
        domain = LinRange(-1.5,1.5,1000)
        lines!(ax3, domain, [U(x,μ) for x in domain], color = (CtpRed,1.0), linewidth = 5.0)
        lines!(ax4, domain, [U(x,μ) for x in domain], color = (CtpRed,1.0), linewidth = 5.0)

        # Plot the reconstructed probability density
        domain = LinRange(-0.35,0.35,5000)
        lines!(ax1, domain, [p(x,linear_coeff) for x in domain], color = CtpGreen, linewidth = 5.0)

        # Plot the reconstructed shifted potentials
        domain = LinRange(-1.5,1.5,1000)
        lines!(ax3, domain, [linear_potential(x) for x in domain], color = CtpGreen, linewidth = 5.0)
end

# Plot the nonlinear solution of the least-squares problem 
function plot_nonlinear_solution(nonlinear_coeff, nonlinear_potential, idx)
        # Plot the reconstructed probability density
        domain = LinRange(-0.35,0.35,5000)
        lines!(ax2, domain, [p(x,nonlinear_coeff) for x in domain], color = (CtpBlue, idx/5.0), linewidth = 5.0)

        # Plot the reconstructed shifted potentials
        domain = LinRange(-1.5,1.5,1000)
        lines!(ax4, domain, [nonlinear_potential(x) for x in domain], color = (CtpBlue, idx/5.0), linewidth = 5.0)
end

# Plot the reconstruction error across the ensemble
function plot_error(error_matrix)
        # Loop over the columns of the error matrix
        for (n, error) in enumerate(eachcol(error_matrix)) 
                # Assemble an empirical distribution out of the data
                u_min = error[argmin(error)]
                u_max = error[argmax(error)]
                range = u_max - u_min
                x = LinRange(u_min - range*0.05, u_max + range*0.05, 20)
                bins = [(x[n+1]+x[n])/2 for n in 1:(length(x)-1)]
                hist = StatsBase.fit(Histogram, error, x)
                pdf = (LinearAlgebra.normalize(hist, mode = :pdf)).weights

                # Differentiate between the linear and nonlinear solutions
                if n == 1
                        # Plot the linear ensemble errors
                        #lines!(ax5, bins, pdf, color = CtpGreen, linewidth = 5.0)
                else
                        # Plot the nonlinear ensemble errors
                        lines!(ax5, bins, pdf, color = (CtpBlue, (n-1)/5.0), linewidth = 5.0)
                end
        end

        # Export the figure
        savefig("polynomial_weights/error.png", fig5)
end
