"""
        ?

???
"""

# Import the necessary packages and local modules
include("inc.jl")

# Import the simulation's scripts
include("./scripts/sim.jl")
include("./scripts/proc.jl")
include("./scripts/plot.jl")

# Define the main algorithm
function main()
        # Define the two runs
        runs = ["bad", "good"]

        # Loop over the two runs
        for n in 1:1
                # Import the solutions of the good and the bad runs
                run = Matrix(CSV.read("../../res/data/"*runs[n]*".csv", DataFrame))
                t = run[:,1] 
                u = run[:,2] 

                # Assemble an empirical distribution out of the data
                u_min = u[argmin(u)]
                u_max = u[argmax(u)]
                range = u_max - u_min
                x = LinRange(u_min - range*0.05, u_max + range*0.05, Nb+1)
                bins = [(x[n+1]+x[n])/2 for n in 1:(length(x)-1)]
                hist = StatsBase.fit(Histogram, u, x)
                pdf = (LinearAlgebra.normalize(hist, mode = :pdf)).weights

                # Solve the linear least-squares problem
                idx = findall(x -> x > 0.0, pdf)
                y = [pdf[n] for n in idx]
                x = [bins[n] for n in idx]
                N=1/sqrt(2*pi*D)
                y = -D.*log.(y./N)
                linear_solution = Polynomials.fit(x, y, Nc).coeffs[2:4]

                # Define the nonlinear model
                model(x, μ) = normalise(f, μ)*f(x, μ)
                #jacob(x, μ) = [(1.0/x)*exp(-μ[2]/x), (-μ[1]/(x^2))*exp(-μ[2]/x)] 
                
                # Compute the linear solution and use it as an initial guess
                guess = linear_solution
               
                # Plot the Jacobian approximation 
                fig = Figure(; size = (1500, 500))
                ax1 = Axis(fig[1,1], xlabel=L"x", ylabel=L"\partial_{c_1}\;f")
                #lines!(ax1, xdata, [jacob(x, μ)[1] for x in xdata], color=:red, label="true", linewidth = 3.0)
                lines!(ax1, bins, approx_jacobian(model, bins, guess)[:,1], color=:orange, label="approximation")
                ax2 = Axis(fig[1,2], xlabel=L"x", ylabel=L"\partial_{c_2}\;f")
                #lines!(ax2, xdata, [jacob(x, μ)[2] for x in xdata], color=:blue, label="true", linewidth = 3.0)
                lines!(ax2, bins, approx_jacobian(model, bins, guess)[:,2], color=:orange, label="approximation")
                ax3 = Axis(fig[1,3], xlabel=L"x", ylabel=L"\partial_{c_3}\;f")
                #lines!(ax3, xdata, [jacob(x, μ)[3] for x in xdata], color=:blue, label="true", linewidth = 3.0)
                lines!(ax3, bins, approx_jacobian(model, bins, guess)[:,3], color=:orange, label="approximation")
                #axislegend(ax1)
                #axislegend(ax2)
                #axislegend(ax3)
                savefig("cubic/jacobian/bad/0.png", fig)

                # Fit the data using gradient descent
                nonlinear_solution = fit_gradient_descent(bins, pdf, model, guess, jac=nothing, maxiter=10000, α=1e-6, tol=1e-10)

                # Shift the solutions to center the potential with the ground truth 
                Vs_linear = shift_potential(U, sqrt(-μ), linear_solution)
                Vs_nonlinear = shift_potential(U, sqrt(-μ), nonlinear_solution)

                # Plot and export the solutions 
                include("./scripts/figs.jl")
                plot_solutions(bins, pdf, linear_solution, nonlinear_solution, Vs_linear, Vs_nonlinear, runs[n])
        end
end

# Execute the main
main()
