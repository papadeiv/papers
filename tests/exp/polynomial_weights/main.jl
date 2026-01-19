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
        #bins_data = Matrix{Float64}(undef, Nb, 2)
        #dist_data = Matrix{Float64}(undef, Nb, 2)

        # Loop over the two runs
        for n in 1:2
                for m in 1:convert(Int64, Ne)
                        # Define empty vector
                        errors = Float64[]

                        # Import the solutions of the good and the bad runs
                        run = Matrix(CSV.read("../../res/data/"*runs[n]*"/$m.csv", DataFrame))
                        t = run[:,1] 
                        u = run[:,2] 

                        # Assemble an empirical distribution out of the data
                        u_min = u[argmin(u)]
                        u_max = u[argmax(u)]
                        range = u_max - u_min
                        Nb = convert(Int64, ceil(range/(3.49*std(u)*(Nt)^(-1.0/3.0))))   # Scott's rule (1985)
                        x = LinRange(u_min - range*0.05, u_max + range*0.05, Nb+1)
                        bins = [(x[n+1]+x[n])/2 for n in 1:(length(x)-1)]
                        #bins_data[:,n] = bins
                        hist = StatsBase.fit(Histogram, u, x)
                        pdf = (LinearAlgebra.normalize(hist, mode = :pdf)).weights
                        #dist_data[:,n] = pdf

                        # Determine the median and the mode of the distribution
                        mode_idx = argmax(pdf) 
                        dx = bins[2] - bins[1]
                        cdf = cumsum(pdf)*dx
                        median_idx = findfirst(>=(0.5), cdf) 

                        # Solve the linear least-squares problem
                        idx = findall(x -> x > 0.0, pdf)
                        y = [pdf[n] for n in idx]
                        x = [bins[n] for n in idx]
                        N = 1/sqrt(2*pi*D)
                        y = -D.*log.(y./N)
                        linear_solution = Polynomials.fit(x, y, Nc).coeffs[2:4]

                        # Shift the solutions to center the potential with the ground truth 
                        Vs_linear = shift_potential(U, sqrt(-μ), linear_solution)
                        push!(errors, get_error(Vs_linear))

                        # Plot the linear solution 
                        include("./scripts/figs.jl")
                        plot_linear_solution(bins, pdf, linear_solution, Vs_linear)

                        # Solve the unweighted nonlinear problem starting from the initial guess
                        lower = [-45.0, -25.0, -40.0] 
                        upper = [45.0, 25.0, 40.0]
                        fit = curve_fit(p, bins, pdf, linear_solution, lower=lower, upper=upper, store_trace=true, lambda=1.0, lambda_increase=10.0, lambda_decrease=0.1, maxIter=50)
                        nonlinear_solution = fit.param 

                        # Shift the solutions to center the potential with the ground truth and compute the error 
                        Vs_nonlinear = shift_potential(U, sqrt(-μ), nonlinear_solution)
                        push!(errors, get_error(Vs_nonlinear))

                        # Plot the nonlinear solution 
                        plot_nonlinear_solution(nonlinear_solution, Vs_nonlinear, 1.0)

                        #############################
                        #     Quadratic weights     #
                        #############################
                        
                        # Define the weights of the nonlinear least-squares problem

                        weights = α .+ β.*(bins .- bins[median_idx]).^2

                        # Solve the weighted nonlinear problem starting from the initial guess
                        fit = curve_fit(p, bins, pdf, weights, linear_solution, lower=lower, upper=upper, store_trace=true, lambda=1.0, lambda_increase=10.0, lambda_decrease=0.1, maxIter=50)
                        nonlinear_solution = fit.param 

                        # Shift the solutions to center the potential with the ground truth and compute the error
                        Vs_nonlinear = shift_potential(U, sqrt(-μ), nonlinear_solution)
                        push!(errors, get_error(Vs_nonlinear))

                        # Plot the nonlinear solution 
                        plot_nonlinear_solution(nonlinear_solution, Vs_nonlinear, 2.0)

                        ###########################
                        #     Quartic weights     #
                        ###########################
                        
                        # Define the weights of the nonlinear least-squares problem
                        weights = α .+ β.*(bins .- bins[median_idx]).^4

                        # Solve the weighted nonlinear problem starting from the initial guess
                        fit = curve_fit(p, bins, pdf, weights, linear_solution, lower=lower, upper=upper, store_trace=true, lambda=1.0, lambda_increase=10.0, lambda_decrease=0.1, maxIter=50)
                        nonlinear_solution = fit.param 

                        # Shift the solutions to center the potential with the ground truth and compute the error
                        Vs_nonlinear = shift_potential(U, sqrt(-μ), nonlinear_solution)
                        push!(errors, get_error(Vs_nonlinear))

                        # Plot the nonlinear solution 
                        plot_nonlinear_solution(nonlinear_solution, Vs_nonlinear, 3.0)

                        ##########################
                        #     Sixtic weights     #
                        ##########################
                        
                        # Define the weights of the nonlinear least-squares problem
                        weights = α .+ β.*(bins .- bins[median_idx]).^6

                        # Solve the weighted nonlinear problem starting from the initial guess
                        fit = curve_fit(p, bins, pdf, weights, linear_solution, lower=lower, upper=upper, store_trace=true, lambda=1.0, lambda_increase=10.0, lambda_decrease=0.1, maxIter=50)
                        nonlinear_solution = fit.param 

                        # Shift the solutions to center the potential with the ground truth and compute the error
                        Vs_nonlinear = shift_potential(U, sqrt(-μ), nonlinear_solution)
                        push!(errors, get_error(Vs_nonlinear))

                        # Plot the nonlinear solution 
                        plot_nonlinear_solution(nonlinear_solution, Vs_nonlinear, 4.0)

                        #########################
                        #     Octic weights     #
                        #########################
                        
                        # Define the weights of the nonlinear least-squares problem
                        weights = α .+ β.*(bins .- bins[median_idx]).^8

                        # Solve the weighted nonlinear problem starting from the initial guess
                        fit = curve_fit(p, bins, pdf, weights, linear_solution, lower=lower, upper=upper, store_trace=true, lambda=1.0, lambda_increase=10.0, lambda_decrease=0.1, maxIter=50)
                        nonlinear_solution = fit.param 

                        # Shift the solutions to center the potential with the ground truth and compute the error
                        Vs_nonlinear = shift_potential(U, sqrt(-μ), nonlinear_solution)
                        push!(errors, get_error(Vs_nonlinear))

                        # Plot the nonlinear solution 
                        plot_nonlinear_solution(nonlinear_solution, Vs_nonlinear, 5.0)

                        # Update the error matrix
                        if n == 1
                                error[m,:] = errors
                        else
                                error[(convert(Int64, Ne)+m),:] = errors
                        end

                        # Loop over the trace of the optimiser
                        println("")
                        println("$(runs[n]) run n.$m: skewness = $(skewness(u))")
                        println("")
                        @printf("                  Optimiser's steps\n")
                        @printf("----------------------------------------------------\n")
                        @printf("Iter     Residual       Lambda       Gradient norm\n")
                        @printf("------   ----------     --------     ---------------\n")
                        for m in 1:length(fit.trace)
                                @printf("%3d    %11.6f    %11.6f    %13.8f\n", fit.trace[m].iteration, fit.trace[m].value, fit.trace[m].metadata["lambda"], fit.trace[m].g_norm)
                        end
                        println("")
                        println("guess = [", (linear_solution)[1], ", ", (linear_solution)[2], ", ", (linear_solution)[3], "]")
                        println("final = [", (nonlinear_solution)[1], ", ", (nonlinear_solution)[2], ", ", (nonlinear_solution)[3], "]")
                        println("")

                        # Export the figure
                        savefig("polynomial_weights/"*runs[n]*"/$m.png", fig1)
                end
        end

        # Plot the reconstruction error ensemble distribution
        plot_error(error)
end

# Execute the main
main()
