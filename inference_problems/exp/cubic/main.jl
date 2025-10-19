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
        for n in 1:2
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
                #linear_solution = [-0.00372292,0.56030212,0.37830491] # Good guess 1
                #linear_solution = [-0.00039618,0.58196175,0.31332258] # Good guess 2

                # Use the above as an initial guess for the nonlinear regression
                lower = [-45.0, -25.0, -40.0] 
                upper = [45.0, 25.0, 40.0]
                fit = curve_fit(p, bins, pdf, linear_solution, lower=lower, upper=upper, store_trace=true)
                nonlinear_solution = fit.param 

                # Loop over the trace of the optimiser
                println("")
                @printf("                  Optimiser's steps\n")
                @printf("----------------------------------------------------\n")
                @printf("Iter     Residual       Lambda       Gradient norm\n")
                @printf("------   ----------     --------     ---------------\n")
                for m in 1:length(fit.trace)
                        @printf("%3d    %11.6f    %11.6f    %13.8f\n", fit.trace[m].iteration, fit.trace[m].value, fit.trace[m].metadata["lambda"], fit.trace[m].g_norm)
                end
                println("")
                #println("guess = [", (linear_solution)[1], ", ", (linear_solution)[2], ", ", (linear_solution)[3], "]")
                #println("final = [", (nonlinear_solution)[1], ", ", (nonlinear_solution)[2], ", ", (nonlinear_solution)[3], "]")

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
