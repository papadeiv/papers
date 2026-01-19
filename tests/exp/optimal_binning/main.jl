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
        # Solve the ensemble problem 
        ensemble = evolve(f, η, Λ, x0, stepsize=dt, steps=Nt, particles=Ne)

        # Loop over the ensemble runs
        @showprogress for n in 1:convert(Int64, Ne)
                # Generate new figures
                include("./scripts/figs.jl")

                # Generate OUP samples
                u = ensemble.state[n]

                # Define the range of the data 
                u_min = u[argmin(u)]
                u_max = u[argmax(u)]
                range = u_max - u_min

                # Define the number of bins
                Nb1 = convert(Int64, ceil(range/(3.49*std(u)*(Nt)^(-1.0/3.0))))   # Scott's rule (1985)
                Nb2 = convert(Int64, ceil(range/(2.0*(quantile(u, 0.75) - quantile(u, 0.25))*(Nt)^(-1.0/3.0))))   # Freedman-Diaconis' rule (1981)
                Nb3 = convert(Int64, ceil(2.0*(Nt)^(1.0/3.0)))   # Rice formula

                # Assemble an empirical distribution of the OUP 
                x = LinRange(u_min - range*0.05, u_max + range*0.05, Nb1)
                bins = [(x[n+1]+x[n])/2 for n in 1:(length(x)-1)]
                hist = StatsBase.fit(Histogram, u, x)
                pdf = (LinearAlgebra.normalize(hist, mode = :pdf)).weights

                # Plot the histograms
                plot_hist(u, bins, pdf, 1)

                # Generate Gaussian drawn samples
                p = Normal(μ, sqrt((σ^2)/(2*θ)))       
                u = rand(p, convert(Int64, Nt))

                # Define the range of the data 
                u_min = u[argmin(u)]
                u_max = u[argmax(u)]
                range = u_max - u_min

                # Define the number of bins
                Nb1 = convert(Int64, ceil(range/(3.49*std(u)*(Nt)^(-1.0/3.0))))   # Scott's rule (1985)
                Nb2 = convert(Int64, ceil(range/(2.0*(quantile(u, 0.75) - quantile(u, 0.25))*(Nt)^(-1.0/3.0))))   # Freedman-Diaconis' rule (1981)
                Nb3 = convert(Int64, ceil(2.0*(Nt)^(1.0/3.0)))   # Rice formula

                # Assemble an empirical distribution of the Gaussian samples
                x = LinRange(u_min - range*0.05, u_max + range*0.05, Nb1)
                bins = [(x[n+1]+x[n])/2 for n in 1:(length(x)-1)]
                hist = StatsBase.fit(Histogram, u, x)
                pdf = (LinearAlgebra.normalize(hist, mode = :pdf)).weights

                # Plot the histograms
                plot_hist(u, bins, pdf, 2)

                # Export the figure
                savefig("optimal_binning/$n.png", fig1)
                savefig("optimal_binning/misfit/$n.png", fig5)
        end
end

# Execute the main
main()
