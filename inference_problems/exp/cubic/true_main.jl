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
        # Import the solutions of the good and the bad runs
        run = Matrix(CSV.read("../../res/data/bad.csv", DataFrame))

        # Loop over the ensemble's sample paths
        for n in 1:Na
                # Extract the trajectory until the tipping and center it 
                t = run[:,1] 
                u = run[:,2] 

                # Solve the nonlinear least-squares problem to fit a local cubic potential
                solution = fit_potential(u, n_bins=Nb, noise=σ, optimiser=β)
                println("guess = [", (solution.guess)[1], ", ", (solution.guess)[2], ", ", (solution.guess)[3], "]")
                println("final = [", (solution.fit)[1], ", ", (solution.fit)[2], ", ", (solution.fit)[3], "]")

                # Perform postprocessing analysis on the solutions
                Vs_linear = shift_potential(U, sqrt(-μ), solution.guess)
                Vs_nonlinear = shift_potential(U, sqrt(-μ), solution.fit)

                # Plot and export the potential reconstruction 
                include("./scripts/figs.jl")
                plot_solutions(t, u, Vs_linear, Vs_nonlinear, n)
        end
end

# Execute the main
main()
