"""
    Main script

Run this file to execute the simulation, analyse and plot the results.
"""

# Import the necessary packages and local modules
include("inc.jl")

# Import the simulation's scripts
include("./scripts/sim.jl")
include("./scripts/plot.jl")

# Define the main algorithm
function main()
        # Generate non-stationary timeseries
        sample_path = evolve(f, η, Λ, x0, steps=Nt, stepsize=dt)
        t = sample_path.time
        μ = sample_path.parameter
        u = (sample_path.state)[1]
        println("μ ∈ [$(μ[1]), $(μ[end])]")
        println("x ∈ [$(u[1]), $(u[end])]")

        # Detrend it 
        PyEMD = pyimport("PyEMD")
        emd = PyEMD.EMD()
        imfs = Array(emd(u))
        display(size(imfs,1))

        # Plot the timeseries decomposition and their distribution
        plot_imfs(t, u, imfs)
end

# Execute the main
main()
