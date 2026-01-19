# Import the necessary packages and local modules
include("inc.jl")

# Import the simulation's scripts
include("./scripts/sim.jl")
include("./scripts/proc.jl")
include("./scripts/plot.jl")

# Define the main algorithm
function main()
        # Initial condition 
        g0 = 0.0
        em0 = 0.0
        ew0 = 0.0
        z0 = 14.77
        u0 = [g0, em0, ew0, z0]

        # Build the discrete-time dynamical system
        galor = DeterministicIteratedMap(f, u0, μ)

        # Extract a forward trajectory
        U, t = trajectory(galor, convert(Int64, N))
        display(Matrix(U))
end

# Execute the main
main()
