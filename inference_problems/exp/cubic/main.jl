"""
    ???

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
        # Solve the 1-dimensional ensemble problem
        t, u = evolve_1d(f, η, μ, u0, δt=δt, Nt=Nt, Ne=Ne)

        display(u)
end

# Execute the main
main()
