include("../../../../inc/IO.jl")
include("../../../../inc/SystemAnalysis.jl")

# Payoff matrix
A = [1.0, 0.0, 0.0, 1.0]

# Control matrix
G = [1.0, 0.0, 0.0, 1.0]

# Define a set of initial conditions (ICs) 
Nx = 100
ICs = collect(LinRange(0.05, 0.95, Nx))

# Define time parameters of the simulation
Nt = 1000
dt = 1e-2

# Define the dynamical system (deterministic drift)
f(x, μ, t) = x*(1 - x)*((μ[1] + μ[4] - μ[2] - μ[3])*x + μ[2] - μ[4])

# Execute the postprocessing and plotting script
include("../plotting/figure_02_postprocessing.jl")
include("../plotting/figure_02_plotting.jl")
