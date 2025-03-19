include("../../../../inc/IO.jl")
include("../../../../inc/SystemAnalysis.jl")

#     Ornstein-Uhlenbeck process
# (linear drift, constant diffusion)

# Specify the parameters of the process
σ = 0.1 
ε = 0.02
α = 0.25
μ = 0.0

# Define the drift and the diffusion of the SDE
f(x, μ) = -(α/ε)*x
g(x) = σ/sqrt(ε)

# Evolve the process from a stable equilibrium as an IC
t, u = evolve_forward(f, g, μ)

# Export the data
writeout(hcat(t, u), "../data/timeseries.csv")
