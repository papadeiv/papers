include("../../../../inc/IO.jl")
include("../../../../inc/SystemAnalysis.jl")

# Define the noise level and timescale separation 
σ = 0.200::Float64
ε = 0.100::Float64

# Define the fast and slow variable vector field
f(x, μ) = -μ - (x^2)
g(t) = ε
η(x) = σ

# Define the IC
x0 = 1.00::Float64
μ0 = -(1.00::Float64)
u0 = [x0, μ0]

# Solve the fast-slow SDE
t, μ, u = evolve_ramped_1d(f, g, η, u0, δt=1e-3, μf=-0.500)

# Export the data 
writeout(hcat(t, μ, u), "../data/slide_09/solution.csv")

# Execute the plotting script
include("../plotting/slide_09_plotting.jl")
