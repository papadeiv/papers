include("../../inc/IO.jl")
include("../../inc/SystemAnalysis.jl")

# Define the noise level and timescale separation 
σ = 0.250::Float64
ε = 0.001::Float64

# Define the fast and slow variable vector field
f(x, μ) = -μ - 2*x + 3*(x^2) - (4/5)*(x^3)
g(t) = ε
η(x) = σ

# Define the IC
x0 = 2.52::Float64
μ0 = 1.20::Float64
u0 = [x0, μ0]

# Solve the fast-slow SDE
t, μ, u = evolve_ramped_1d(f, g, η, u0, δt=1e-2, μf=1.43)

# Export the data 
writeout(hcat(t, μ, u), "../data/NLS/solution.csv")

# Execute the postprocessing and plotting scripts
include("../postprocessing/NLS:proc.jl")
include("../plotting/NLS:plot.jl")
