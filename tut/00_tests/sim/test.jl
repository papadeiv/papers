include("../../../inc/IO.jl")
include("../../../inc/SystemAnalysis.jl")

#=
# Define the noise level and timescale separation 
σ = 0.125::Float64
ε = 0.001::Float64

# Define the fast and slow variable vector field
f(x, μ) = -μ - 2*x + 3*(x^2) - (4/5)*(x^3)
g(t) = ε
η(x) = σ

# Define the IC
x0 = 2.60::Float64
μ0 = 1.00::Float64
u0 = [x0, μ0]

# Solve the fast-slow SDE
t, μ, u = evolve_ramped_1d(f, g, η, u0, δt=5e-2, μf=1.675)

# Export the data 
writeout(hcat(t, μ, u), "../data/sim_03/solution.csv")
=#

# Execute the postprocessing and plotting scripts
include("../postprocessing/sim_03_postprocessing.jl")
include("../plotting/sim_03_plotting.jl")
