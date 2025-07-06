include("../../../../inc/IO.jl")
include("../../../../inc/SystemAnalysis.jl")

#=
# Define the dynamical system (drift)
f(x, μ) = -μ - 2*x + 3*(x^2) - (4/5)*(x^3)

# Define the dynamical system (diffusion)
σ = 0.150::Float64
η(x) = σ

# Define the dynamical system (shift)
ε = 0.001::Float64
Λ(t) = ε

# Define the IC
x0 = 2.60::Float64
μ0 = 1.00::Float64
u0 = [x0, μ0]

# Solve the fast-slow SDE
t, μ, u = evolve_ramped_1d(f, Λ, η, u0, δt=5e-2, μf=1.675)

# Export the data 
writeout(hcat(t, μ, u), "../data/slide_04/solution.csv")
=#

# Execute the postprocessing and plotting scripts
#include("../postprocessing/slide_04_postprocessing.jl")
include("../plotting/slide_04_plotting.jl")
