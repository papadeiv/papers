include("../../../../inc/IO.jl")
include("../../../../inc/SystemAnalysis.jl")

# Define the noise level and timescale separation 
σ = 0.250::Float64
ε = 0.010::Float64

# Define the fast and slow variable vector field
f(x, μ) = -μ - 2*x + 3*(x^2) - (4/5)*(x^3)
g(t) = ε
η(x) = σ

# Define the IC
x0 = 2.88::Float64
μ0 = 0.00::Float64
u0 = [x0, μ0]

# Solve the fast-slow SDE
t, μ, u = evolve_ramped_1d(f, g, η, u0, δt=1e-1, μf=1.800)

# Export the data 
writeout(hcat(t, μ, u), "../data/slide_08/solution.csv")

# Execute the postprocessing and plotting scripts
include("../postprocessing/slide_08_postprocessing.jl")
include("../plotting/slide_08_plotting.jl")
