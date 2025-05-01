include("../../../../inc/IO.jl")
include("../../../../inc/SystemAnalysis.jl")

# Define the noise level and timescale separation 
σ = 0.050::Float64
ε = 0.001::Float64

# Define the fast and slow variable vector field
f(x, μ) = -μ - 2*x + 3*(x^2) - (4/5)*(x^3)
η(x) = σ
g(t) = ε

# Define the IC
x0 = 2.41::Float64
μ0 = 1.40::Float64
u0 = [x0, μ0]

# Solve the fast-slow SDE
t, μ, u = evolve_ramped_1d(f, g, η, u0, δt=1e-1, μf=1.7)

# Export the data 
writeout(t, "../data/figure_04/time.csv")
writeout(μ, "../data/figure_04/μ.csv")
writeout(u, "../data/figure_04/u.csv")

# Execute the postprocessing and plotting scripts
include("../postprocessing/figure_04_postprocessing.jl")
include("../plotting/figure_04_plotting.jl")
