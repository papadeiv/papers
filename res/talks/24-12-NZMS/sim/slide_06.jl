include("../../../../inc/IO.jl")
include("../../../../inc/SystemAnalysis.jl")

# Define the noise level and timescale separation 
σ = 0.1000::Float64
ε = 0.0001::Float64

# Define the fast and slow variable vector field
f(x, μ) = -μ - 2*x + 3*(x^2) - (4/5)*(x^3)
g(x) = σ
h(t) = ε

# Define the IC
x0 = 2.41::Float64
μ0 = 1.40::Float64
u0 = [x0, μ0]

# Solve the fast-slow SDE
t, μ, u = evolve_ramped_1d(f, g, h, u0, δt=1e-1, μf=1.7)

# Export the data 
writeout(t, "../data/slide_06/time.csv")
writeout(μ, "../data/slide_06/μ.csv")
writeout(u, "../data/slide_06/u.csv")

# Execute the postprocessing and plotting scripts
include("../postprocessing/slide_06_postprocessing.jl")
include("../plotting/slide_06_plotting.jl")
