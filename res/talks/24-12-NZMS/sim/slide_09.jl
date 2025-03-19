include("../../../../inc/IO.jl")
include("../../../../inc/SystemAnalysis.jl")

# Define the noise level in the system 
σ = 0.700::Float64
g(x) = σ

# Define the timescale constant for the parameter ramp
ε = 0.010::Float64
h(t) = ε

# Define the normal form
f(x, μ) = -μ - 2*x + 3*(x^2) - (4/5)*(x^3)

# Define the IC
x0 = 3.33::Float64
μ0 = -(3.0::Float64)
u0 = [x0, μ0]

# Solve the fast-slow SDE
t, μ, u = evolve_ramped_1d(f, g, h, u0, δt=1e-2, μf=4.0)

# Export the data 
writeout(t, "../data/slide_09/time.csv")
writeout(μ, "../data/slide_09/μ.csv")
writeout(u, "../data/slide_09/u.csv")

# Execute the postprocessing and plotting scripts
include("../postprocessing/slide_09_postprocessing.jl")
include("../plotting/slide_09_plotting.jl")
