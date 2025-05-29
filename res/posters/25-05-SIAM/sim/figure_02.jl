include("../../../../inc/IO.jl")
include("../../../../inc/SystemAnalysis.jl")

# Define the noise level in the system 
σ = 0.300::Float64
η(x) = σ

# Define the timescale constant for the parameter ramp
ε = 0.010::Float64
g(t) = ε

# Define the normal form
f(x, μ) = -μ - 2*x + 3*(x^2) - (4/5)*(x^3)

# Define the IC
x0 = 3.33::Float64
μ0 = -(3.0::Float64)
u0 = [x0, μ0]

# Solve the fast-slow SDE
t, μ, u = evolve_ramped_1d(f, g, η, u0, Nt=convert(Int64,1e3), μf=4.0)

# Export the data 
writeout(t, "../data/figure_02/time.csv")
writeout(μ, "../data/figure_02/μ.csv")
writeout(u, "../data/figure_02/u.csv")

# Execute the postprocessing and plotting scripts
include("../postprocessing/figure_02_postprocessing.jl")
include("../plotting/figure_02_plotting.jl")
