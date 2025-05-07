include("../../../../inc/IO.jl")
include("../../../../inc/SystemAnalysis.jl")

#############
# R-tipping #
#############

# Define the noise level in the system 
σ = 0.000::Float64
η(x) = σ

# Define the timescale constant for the parameter ramp
ε = 1.100::Float64
g(t) = ε

# Define the normal form
f(x, μ) = (x + μ)^2 - 1.0

# Define the IC
x0 = -(5.0::Float64)
μ0 = 5.0::Float64
u0 = [x0, μ0]

# Solve the fast-slow SDE
t, μ, u = evolve_ramped_1d(f, g, η, u0, δt=1e-3, Nt=convert(Int64,8e3), saveat=0.01)

# Export the data 
writeout(t, "../data/figure_03/R-tipping/time.csv")
writeout(μ, "../data/figure_03/R-tipping/μ.csv")
writeout(u, "../data/figure_03/R-tipping/u.csv")

#############
# N-tipping #
#############

# Define the noise level in the system 
σ = 0.750::Float64
η(x) = σ

# Define the normal form
f(x, μ) = -μ - 2*x + 3*(x^2) - (4/5)*(x^3)

# Define the IC
x0 = 2.6::Float64
μ0 = 1.0::Float64

# Solve the fast-slow SDE
t, u = evolve_forward_1d(f, η, μ0, IC=x0, δt=1e-3, Nt=convert(Int64,1e5), saveat=0.05)

# Export the data 
writeout(t, "../data/figure_03/N-tipping/time.csv")
writeout(u, "../data/figure_03/N-tipping/u.csv")

# Execute the plotting scripts
include("../plotting/figure_03_plotting.jl")
