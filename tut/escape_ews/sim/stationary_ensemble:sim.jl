include("../../../inc/IO.jl")
include("../../../inc/SystemAnalysis.jl")

# Define the dynamical system (drift) 
f(x, μ) = -μ - 2*x + 3*(x^2) - (4/5)*(x^3)

# Define the dynamical system (diffusion)
σ = 0.125::Float64
η(x) = σ

# Define the IC and parameter value
x0 = 2.52::Float64
μ0 = 1.20::Float64

# Define simulation parameters
δt = 1e-3
Nt = convert(Int64, 5e3)
Ne = convert(Int64, 5e3)

# Solve the ensemble problem 
t, u = evolve_ensemble(f, η, μ0, IC=x0, δt=δt, Nt=Nt, Ne=Ne)

# Export the data 
writeout(t, "../data/stationary_ensemble/time.csv")
writeout(u, "../data/stationary_ensemble/solutions.csv")

# Execute the postprocessing and plotting scripts
include("../proc/stationary_ensemble:proc.jl")
include("../plot/stationary_ensemble:plot.jl")
