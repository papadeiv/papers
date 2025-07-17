include("../../../inc/IO.jl")
include("../../../inc/SystemAnalysis.jl")

# Define the dynamical system (diffusion) 
η(x) = σ
σ = 0.250::Float64

# Define the dynamical system (drift) 
f(x, μ) = -μ - 2*x + 3*(x^2) - (4/5)*(x^3)

# Define the IC
x0 = 2.60::Float64
μ0 = 1.00::Float64

# Solve the SDE
Nt = convert(Int64, 1e6)
t, u = evolve_forward_1d(f, η, μ0, IC=x0, δt=1e-3, Nt=Nt)

# Export the data 
writeout(hcat(t, u), "../data/MLE/solution.csv")

# Execute the postprocessing and plotting scripts
include("../proc/MLE_potential:proc.jl")
#include("../plot/MLE_potential:plot.jl")
