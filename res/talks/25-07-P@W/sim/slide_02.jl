include("../../../../inc/IO.jl")
include("../../../../inc/SystemAnalysis.jl")

#=
# Define the diffusion 
σ = 0.50::Float64
η(x) = σ

# Define the drift 
f(x, μ) = -μ - 2*x + 3*(x^2) - (4/5)*(x^3)

# Define the IC
x0 = -(0.07249::Float64)
μ0 = 0.6::Float64

# Define time simulation settings
δt = 5e-3
Nt = convert(Int64,5e6)

# Solve the fast-slow SDE
t, u = evolve_forward_1d(f, η, μ0, IC=x0, Nt=Nt, δt=δt)

# Export the data 
writeout(hcat(t, u), "../data/slide_02/solution.csv")

# Execute the postprocessing and plotting scripts
=#
include("../postprocessing/slide_02_postprocessing.jl")
include("../plotting/slide_02_plotting.jl")
