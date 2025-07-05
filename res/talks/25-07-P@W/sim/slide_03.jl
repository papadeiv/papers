include("../../../../inc/IO.jl")
include("../../../../inc/SystemAnalysis.jl")

#=
# Define the diffusion
σ = 0.4::Float64
η(x) = σ

# Define the drift 
f(x, μ) = -μ - 2*x + 3*(x^2) - (4/5)*(x^3)

# Define the shift 
ε = 1e-3
g(t) = ε

# Define the IC
x0 = 2.73::Float64
μ0 = 0.6::Float64
u0 = [x0, μ0]

# Define time simulation settings
δt = 5e-3
μf = 1.6::Float64

# Solve the fast-slow SDE
t, μ, u = evolve_ramped_1d(f, g, η, u0, δt=δt, μf=μf)

# Export the data 
writeout(hcat(t, μ, u), "../data/slide_03/solution.csv")

# Execute the postprocessing and plotting scripts
include("../postprocessing/slide_03_postprocessing.jl")
=#
include("../plotting/slide_03_plotting.jl")
