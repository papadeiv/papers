include("../../../../inc/IO.jl")
include("../../../../inc/SystemAnalysis.jl")

# Specify the settings of the ensemble problem
Nt = convert(Int64, 1e3)
Ne = convert(Int64, 3e3)

# Define the noise level in the system 
σ = 1.000::Float64
η(x) = σ

# Define the normal form
f(x, μ) = -μ - 2*x + 3*(x^2) - (4/5)*(x^3)

# Define the IC
x0 = 0.00::Float64

# Specify the parameter value
μ = 0.000::Float64

# Solve the ensemble problem 
t, u = evolve_ensemble(f, η, μ, IC=x0, δt=5e-3, Nt=Nt, Ne=Ne)

# Export the data 
writeout(t, "../data/figure_05/time.csv")
writeout(u, "../data/figure_05/solutions.csv")

# Execute the postprocessing and plotting scripts
include("../postprocessing/figure_05_postprocessing.jl")
include("../plotting/figure_05_plotting.jl")
