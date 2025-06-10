include("../../../../inc/IO.jl")
include("../../../../inc/SystemAnalysis.jl")

# Number of timesteps of each particle
Nt = convert(Int64,1e5)
# Number of (fixed) parameter values
Nμ = convert(Int64,2e0)

# Define the normal form
f(x, μ) = -μ - 2*x + 3*(x^2) -(4/5)*x^3

# Define the noise level in the system 
σ = 0.200::Float64
g(x) = σ

# Specify the parameter range
μ_min = 0.62::Float64
μ_max = 0.63::Float64

# Specify the initial conditions
x0 = 2.73::Float64.*ones(Nμ)

# Propagate the saddle-node normal form forward in time  
t, μ, u = evolve_fixed_1d(f, g, [μ_min, μ_max], x0, Nt=Nt, Nμ=Nμ) 

# Export the data
writeout(hcat(t, u), "../data/solution.csv")
writeout(μ, "../data/parameter.csv")

# Execute the additional scripts for this simulation
include("../postprocessing/02_postprocessing.jl")
include("../plotting/02_plotting.jl")
