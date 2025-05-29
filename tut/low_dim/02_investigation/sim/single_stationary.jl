include("../../../../inc/IO.jl")
include("../../../../inc/SystemAnalysis.jl")

# Number of timesteps for the SDE 
Nt = convert(Int64,1e5)
# Number of (fixed) parameter values
Nμ = convert(Int64,1e3)

# Define the normal form
f(x, μ) = -μ - 2*x + 3*(x^2) -(4/5)*x^3

# Define the noise level in the system 
σ = 0.200::Float64
g(x) = σ

# Specify the parameter range
μ_min = 0.400::Float64
μ_max = 1.450::Float64

# Specify the initial conditions
x0 = [sort(get_equilibria(f, μ, domain=[-20,20])[1], rev=true)[1] for μ in LinRange(μ_min, μ_max, Nμ)] 

# Propagate the saddle-node normal form forward in time  
t, μ, u = evolve_fixed_1d(f, g, [μ_min, μ_max], x0, Nt=Nt, Nμ=Nμ)

# Export the data
writeout(hcat(μ, x0), "../data/equilibria.csv")
writeout(u, "../data/solutions.csv")

# Execute the additional scripts for this simulation
include("../postprocessing/single_stationary_postprocessing.jl")
include("../plotting/single_stationary_plotting.jl")
