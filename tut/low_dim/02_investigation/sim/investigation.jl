include("../../../../inc/IO.jl")
include("../../../../inc/SystemAnalysis.jl")

# Number of timesteps for the SDE 
Nt = convert(Int64,1e5)

# Define the normal form
f(x, μ) = -μ - 2*x + 3*(x^2) -(4/5)*x^3

# Define the noise level in the system 
σ = 0.200::Float64
g(x) = σ

# Specify the parameter range
μ_1 = 0.6617117117117117::Float64
μ_2 = 0.9244744744744745::Float64
μ_3 = 1.1872372372372372::Float64

# Specify the initial conditions
x0 = [sort(get_equilibria(f, μ, domain=[-20,20])[1], rev=true)[1] for μ in [μ_1, μ_2, μ_3]] 

# Propagate the saddle-node normal form forward in time  
t, μ, u = evolve_fixed_1d(f, g, [μ_1, μ_2, μ_3], x0, Nt=Nt)

# Export the data
writeout(hcat(μ, x0), "../data/equilibria.csv")
writeout(u, "../data/solutions.csv")

# Execute the additional scripts for this simulation
include("../postprocessing/investigation_postprocessing.jl")
include("../plotting/investigation_plotting.jl")
