include("../../../../inc/IO.jl")
include("../../../../inc/SystemAnalysis.jl")

# Specify the settings of the ensemble problem
Nt = convert(Int64, 1e4)
Ne = convert(Int64, 3e3)
Nμ = convert(Int64, 2e2)

# Define the noise level in the system 
σ = 0.200::Float64
η(x) = σ

# Define the normal form
f(x, μ) = -μ - 2*x + 3*(x^2) - (4/5)*(x^3)

# Define the IC
x0 = 0.00::Float64

# Specify the parameter range
μ_min = 0.400::Float64
μ_max = 1.250::Float64

# Specify the initial conditions
x0 = [sort(get_equilibria(f, μ, domain=[-20,20])[1], rev=true)[1] for μ in LinRange(μ_min, μ_max, Nμ)] 

# Propagate the saddle-node normal form forward in time  
t, μ, u = evolve_ensemble_fixed(f, η, [μ_min, μ_max], x0, Nt=Nt, Nμ=Nμ, Ne=Ne)

# Export the data 
writeout(hcat(μ, x0), "../data/figure_06/equilibria.csv")
for n in 1:Nμ
        writeout(u[n], "../data/figure_06/solutions/$n.csv")
end

# Execute the postprocessing and plotting scripts
include("../postprocessing/figure_06_postprocessing.jl")
include("../plotting/figure_06_plotting.jl")
