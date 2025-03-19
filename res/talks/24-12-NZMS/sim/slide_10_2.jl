include("../../../../inc/IO.jl")
include("../../../../inc/SystemAnalysis.jl")

# Number of timesteps of each particle
Nt = convert(Int64,1e4)
# Number of (fixed) parameter values
Nμ = convert(Int64,1e2)
# Number of particles in the ensemble
Ne = convert(Int64,2e3)

# Define the noise level in the system 
σ = 0.250::Float64
g(x) = σ

# Define the normal form
f(x, μ) = -μ - 2*x + 3*(x^2) - (4/5)*(x^3)

# Define the parameter range
μ_range  = LinRange(-2.00, 1.25, Nμ)

# Define array to store the ICs
u0 = Vector{Float64}(undef, Nμ)
μ0 = Vector{Float64}(undef, Nμ)

# Define the ICs for each parameter
for n in 1:Nμ
        μ0[n] = μ_range[n]
        equilibria = get_equilibria(f, μ_range[n], domain=[1,4])
        stable = equilibria[1]
        u0[n] = stable[1]
end

# Solve the ensemble problem 
t, p, u = evolve_ensemble_fixed(f, g, μ0, u0, Nt=Nt, Ne=Ne, δt=1e-1)

# Export the data 
writeout(μ0, "../data/slide_10_2/μ.csv")

# Loop over the parameter values
for n in 1:length(u)
        writeout(u[n], "../data/slide_10_2/ensemble/$n.csv")
end

# Execute the postprocessing and plotting scripts
include("../postprocessing/slide_10_2_postprocessing.jl")
include("../plotting/slide_10_2_plotting.jl")
