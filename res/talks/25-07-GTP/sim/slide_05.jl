include("../../../../inc/IO.jl")
include("../../../../inc/SystemAnalysis.jl")

# Define the noise level in the system 
σ = 0.125::Float64
g(x) = σ

# Define the normal form
f(x, μ) = -μ - 2*x + 3*(x^2) - (4/5)*(x^3)

# Define the IC
x0 = 2.730::Float64
μ0 = 0.625::Float64 # Maxwell point

# Solve the fast-slow SDE
printstyled("Solving the autonomous SDE\n"; bold=true, underline=true, color=:light_blue)
t, u = evolve_forward_1d(f, g, μ0, IC=x0, δt=1e-3, Nt=convert(Int64,5e4), saveat=0.001)

# Export the data 
writeout(hcat(t,u), "../data/slide_05/solution.csv")

# Execute the postprocessing and plotting scripts
include("../postprocessing/slide_05_postprocessing.jl")
include("../plotting/slide_05_plotting.jl")
