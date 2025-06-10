include("../../../../inc/IO.jl")
include("../../../../inc/SystemAnalysis.jl")

#############
# B-tipping #
#############

# Define the noise level in the system 
σ = 0.100::Float64
η(x) = σ

# Define the timescale constant for the parameter ramp
ε = 0.100::Float64
g(t) = ε

# Define the normal form
f(x, μ) = -μ - 2*x + 3*(x^2) - (4/5)*(x^3)

# Define the IC
x0 = 3.33::Float64
μ0 = -(3.0::Float64)
u0 = [x0, μ0]

# Solve the fast-slow SDE
t, μ, u = evolve_ramped_1d(f, g, η, u0, δt=1e-3, saveat=0.10, μf=4.75)

# Export the data 
writeout(hcat(t,μ,u), "../data/slide_04/B-tipping.csv")

#############
# N-tipping #
#############

# Define the noise level in the system 
σ = 2.000::Float64
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
writeout(hcat(t,u), "../data/slide_04/N-tipping.csv")

# Execute the postprocessing and plotting scripts
include("../postprocessing/slide_04_postprocessing.jl")
include("../plotting/slide_04_plotting.jl")
