include("../../../../inc/IO.jl")
include("../../../../inc/SystemAnalysis.jl")

# Define the stationary parameters of the process
σ = 0.1        # Level of (additive) noise
ε = 0.01       # Timescale of the slow ramping of the forcing parameter
r = 1.0        # Growth rate
k = 10.0       # Carrying capacity
Vh = 1.0       # Half-grazing biomass

# Define the IC for the dynamic variables (state and forcing parameter)
V0 = 8.0       # Vegetation biomass (state)
c0 = 1.5       # Grazing rate (forcing parameter)
x0 = [V0, c0]

# Define the forcing parameter endstate
cT = 3.52

# Define the vegetation dynamics model (May 1977) 
f(x, μ) = r*x*(1.0-(1.0/k)*x) - μ*((x^2)/((x^2) + (Vh^2)))
g(x) = ε
η(x) = σ

# Solve the fast-slow SDE
t, μ, x = evolve_ramped_1d(f, g, η, x0, δt = 1e-3, μf = cT)

# Export the data 
writeout(hcat(t, μ, x), "../data/slide_03/May77.csv")

# Execute the plotting script
include("../plotting/slide_03_1_plotting.jl")
