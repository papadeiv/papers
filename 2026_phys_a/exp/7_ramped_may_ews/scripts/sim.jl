"""
    Simulation script

Storage of the definitions of the system alongside all the settings of the problem.
"""

# System parameters
r = 1.0                                           # Growth rate
k = 10.0                                          # Carrying capacity
h = 1.0                                           # Half-grazing biomass
μ0 = 1.80                                         # Initial value of the bifurcation parameter
μf = 2.60                                         # Final value of the bifurcation parameter
ε = 5e-5                                          # Timescale separation
σ = 0.100                                         # Noise level (additive)
D = (σ^2)/2.0                                     # Diffusion level (additive) 

# Dynamical system  
f(x, μ) = r*x*(1-x/k)-μ*((x^2)/((x^2)+(h^2)))     # Drift
Λ(t) = ε                                          # Shift/Ramp
η(x) = σ                                          # Diffusion

# Simulation parameters
dt = 1e-1                                         # Timestep
Ne = 1e2                                          # Number of particles in the ensemble 
