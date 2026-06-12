"""
    Simulation script

Storage of the definitions of the system alongside all the settings of the problem.
"""

# System parameters
μ_set = collect(range(-0.9, stop=-0.3, step=0.1)) # Bifurcation parameter value
α_set = LinRange(0, 1e-2, 1000)                   # Regularisation coefficient
ε = 0.0                                           # Timescale separation
σ = 0.100                                         # Noise level (additive)
D = (σ^2)/2.0                                     # Diffusion level (additive) 

# Dynamical system  
f(x, μ) = -μ - x^2                                # Drift
Λ(t) = ε                                          # Shift/Ramp
η(x) = σ                                          # Diffusion

# Simulation parameters
dt = 1e-1                                         # Timestep
Nt = 1e3                                          # Total number of steps
Ne = 2e1                                          # Number of particles in the ensemble 
