"""
    Simulation script

Storage of the definitions of the system alongside all the settings of the problem.
"""

# System parameters
μ0 = -1.0                                         # Initial value of the bifurcation parameter
μf = 0.0                                          # Final value of the bifurcation parameter
ε = 9e-5#1e-2                                     # Timescale separation
σ = 0.100                                         # Noise level (additive)
D = (σ^2)/2.0                                     # Diffusion level (additive) 

# Dynamical system  
f(x, μ) = -μ - x^2                                # Drift
Λ(t) = ε                                          # Shift/Ramp
η(x) = σ                                          # Diffusion

# Simulation parameters
dt = 1e-1#1e-3                                    # Timestep
Ne = 1e0                                          # Number of particles in the ensemble 
