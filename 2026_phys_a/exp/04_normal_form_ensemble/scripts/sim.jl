"""
    Simulation script

Storage of the definitions of the system alongside all the settings of the problem.
"""

# System parameters
μ_set = [-0.9, -0.6, -0.3]                    # Bifurcation parameter value
ε = 0.0                                       # Timescale separation
σ = 0.100                                     # Noise level (additive)
D = (σ^2)/2.0                                 # Diffusion level (additive) 

# Dynamical system  
f(x, μ) = -μ - x^2                            # Drift
Λ(t) = ε                                      # Shift/Ramp
η(x) = σ                                      # Diffusion

# Simulation parameters
dt = 1e-1                                     # Timestep
Nt = 1e4                                      # Total number of steps
Ne = 5e3                                      # Number of particles in the ensemble 
Nb = 25                                       # Number of bins of the histogram
