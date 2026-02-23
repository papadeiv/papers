"""
    Simulation script

Storage of the definitions of the system alongside all the settings of the problem.
"""

# System parameters
M0 = collect(range(0.0, stop=1.3, step=0.1))  # Collection of initial parameter values
ε = 1e-6                                      # Slow timescale
σ = 0.100::Float64                            # Noise level (additive)
D = (σ^2)/2.0                                 # Diffusion level (additive) 

# Dynamical system  
f(x, μ) = -μ + 4*x - 4*x^3                    # Drift
Λ(t) = ε                                      # Shift/Ramp
η(x) = σ                                      # Diffusion

# Simulation parameters
δt = 5e-2                                     # Timestep
Nt = 2e5                                      # Total number of steps
Ne = 2e2                                      # Number of particles in the ensemble 
