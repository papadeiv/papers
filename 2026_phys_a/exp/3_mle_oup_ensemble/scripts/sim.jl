"""
    Simulation script

Storage of the definitions of the system alongside all the settings of the problem.
"""

# System parameters
ε = 0.0                                       # Timescale separation
σ = 0.100                                     # Noise level (additive)
D = (σ^2)/2.0                                 # Diffusion level (additive) 
a = 0.0                                       # Stable equilibrium
α = -1.0                                       # Bifurcation parameter 

# Dynamical system  
f(x, α) = α*(x - a)                           # Drift
Λ(t) = ε                                      # Shift/Ramp
η(x) = σ                                      # Diffusion

# Simulation parameters
dt = 1e-3                                     # Timestep
Nt = 2e5                                      # Total number of steps
Ne = 1e3                                      # Number of particles
Nb = 50                                       # Number of bins of the histogram
N_exec = 20                                   # Number of executions of the main algorithm
