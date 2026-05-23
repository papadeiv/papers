"""
    Simulation script

Storage of the definitions of the system alongside all the settings of the problem.
"""

# System parameters
ε = 0.0                                       # Timescale separation
σ = 1.00                                      # Noise level (additive)
D = (σ^2)/2.0                                 # Diffusion level (additive) 
a = 0.0                                       # Stable equilibrium
α = [10.0, 0.5]                               # Bifurcation parameter 

# Dynamical system  
f(x, α) = α*(a - x)                           # Drift
Λ(t) = ε                                      # Shift/Ramp
η(x) = σ                                      # Diffusion

# Simulation parameters
dt = 1e-1                                     # Timestep
Nt = 1e3                                      # Total number of steps
Ne = 1e3                                      # Number of particles
Nb = 25                                       # Number of bins of the histogram
N_exec = 10                                   # Number of executions of the main algorithm
