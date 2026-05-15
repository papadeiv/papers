"""
    Simulation script

This is where we store the definition of the system alongside all the settings of the problem.
"""

# System parameters
ε = 0.0                                       # Timescale separation
σ = 0.100                                     # Noise level (additive)
D = (σ^2)/2.0                                 # Diffusion level (additive) 
θ = 1.000                                     # OUP parameter             
μ = 0.000                                     # OUP mean

# Dynamical system 
f(x, μ) = -θ*(x - μ)                          # Drift
Λ(t) = ε                                      # Shift/Ramp
η(x) = σ                                      # Diffusion

# Initial condition 
equilibria = get_equilibria(f, μ, domain=[-10,10])
x0 = [equilibria.stable[1], μ]

# Simulation parameters
dt = 5e-2                                     # Timestep
Nt = 5e4                                      # Total number of steps
Ne = 1e3                                      # Number of particles in the ensemble 
