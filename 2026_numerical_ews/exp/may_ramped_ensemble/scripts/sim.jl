"""
    Simulation script

Storage of the definitions of the system alongside all the settings of the problem.
"""

# System parameters
r = 1.0                                       # Growth rate
k = 10.0                                      # Carrying capacity
Vh = 1.0                                      # Half-grazing biomass
M0 = collect(range(1.8, stop=2.5, step=0.1))  # Collection of initial bifurcation parameters (grazing rate) 
ε = 1e-6                                      # Slow timescale
σ = 0.100::Float64                            # Noise level (additive)
D = (σ^2)/2.0                                 # Diffusion level (additive) 

# Dynamical system  
f(x, μ) = r*x*(1.0-(1.0/k)*x) - μ*((x^2)/((x^2) + (Vh^2)))  # Drift
Λ(t) = ε                                                    # Shift/Ramp
η(x) = σ                                                    # Diffusion

# Simulation parameters
δt = 5e-2                                     # Timestep
Nt = 5e5                                      # Total number of steps
Ne = 2e2                                      # Number of particles in the ensemble 
