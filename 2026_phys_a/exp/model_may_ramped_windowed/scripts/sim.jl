"""
    Simulation script

Storage of the definitions of the system alongside all the settings of the problem.
"""

# System parameters
r = 1.0                                       # Growth rate
k = 10.0                                      # Carrying capacity
h = 1.0                                       # Half-grazing biomass
c0 = 1.00                                     # Initial value of the bifurcation parameter 
cf = 3.0                                      # Final value of the bifurcation parameter
ε = 5e-5                                      # Slow timescale
σ = 0.300                                     # Noise level (additive)
D = (σ^2)/2.0                                 # Diffusion level (additive) 

# Dynamical system  
f(x, μ) = r*x*(1 - x/k) - μ*((x^2)/((x^2) + (h^2)))         # Drift
Λ(t) = ε                                                    # Shift/Ramp
η(x) = x > 2.0 ? σ  : σ/10.0                                # Diffusion

# Simulation parameters
δt = 1e-1                                     # Timestep
Ne = 1e2                                      # Number of particles in the ensemble 
