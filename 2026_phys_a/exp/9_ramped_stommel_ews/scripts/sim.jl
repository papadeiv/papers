"""
    Simulation script

Storage of the definitions of the system alongside all the settings of the problem.
"""

# System parameters
η1 = 3.00                                         # Temperature meridional difference 
η3 = 0.30                                         # Timescale ratio between temperature and salinity 
μ0 = 0.00                                         # Initial value of the bifurcation parameter
μf = 1.20                                         # Final value of the bifurcation parameter
ε = 5e-5                                          # Timescale separation
σ = 0.050                                         # Noise level (additive)
D = (σ^2)/2.0                                     # Diffusion level (additive) 

# Dynamical system  
f1(x, y, μ) = η1 - x - abs(x - y)*x               # Temperature (drift)
f2(x, y, μ) = μ - η3*y - abs(x - y)*y             # Salinity (drift)
f = [f1, f2]                                      # Drift
Λ(t) = ε                                          # Shift/Ramp
g(x, y) = σ                                       # Noise level (additive) 
η = [g, g]                                        # Diffusion

# Simulation parameters
dt = 1e-1                                         # Timestep
Ne = 1e2                                          # Number of particles in the ensemble 
