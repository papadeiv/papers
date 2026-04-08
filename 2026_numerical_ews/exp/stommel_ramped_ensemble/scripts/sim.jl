"""
    Simulation script

Storage of the definitions of the system alongside all the settings of the problem.
"""

# System parameters
η1 = 3.0                                           # Temperature meridional difference 
M0 = collect(range(0.301, stop=0.751, step=0.05))  # Collection of initial bifurcation parameters (η2) 
η3 = 0.1                                           # Timescale ratio between temperature and salinity 
ε = 1e-6                                           # Slow timescale
σ = 0.100::Float64                                 # Noise level (additive)
D = (σ^2)/2.0                                      # Diffusion level (additive) 

# Dynamical system: drift 
f1(x, y, μ) = η1 - x - abs(x - y)*x
f2(x, y, μ) = μ - η3*y - abs(x - y)*y
f = [f1, f2]

# Dynamical system: shift
Λ(t) = ε

# Dynamical system: diffusion
g(x, y) = σ
η = [g, g]

# Simulation parameters
δt = 5e-2                                     # Timestep
Nt = 5e5                                      # Total number of steps
Ne = 2e2                                      # Number of particles in the ensemble 
