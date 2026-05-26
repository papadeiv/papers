"""
    Simulation script

Storage of the definitions of the system alongside all the settings of the problem.
"""

# System parameters
μ_set = [-0.3, 0.0, 0.1]                      # Set of bifurcation parameter values
ε = 0.0                                       # Timescale separation
σ = 0.150                                    # Noise level (additive)
D = (σ^2)/2.0                                 # Diffusion level (additive) 

# Polynomial coefficients of the true potential
c0 = 0.0
c2 = -1.0 
c3 = 0.0
c4 = 2.0

# Dynamical system  
f(x, μ) = -(μ + 2*c2*x + 3*c3*x^2 + 4*c4*x^3) # Drift
Λ(t) = ε                                      # Shift/Ramp
η(x) = σ                                      # Diffusion

# Scalar potential of the conservative system 
U(x, μ) = c0 + μ*x + c2*x^2 + c3*x^3 + c4*x^4

# Simulation parameters
dt = 5e-2                                     # Timestep
Nt = 1e5                                      # Total number of steps
