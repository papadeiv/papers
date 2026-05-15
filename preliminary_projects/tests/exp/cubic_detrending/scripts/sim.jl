"""
    Simulation script

Storage of the definitions of the system alongside all the settings of the problem.
"""

# System parameters
μ0 = -1.00                                    # Initial parameter value
ε = 1e-3                                      # Slow timescale
σ = 0.100                                     # Noise level (additive)
D = (σ^2)/2.0                                 # Diffusion level (additive) 

# Dynamical system  
f(x, μ) = -μ - x^2                            # Drift
Λ(t) = ε                                      # Shift/Ramp
η(x) = σ                                      # Diffusion

# Initial condition 
x0 = [sqrt(-μ0), μ0]

# Time parameters
dt = 5e-2
Nt = 1e4

# Relative width of the sliding window 
window_size = 0.4
