"""
    Simulation script

Storage of the definitions of the system alongside all the settings of the problem.
"""

# Bifurcation parameter
μ_set = collect(range(start=0.0, step=0.25, stop=2.0)) 

# Dynamical system
f(x, μ) = -μ + 4*x - 4*x^3                    # Drift
Λ(t) = 0.0                                    # Shift/Ramp
η(x) = 0.0                                    # Diffusion

# Simulation parameters
x0 = 0.5                                      # Initial condition
dt = 1e-2                                     # Timestep
Nt = 6e2                                      # Total number of steps
