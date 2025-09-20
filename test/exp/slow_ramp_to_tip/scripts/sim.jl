"""
    Simulation script

This is where we store the definition of the system alongside all the settings of the problem.
"""

# System parameters
μ0 = 1.30::Float64                            # Initial parameter value
μf = 1.40::Float64                            # Final parameter value
ε = 1e-3                                      # Slow timescale
σ = 0.250::Float64                            # Noise level (additive)
D = (σ^2)/2.0                                 # Diffusion level (additive) 

# Dynamical system  
f(x, μ) = -μ - 2*x + 3*(x^2) - (4/5)*(x^3)    # Drift
g(t) = ε                                      # Shift/Ramp
η(x) = σ                                      # Diffusion

# Initial condition 
equilibria = get_equilibria(f, μ0, domain=[-10,10])
x0 = [equilibria.stable[2], μ0]

# Timestep 
δt = 5e-3
