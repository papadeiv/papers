"""
    Simulation script

This is where we store the definition of the system alongside all the settings of the problem.
"""

# System parameters
μ = 1.20::Float64                             # Bifurcation parameter value
σ = 0.125::Float64                            # Noise level (additive)
D = (σ^2)/2.0                                 # Diffusion level (additive) 

# Dynamical system  
f(x, μ) = -μ - 2*x + 3*(x^2) - (4/5)*(x^3)    # Drift
η(x) = σ                                      # Diffusion

# Equilibria
equilibria = get_equilibria(f, μ, domain=[-10,10])
stable = equilibria[1]
unstable = equilibria[2]

# Initial condition 
x0 = stable[2]

# Time parameters
δt = 1e-3                                     # Timestep
Nt = convert(Int64, 1e5)                      # Total number of steps

# Ensemble parameters
Ne = convert(Int64, 1e3)
