"""
    Simulation script

This is where we store the definition of the system alongside all the settings of the problem.
"""

# System parameters
μ = [0.0, 3.0, 1.0]                           # Bifurcation parameter value
σ = 0.125::Float64                            # Noise level (additive)
D = (σ^2)/2.0                                 # Diffusion level (additive) 

# Dynamical system  
f(x, μ) = μ[1] + 2.0*μ[2]*x + 3.0*μ[3]*(x^2)  # Drift
η(x) = σ                                      # Diffusion

# Initial condition 
equilibria = get_equilibria(f, μ, domain=[-10,10])
x0 = equilibria.stable[1]

# Time parameters
δt = 1e-3                                     # Timestep
Nt = convert(Int64, 1e5)                      # Total number of steps

# Ensemble parameters
Ne = convert(Int64, 1e3)
