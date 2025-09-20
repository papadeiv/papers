"""
    Simulation script

This is where we store the definition of the system alongside all the settings of the problem.
"""

# System parameters
μ = 1.00::Float64                             # Centre 
θ = 1.00::Float64                             # Curvature 
σ = 0.125::Float64                            # Noise level (additive)
D = (σ^2)/2.0                                 # Diffusion level (additive) 

# Dynamical system  
f(x, λ) = -λ*(μ - x)                          # Drift
η(x) = σ                                      # Diffusion

# Initial condition 
equilibria = get_equilibria(f, θ, domain=[-10,10])
x0 = equilibria.stable[2]

# Time parameters
δt = 1e-3                                     # Timestep
Nt = convert(Int64, 1e5)                      # Total number of steps

# Ensemble parameters
Ne = convert(Int64, 1e3)
