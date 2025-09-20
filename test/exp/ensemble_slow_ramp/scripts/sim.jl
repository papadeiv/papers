"""
    Simulation script

This is where we store the definition of the system alongside all the settings of the problem.
"""

# System parameters
μ0 = collect(1.2:0.1:1.3)                       # List of initial parameter values
δμ = 5e-3                                       # Range of the parameter's ramp
μf = round.(μ0 .+ δμ, digits=3)                 # List of final parameter value
ε = 1e-3                                        # Slow timescale
σ = 0.250::Float64                              # Noise level (additive)
D = (σ^2)/2.0                                   # Diffusion level (additive) 

# Dynamical system  
f(x, μ) = -μ - 2*x + 3*(x^2) - (4/5)*(x^3)      # Drift
g(t) = ε                                        # Shift/Ramp
η(x) = σ                                        # Diffusion

# Initial conditions
equilibria = [(get_equilibria(f, μ, domain=[-10,10])).stable[2] for μ in μ0]
x0 = hcat(equilibria, μ0)

# Timestep
δt = 1e-3

# Ensemble parameters
Ne = convert(Int64, 1e3)
