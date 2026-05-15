"""
    Simulation script

Storage of the definitions of the system alongside all the settings of the problem.
"""

# System parameters
μ0 = 0.00::Float64                            # Initial parameter value
μf = 1.50::Float64                            # Final parameter value
ε = 1e-3                                      # Slow timescale
σ = 0.100::Float64                            # Noise level (additive)
D = (σ^2)/2.0                                 # Diffusion level (additive) 

# Dynamical system  
f(x, μ) = -μ + 4*x - 4*x^3                    # Drift
Λ(t) = ε                                      # Shift/Ramp
η(x) = σ                                      # Diffusion

# Initial condition 
equilibria = get_equilibria(f, μ0, domain=[-10,10])
x0 = [equilibria.stable[2], μ0]

# Timestep 
δt = 1e-2

# Counter for the parallelized for-loop
counter = Base.Threads.Atomic{Int}(0)
print_lock = ReentrantLock()
