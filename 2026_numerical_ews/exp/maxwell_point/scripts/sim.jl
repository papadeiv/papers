"""
    Simulation script

Storage of the definitions of the system alongside all the settings of the problem.
"""

# System parameters
μ = 0.0                                       # Bifurcation parameter value
ε = 0.0                                       # Timescale separation
σ = 1.5                                       # Noise level (additive)
D = (σ^2)/2.0                                 # Diffusion level (additive) 

# Dynamical system  
f(x, μ) = -μ + 4*x - 4*x^3                    # Drift
Λ(t) = ε                                      # Shift/Ramp
η(x) = σ                                      # Diffusion

# Scalar potential of the conservative system 
U(x, μ) = 1 + μ*x - 2*x^2 + x^4                             # Ground truth
V(x, c) = c[1]*x + c[2]*(x^2) + c[3]*(x^3) + c[4]*(x^4)     # Arbitrary

# Stationary equilibrium distribution
Z = 2.07039469152
ρ(x, μ) = (1/Z)*exp(-(U(x,μ)/D))

# Target (arbitrary) distribution 
g(x, c) = exp(-(1.0::Float64/D)*(V(x, c)))
N(c) = normalise(g, c, (-10,10))
p(x, c) = N(c)*exp.(-(1.0::Float64/D).*(c[1]*x .+ c[2]*(x.^2) .+ c[3]*(x.^3) .+ c[4]*(x.^4)))

# Initial condition 
equilibria = get_equilibria(f, μ, domain=[-10,10])
x0 = [equilibria.stable[1], μ]

# Simulation parameters
dt = 5e-2                                     # Timestep
Nt = 2e5                                      # Total number of steps
Ne = 2e2                                      # Number of particles in the ensemble 
β = 1e-3                                      # Std of the guess perturbation 
Na = convert(Int64, 1e4)                      # Number of attempts per guess 
