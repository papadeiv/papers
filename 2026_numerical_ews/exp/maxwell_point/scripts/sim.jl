"""
    Simulation script

Storage of the definitions of the system alongside all the settings of the problem.
"""

# System parameters
μ = 0.0                                       # Bifurcation parameter value
λ = -1.0                                      # Bifurcation parameter value
ε = 0.0                                       # Timescale separation
σ = 0.25                                      # Noise level (additive)
D = (σ^2)/2.0                                 # Diffusion level (additive) 

# Polynomial coefficients of the true potential
c0 = 0.3
c1 = μ
c2 = λ 
c3 = 0.0
c4 = 2.0

# Dynamical system  
f(x, μ) = -(c1 + 2*c2*x + 3*c3*x^2 + 4*c4*x^3)# Drift
Λ(t) = ε                                      # Shift/Ramp
η(x) = σ                                      # Diffusion

# Scalar potential of the conservative system 
U(x, μ) = c0 + c1*x + c2*x^2 + c3*x^3 + c4*x^4              # Ground truth
V(x, c) = c[1]*x + c[2]*(x^2) + c[3]*(x^3) + c[4]*(x^4)     # Arbitrary

# Stationary equilibrium distribution
P(x, μ) = exp(-(U(x,μ)/D))
Z = normalise(P, μ, (-10,10))
display(Z)
ρ(x, μ) = Z*P(x, μ)

# Target (arbitrary) distribution 
g(x, c) = exp(-(1.0::Float64/D)*(V(x, c)))
N(c) = normalise(g, c, (-10,10))
p(x, c) = N(c)*exp.(-(1.0::Float64/D).*(c[1]*x .+ c[2]*(x.^2) .+ c[3]*(x.^3) .+ c[4]*(x.^4)))

# Initial condition 
equilibria = get_equilibria(f, μ, domain=[-10,10])
x0 = [equilibria.stable[1], μ]

# Simulation parameters
dt = 5e-2                                     # Timestep
Nt = 5e6                                      # Total number of steps
Ne = 1e1                                      # Number of particles in the ensemble 
