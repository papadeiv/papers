"""
    Simulation script

Storage of the definitions of the system alongside all the settings of the problem.
"""

# System parameters
r = 1.0                                       # Growth rate
k = 10.0                                      # Carrying capacity
h = 1.0                                       # Half-grazing biomass
c0 = 1.00                                     # Initial value of the bifurcation parameter for the continuation 
cf = 3.00                                     # Final value of the bifurcation parameter for the continuation 
c = collect(range(1.2, 2.8, 5))               # Slices of the parameter space
D_set = collect(range(1e-2, 1, 30))           # Set of diffusion values

# Dynamical system  
f(x, μ) = r*x*(1 - x/k) - μ*((x^2)/((x^2) + (h^2)))
