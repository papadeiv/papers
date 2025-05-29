include("../../../../inc/IO.jl")
include("../../../../inc/SystemAnalysis.jl")

# Parameters of the payoff matrix
a = 1.00::Float64
b = 1.00::Float64
c = 1.00::Float64
d = 1.00::Float64

# Parameters of the control matrix
G11 = 1.00::Float64
G12 = 1.00::Float64
G21 = 1.00::Float64
G22 = 1.00::Float64

# Define the timescale constant for the parameter ramp
ε = 0.010::Float64
g(t) = ε

# Define the dynamical system (deterministic drift) 
f1(x, y) = x*(1.00::Float64 - x)*((a + d - b - c)*x + b - d + (G11 - G21)*y*x + (G12 - G22)*y*(1.00::Float64 - x))
f2(x, y) = g(x)*y

# Define the dynamical system (stochastic diffusion) 
g1(x, μ) = 0.0::Float64
g2(x, μ) = 0.0::Float64

#=
# Define the IC
x1 = 3.33::Float64
x2 = -(3.0::Float64)
x0 = [x1, x2]

# Solve the fast-slow SDE
t, μ, u = evolve_fixed_2d(f1, f2, x0, Nt=convert(Int64,1e3))

# Export the data 
writeout(t, "../data/figure_01/time.csv")
writeout(μ, "../data/figure_01/μ.csv")
writeout(u, "../data/figure_01/u.csv")

# Execute the postprocessing and plotting scripts
include("../postprocessing/figure_01_postprocessing.jl")
include("../plotting/figure_01_plotting.jl")
=#
