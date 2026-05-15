include("../../../../inc/IO.jl")
include("../../../../inc/SystemAnalysis.jl")

printstyled("Solving the non-autonomous SDEs\n"; bold=true, underline=true, color=:light_blue)

# Define the dynamical system (drift of the state variable) 
a = -(0.25::Float64)
b = 1.20::Float64
c = -(0.40::Float64)
d = -(0.30::Float64)
e = 3.00::Float64
K = 2.00::Float64
f(x, λ) = -((x + a + b*λ)^2 + c*tanh(λ - d))*(x - K/(cosh(e*λ)))

# Define the dynamical system (drift of the bifurcation parameter)
C = [-1.6, -2.0, 2.6, 2.0]
g(t) = (C[1] + 2*C[2]*(tanh(t)) + 3*C[3]*(tanh(t))^2 + 4*C[4]*(tanh(t))^3)*(sech(t))^2 

# Define the dynamical system (diffusion of the state variable) 
σ = 0.00::Float64
η(x) = σ

# Define the the past limit system
T = 10.0
t∞ = [-T, +T]
x∞ = 0.26
λ∞ = -1.0
u∞ = [x∞, λ∞]

# Solve the fast-slow SDE
t, λ, x = evolve_shifted_1d(f, g, η, u∞, t∞, Nt=convert(Int64,1e4))

# Export the data 
writeout(hcat(t, λ, x), "../data/figure_02/solution.csv")

# Execute the postprocessing and plotting scripts
include("../postprocessing/figure_02.b_postprocessing.jl")
include("../plotting/figure_02.b_plotting.jl")
