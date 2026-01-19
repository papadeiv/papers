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
C = [1.9, -0.6, -7.0, -2.0, +1.7, +2.4, +3.2, +0.4, +2.2]
g(t) = (C[1] + 2*C[2]*(tanh(t)) + 3*C[3]*(tanh(t))^2 + 4*C[4]*(tanh(t))^3 + 5*C[5]*(tanh(t))^4 + 6*C[6]*(tanh(t))^5 + 7*C[7]*(tanh(t))^6 + 8*C[8]*(tanh(t))^7 + 9*C[9]*(tanh(t))^8)*(sech(t))^2 

# Define the dynamical system (diffusion of the state variable) 
σ = 0.00::Float64
η(x) = σ

# Define the the past limit system
T = 10.0
t∞ = [-T, +T]
x∞ = 0.01742
λ∞ = -2.0
u∞ = [x∞, λ∞]

# Solve the fast-slow SDE
t, λ, x = evolve_shifted_1d(f, g, η, u∞, t∞, Nt=convert(Int64,1e4))

# Export the data 
writeout(hcat(t, λ, x), "../data/figure_03/solution.csv")

# Execute the postprocessing and plotting scripts
include("../postprocessing/figure_03_postprocessing.jl")
include("../plotting/figure_03_plotting.jl")
