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

# Define the coefficients of the parameter's shift 
C = [-1.6, -2.0, 2.6, 2.0]

# Define the dynamical system (diffusion of the state variable) 
σ = 0.00::Float64
η(x) = σ

# Define the the past limit system
T = 15.0
t∞ = [-T, +T]
x∞ = +0.26630
λ∞ = -0.99512
u∞ = [x∞, λ∞]

# Define the dynamical system (drift of the bifurcation parameter)
G(t, n) = (n*(t)^(n-1))/(1 + t^2)^(n/2 +1)
g(t) = C[1]*G(t,1) + C[2]*G(t,2) + C[3]*G(t,3) + C[4]*G(t,4) 

# Solve the fast-slow SDE
t, λ, x = evolve_shifted_1d(f, g, η, u∞, t∞, Nt=convert(Int64,1e4))

# Export the data 
writeout(hcat(t, λ, x), "../data/figure_05/solution.csv")

# Execute the postprocessing and plotting scripts
include("../postprocessing/figure_05_postprocessing.jl")
include("../plotting/figure_05_plotting.jl")
