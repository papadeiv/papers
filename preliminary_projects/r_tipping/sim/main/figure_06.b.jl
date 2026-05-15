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
G(t, n) = n*(((sinh(t))^(n-1))/((cosh(t))^(n+1))) 
C = [0.0, 0.0, 0.0, 0.0, -1.6, 0.0, 2.0, 0.0, 1.1, 0.0, -0.5]
g(t) = C[1]*G(t,1) + C[2]*G(t,2) + C[3]*G(t,3) + C[4]*G(t,4) + C[5]*G(t,5) + C[6]*G(t,6) + C[7]*G(t,7) + C[8]*G(t,8) + C[9]*G(t,9) + C[10]*G(t,10) + C[11]*G(t,11)

# Define the dynamical system (diffusion of the state variable) 
σ = 0.00::Float64
η(x) = σ

# Define the the past limit system
T = 10.0
t∞ = [-T, +T]
x∞ = +0.26281
λ∞ = -1.0
u∞ = [x∞, λ∞]

# Solve the fast-slow SDE
t, λ, x = evolve_shifted_1d(f, g, η, u∞, t∞, Nt=convert(Int64,1e4))

# Export the data 
writeout(hcat(t, λ, x), "../data/figure_06/solution.csv")

# Execute the postprocessing and plotting scripts
include("../postprocessing/figure_06.b_postprocessing.jl")
include("../plotting/figure_06.b_plotting.jl")
