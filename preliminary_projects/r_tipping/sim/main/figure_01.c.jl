include("../../../../inc/IO.jl")
include("../../../../inc/SystemAnalysis.jl")

printstyled("Solving the non-autonomous ensemble problem\n"; bold=true, underline=true, color=:light_blue)

# Define the dynamical system (drift of the state variable) 
a = -(0.25::Float64)
b = 1.20::Float64
c = -(0.40::Float64)
d = -(0.30::Float64)
e = 3.00::Float64
K = 2.00::Float64
f(x, λ) = -((x + a + b*λ)^2 + c*tanh(λ - d))*(x - K/(cosh(e*λ)))

# Define the parameter's rates
ε = [0.10,0.11,0.12,0.13,0.14] 
Nε = length(ε)

# Define the dynamical system (diffusion of the state variable) 
σ = 0.00::Float64
η(x) = σ

# Define the the past limit system
T = 200.0
t∞ = [-T, +T]
x∞ = 0.2
λ∞ = -1.0
u∞ = [x∞, λ∞]

# Define the parameter shift
Λ(t) = tanh(ε[n]*t)

# Compute the (time) integration constant for the specified past limit system
C = λ∞ - Λ(t∞[1])
        
# Define the dynamical system (drift of the bifurcation parameter)
g(t) = ε[n]*(sech(ε[n]*t))^2

# Solve the fast-slow SDE
local t, λ, x = evolve_shifted_1d(f, g, η, u∞, t∞, Nt=convert(Int64,1e4))

# Export the data 
writeout(hcat(t, λ, x), "../data/figure_01/solutions/$n.csv")

# Export the parameter rates
writeout(ε, "../data/figure_01/rates.csv")

# Execute the postprocessing and plotting scripts
include("../postprocessing/figure_01.c_postprocessing.jl")
include("../plotting/figure_01.c_plotting.jl")
