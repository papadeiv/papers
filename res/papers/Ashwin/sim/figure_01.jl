include("../../../../inc/IO.jl")
include("../../../../inc/SystemAnalysis.jl")

# Define the dynamical system (deterministic drift) 
a = -(0.25::Float64)
b = 1.20::Float64
c = -(0.40::Float64)
d = -(0.30::Float64)
e = 3.00::Float64
K = 2.00::Float64
f(x, λ) = -((x + a + b*λ)^2 + c*tanh(λ - d))*(x - K/(cosh(e*λ)))

# Define the parameter's rates
ε = collect(LinRange(0.01, 0.15, 1000)) 
Nε = length(ε)

# Define the dynamical system (stochastic diffusion) 
σ = 0.00::Float64
η(x) = σ

# Define the state in the past limit system
x∞ = 0.2
λ∞ = -1.0
u∞ = [x∞, λ∞]

# Define the compactified time-domain
T = 200.0
t∞ = [-T, +T]

# Loop over the rates
printstyled("Solving the non-autonomous SDEs for different rates\n"; bold=true, underline=true, color=:light_blue)
@showprogress for n in 1:Nε 
        # Define the dynamical system (parameter shift)
        Λ(t) = ε[n]*(sech(ε[n]*t))^2

        # Solve the fast-slow SDE
        local t, λ, x = evolve_shifted_1d(f, Λ, η, u∞, t∞, Nt=convert(Int64,1e4))

        # Export the data 
        writeout(hcat(t, λ, x), "../data/figure_01/solutions/$n.csv")
end

# Export the parameter rates
writeout(ε, "../data/figure_01/rates.csv")

# Execute the postprocessing and plotting scripts
include("../postprocessing/figure_01_postprocessing.jl")
include("../plotting/figure_01_plotting.jl")
