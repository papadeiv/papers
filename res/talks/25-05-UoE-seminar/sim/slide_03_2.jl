include("../../../../inc/IO.jl")
include("../../../../inc/SystemAnalysis.jl")

# Define the normal form
f(x, μ) = -μ - 2*x + 3*(x^2) - (4/5)*(x^3)
# Define the parameter range for the sweep
Nμ = convert(Int64,5e1)
μ0 = LinRange(1.000::Float64, 1.800::Float64, Nμ)
x0 = (1.800::Float64).*ones(Float64, Nμ)

# Solve the parameter sweep of the deterministic case
η(x) = 0.0::Float64 
t1, μ1, u1 = evolve_fixed_1d(f, η, [μ0[1],μ0[end]], x0, Nt=150, Nμ=Nμ)

# Define the timescale constant for the parameter ramp
ε = 0.001::Float64
g(t) = ε
# Define the IC for the non-autonomous case
u0 = [sort(get_equilibria(f, μ0[25], domain=[-20,20])[1], rev=true)[1], μ0[25]]

# Solve the non-autonomous deterministic case
t2, μ2, u2 = evolve_ramped_1d(f, g, η, u0, Nt=convert(Int64,1e3), μf=μ0[end])

# Define the noise level in the system 
σ = 0.100::Float64
η(x) = σ

# Solve the fast-slow SDE
t3, μ3, u3 = evolve_ramped_1d(f, g, η, u0, Nt=convert(Int64,1e3), μf=μ0[end])

# Export the data 
writeout(t1, "../data/slide_03/t1.csv")
writeout(u1, "../data/slide_03/u1.csv")
writeout(hcat(t2, μ2, u2, u3), "../data/slide_03/u2&3.csv")

# Execute the plotting script
include("../plotting/slide_03_2_plotting.jl")
