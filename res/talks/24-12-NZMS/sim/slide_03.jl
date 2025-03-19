include("../../../../inc/IO.jl")
include("../../../../inc/SystemAnalysis.jl")

# Define the stationary parameters of the process
σ = 0.100::Float64
ε = 0.010::Float64

# Define a 2-dimensional bistable saddle-node normal form
f1(x, y, μ) = -μ - 2*x + 3*(x^2) - (4/5)*(x^3) + y
f2(x, y, μ) = -μ - y
g(x, y) = σ
h(t) = ε

# Define the IC
x0 = 1.25::Float64
y0 = 0.75::Float64
μ0 = 0.25::Float64
u0 = [x0, y0, μ0]

# Solve the fast-slow SDE
t, μ, u = evolve_ramped_2d(f1, f2, g, g, h, u0, δt=1e-3, saveat=0.10, μf=1.25)
x = u[:,1]
y = u[:,2]

# Export the data 
writeout(t, "../data/slide_03/time.csv")
writeout(μ, "../data/slide_03/μ.csv")
writeout(x, "../data/slide_03/x1.csv")
writeout(y, "../data/slide_03/x2.csv")

# Execute the plotting scripts
include("../plotting/slide_03_plotting.jl")
