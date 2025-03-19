include("../../../../inc/IO.jl")
include("../../../../inc/SystemAnalysis.jl")

# Define a 2-dimensional bistable saddle-node normal form
f1(x, y, μ) = -μ - 2*x + 3*(x^2) - (4/5)*(x^3) + y
f2(x, y, μ) = -μ - y
g(x, y) = 0.0::Float64 

# Define the parameter range
μ_values = [0.0, 0.8108]

# Define the IC
x0 = 1.25::Float64
y0 = 0.75::Float64
u0 = [x0, y0]

# Solve the fast-slow SDE
t, μ, u = evolve_fixed_2d(f1, f2, g, g, μ_values, u0, δt=1e-3, Nt=convert(Int64,1e4), Nμ=30::Int64)
x = u[1]
y = u[2]

# Export the data 
writeout(t, "../data/slide_02/time.csv")
writeout(μ, "../data/slide_02/μ.csv")
writeout(x, "../data/slide_02/x1.csv")
writeout(y, "../data/slide_02/x2.csv")

# Execute the plotting scripts
include("../plotting/slide_02_plotting.jl")
