include("./solve_nonlinear_burgers_1d.jl")
include("../../../../inc/IO.jl")
using LinearAlgebra

# ==========================================
# | Discretisation settings
# ==========================================

Ω = [-1,1]    # Spatial domain
BC = [0,0]    # Boundary conditions (BCs)
Nh = convert(Int64,1e3)   # Number of DOFs
T = 1.0    # Time horizon
δt = 1e-3   # Timestep
Nt = Int(T/δt)    # Number of timesteps
μ = 0.001   # Training set

# ==========================================
# | Offline phase
# ==========================================
 
println("--- Offline Stage: Full-Order Model for Snapshot Generation ---")

X = Matrix{Float64}(undef, Nh-1, Nt+1)

# Assemble and solve the FOM
output = solve_nonlinear_burgers_1d(μ, Ω, Nh, BC, T, δt)
X = output[1]

# ==========================================
# | Wrap-up
# ==========================================

# Export the data
writeout(X, "../data/burgers_eq/snapshots.csv")

# Execute the postprocessing and plotting scripts
include("../postprocessing/fig:burgers_pod.jl")
include("../plotting/fig:burgers_pod.jl")
