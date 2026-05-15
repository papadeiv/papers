include("./solve_heat_equation_1d.jl")
include("../../../../inc/IO.jl")
using LinearAlgebra

# ==========================================
# | Discretisation settings
# ==========================================

Ω = [-1,1]    # Spatial domain
BC = [1,2]    # Boundary conditions (BCs)
Nh = convert(Int64,1e3)   # Number of DOFs
T = 1.0    # Time horizon
δt = 1e-3   # Timestep
Nt = Int(T/δt)    # Number of timesteps

# ==========================================
# | Training settings
# ==========================================

P1 = collect(LinRange(0.025, 0.1, 1000))   # Training set of problem 1
P2 = collect(LinRange(0.1, 0.25, 10))    # Training set of problem 2

# ==========================================
# | Offline phase: problem 1
# ==========================================
 
Nμ = length(P1)   # Number of training values of μ
Nf = Nμ   # Number of snapshots collected

X1 = Matrix{Float64}(undef, Nh-1, Nf)   # Snapshot matrix for problem 1

# Loop over values in the training set
@showprogress for n in 1:Nμ
        μ = P1[n]    # Extract the parameter value from the training set

        # Assemble and solve the FOM
        output = solve_heat_equation_1d(μ, Ω, Nh, BC, T, δt)
        Uh = output[1]

        # Fill-in the snapshot matrix
        X1[:,n] = Uh[end]
end

# ==========================================
# | Offline phase: problem 2
# ==========================================

Nμ = length(P2)   # Number of training values of μ

X2 = Matrix{Float64}(undef, Nh-1, Nf+1)   # Snapshot matrix for problem 1

# Define a set of indices to subsample the desired timesteps (1 every 10)
T_sub = collect(StepRange(1, Int64(10), 1001)) 

# Loop over values in the training set
@showprogress for n in 1:Nμ
        μ = P2[n]    # Extract the parameter value from the training set

        # Assemble and solve the FOM
        output = solve_heat_equation_1d(μ, Ω, Nh, BC, T, δt)
        Uh = output[1]

        # Fill-in the snapshot matrix 
        τ = (n-1)*100 + 1 
        for t in T_sub
                X2[:,τ] = Uh[t]
                τ += 1
        end
end

# ==========================================
# | Wrap-up
# ==========================================

# Export the data
writeout(X1, "../data/heat_eq/problem_1/snapshots.csv")
writeout(X2, "../data/heat_eq/problem_2/snapshots.csv")

# Execute the postprocessing and plotting scripts
include("../postprocessing/fig:heat_pod.jl")
include("../plotting/fig:heat_pod.jl")
