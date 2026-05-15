include("../sim/solve_nonlinear_burgers_1d.jl")
include("../../../../inc/IO.jl")
using LinearAlgebra, SparseArrays

# ==========================================
# | Proper Orthogonal Decomposition (POD)
# ==========================================

Ω = [-1,1]    # Spatial domain 
BC = [0,0]    # Boundary conditions (BCs)
Nh = convert(Int64,1e3)   # Number of DOFs 
T = 1.0    # Time horizon 
δt = 1e-3   # Timestep
Nt = Int(T/δt)    # Number of timesteps
μ_test = 0.0002   # Test value of the parameter

nodes = LinRange(-1.0, 1.0, Nh+1)
internal_indices = 2:(Nh)

# Import the data from csv
X = readin("../data/burgers_eq/snapshots.csv")

# ==========================================
# | Offline phase
# ==========================================

println("--- Extraction of POD Basis ---")

# Perform SVD on X
V, S, _ = svd(X)

# Truncate to the optimal number of singular vectors
energy = cumsum(S.^2)/sum(S.^2)
r = findfirst(e -> e > 0.999, energy)
println("Extracted $r modes to retain 99% of the energy.")

# Build the basis of the ROM
Vr = V[:, 1:r]

# Export the data
writeout(S, "../data/burgers_eq/singular_values.csv")
writeout(V, "../data/burgers_eq/singular_vector.csv")
writeout(Vr, "../data/burgers_eq/pod_basis.csv")

# ==========================================
# | Online phase
# ==========================================

println("\n--- Online Stage: Reduced-Order Model ---")

# Compute the FOM solution
output = solve_nonlinear_burgers_1d(μ_test, Ω, Nh, BC, T, δt)
u_fom = output[1]

# Extract the full-order matrices
M_int = output[2]
K_int_int = output[3]
U0 = output[4]
h = output[5]

# Compute the reduced-order matrices
Mr = transpose(Vr)*M_int*Vr
Kr = transpose(Vr)*K_int_int*Vr 

# Compute the reduced-order initial condition
U_reduced_initial = transpose(Vr)*U0

# Storage for the reconstructed full-order solution
u_bar = Matrix{Float64}(undef, r, Nt+1)
u_bar[:,1] = U_reduced_initial
u_rom = Matrix{Float64}(undef, Nh+1, Nt+1)
u_rom[:,1] = [BC[1]; U0; BC[2]]   # Store the full IC

# Time loop for the reduced-order model
U_current_reduced = copy(U_reduced_initial)
for n in 1:Nt
        # Call the newton_raphson function with the reduced mass matrix (Mr)
        # and explicitly pass the Vr matrix for reconstruction and projection.
        U_next_reduced = newton_raphson(U_current_reduced, U_current_reduced, Mr, Kr, μ_test, δt, h, Nh+1, internal_indices; Phi=Vr)
    
        # Reconstruct the full solution from the reduced solution for plotting
        U_reconstructed_int = Vr*U_next_reduced
        U_reconstructed_full = [BC[1]; U_reconstructed_int; BC[2]] # Add homogeneous BCs back

        # Update the solutions
        u_bar[:, n+1] = U_next_reduced
        u_rom[:, n+1] = U_reconstructed_full

        # Update the iterative previous solution
        global U_current_reduced = U_next_reduced
end

# Export the data
writeout(u_bar, "../data/burgers_eq/coefficients.csv")
writeout(u_fom, "../data/burgers_eq/FOM_solutions.csv")
writeout(u_rom, "../data/burgers_eq/ROM_solutions.csv")
