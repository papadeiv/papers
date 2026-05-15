include("../sim/solve_heat_equation_1d.jl")
include("../../../../inc/IO.jl")
using LinearAlgebra, SparseArrays

# ==========================================
# | Proper Orthogonal Decomposition (POD)
# ==========================================

Ω = [-1,1]    # Spatial domain 
BC = [1,2]    # Boundary conditions (BCs)
Nh = convert(Int64,1e3)   # Number of DOFs 
T = 1.0    # Time horizon 
δt = 1e-3   # Timestep
Nt = Int(T/δt)    # Number of timesteps
μ_test = 0.15   # Test value of the parameter

# Import the data from csv
X1 = readin("../data/heat_eq/problem_1/snapshots.csv")
X2 = readin("../data/heat_eq/problem_2/snapshots.csv")

# ==========================================
# | Offline phase: problem 1
# ==========================================

# Perform SVD on X
V, S, _ = svd(X1)

# Truncate to the optimal number of singular vectors
energy = cumsum(S.^2)/sum(S.^2)
r = findfirst(e -> e > 0.999, energy)

# Build the basis of the ROM
Vr = V[:, 1:r]

# Export the data
writeout(S, "../data/heat_eq/problem_1/singular_values.csv")
writeout(V, "../data/heat_eq/problem_1/singular_vector.csv")
writeout(Vr, "../data/heat_eq/problem_1/pod_basis.csv")

# ==========================================
# | Online phase: problem 1
# ==========================================

# Compute the FOM solution
t = time()
output = solve_heat_equation_1d(μ_test, Ω, Nh, BC, T, δt)
Δt_FOM = time() - t

# Extract the steady-state
Uh = output[1]
u_fom = Uh[end]

# Extract the FOM model 
K_int_int = output[3]
K_int_bc = output[4]

# Compute the ROM solution
t = time()
Ar = transpose(Vr)*K_int_int*Vr
br = -transpose(Vr)*K_int_bc*BC
u_bar = Ar\br
u_rom = Vr*u_bar
Δt_ROM = time() - t

# Export the data
display(1 - (Δt_ROM/Δt_FOM))
writeout(u_bar, "../data/heat_eq/problem_1/coefficients.csv")
writeout(hcat([BC[1]; u_fom; BC[2]], [BC[1]; u_rom; BC[2]]), "../data/heat_eq/problem_1/solutions.csv")

# ==========================================
# | Offline phase: problem 2
# ==========================================

# Perform SVD on X
V, S, _ = svd(X2)

# Truncate to the optimal number of singular vectors
energy = cumsum(S.^2) / sum(S.^2)
r = findfirst(e -> e > 0.999999, energy)

# Build the basis of the ROM
Vr = V[:, 1:r]

# Export the data
writeout(S, "../data/heat_eq/problem_2/singular_values.csv")
writeout(V, "../data/heat_eq/problem_2/singular_vector.csv")
writeout(Vr, "../data/heat_eq/problem_2/pod_basis.csv")

# ==========================================
# | Online phase: problem 2
# ==========================================

# Compute the FOM solution
output = solve_heat_equation_1d(μ_test, Ω, Nh, BC, T, δt)
Uh = output[1]

# Format the FOM solution for export
u_fom = Matrix{Float64}(undef, Nh+1, Nt+1)
for n in 1:length(Uh)
        u_fom[:,n] = [BC[1]; Uh[n]; BC[2]]
end

# Compute the reduced-order matrices and forcing vector
M_int = output[2]
K_int_int = output[3]
K_int_bc = output[4]
U0 = output[5]
Mr = transpose(Vr)*M_int*Vr 
Kr = transpose(Vr)*K_int_int*Vr
Fr = -transpose(Vr)*K_int_bc*BC

# Compute the reduced-order initial condition
U_reduced_initial = transpose(Vr)*U0

# The LHS matrix for the reduced-order system (Implicit Euler)
A_reduced = Mr + δt*Kr

# Storage for the reconstructed full-order solution at each online time step
u_bar = Matrix{Float64}(undef, r, Nt+1)
u_bar[:,1] = U_reduced_initial
u_rom = Matrix{Float64}(undef, Nh+1, Nt+1)
u_rom[:,1] = [BC[1]; U0; BC[2]] # Store the full IC

# Time loop for the reduced-order model
U_current_reduced = copy(U_reduced_initial)
@showprogress for n in 1:Nt
        # RHS for the reduced-order system
        b_reduced = Mr*U_current_reduced + δt*Fr
    
        # Solve the small, reduced-order linear system
        U_next_reduced = A_reduced\b_reduced
    
        # Reconstruct the full solution for plotting
        U_reconstructed_int = Vr*U_next_reduced
        U_reconstructed_full = [BC[1]; U_reconstructed_int; BC[2]]

        # Update the solutions
        u_bar[:, n+1] = U_next_reduced
        u_rom[:, n+1] = U_reconstructed_full

        # Update the iterative previous solution
        global U_current_reduced = U_next_reduced
end

# Export the data
writeout(u_bar, "../data/heat_eq/problem_2/coefficients.csv")
writeout(u_fom, "../data/heat_eq/problem_2/FOM_solutions.csv")
writeout(u_rom, "../data/heat_eq/problem_2/ROM_solutions.csv")
