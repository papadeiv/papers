include("../../../../inc/IO.jl")
using LinearAlgebra, SparseArrays

# Domain and FEM discretization
const N = 200                      # Number of grid points
const x = LinRange(-1, 1, N)
const h = step(x)

# Import the data from csv
A = readin("../data/stiffness.csv")
X = readin("../data/snapshots.csv")

# Perform SVD on X
V, S, _ = svd(X)

# Truncate to the optimal number of singular vectors
energy = cumsum(S.^2) / sum(S.^2)
r = findfirst(e -> e > 0.999, energy)

# Build the basis of the ROM
Vr = V[:, 1:r]

# Export the data
writeout(S, "../data/singular_values.csv")
writeout(V, "../data/singular_vector.csv")
writeout(Vr, "../data/pod_basis.csv")

# Define parameterized source
function f(x, mu)
    return sin(pi * x) + mu * cos(2pi * x)
end

function assemble_b(f::Function)
    b = f.(x)
    b[1] = 0.0
    b[end] = 0.0
    return b
end

# Define parameter value at online stage 
μ_test = 0.73

# Compute the source term and the FOM solution
b = assemble_b(x -> f(x, μ_test))
t = time()
u_fom = A\b
Δt_FOM = time() - t

# Compute the ROM
b_r = transpose(Vr)*b
A_r = transpose(Vr)*A*Vr 

# Compute the reduced coefficients and thus the ROM solution 
t = time()
u_bar = A_r\b_r
Δt_ROM = time() - t
u_rom = Vr*u_bar

display(1 - (Δt_ROM/Δt_FOM))

# Export the data
writeout(u_bar, "../data/coefficients.csv")
writeout(hcat(u_fom, u_rom), "../data/solutions.csv")
