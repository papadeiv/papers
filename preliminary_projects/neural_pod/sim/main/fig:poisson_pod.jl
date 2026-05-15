include("../../../../inc/IO.jl")
using LinearAlgebra, SparseArrays

# Domain and FEM discretization
const N = 200                      # Number of grid points
const x = LinRange(-1, 1, N)
const h = step(x)

# FEM assembly: stiffness matrix and load vector
function assemble_A(N, h)
    A = spzeros(N, N)
    for i in 2:N-1
        A[i, i-1] = -1.0
        A[i, i]   = 2.0
        A[i, i+1] = -1.0
    end
    A = A / h^2
    A[1, :] .= 0.0; A[end, :] .= 0.0
    A[:, 1] .= 0.0; A[:, end] .= 0.0
    A[1,1] = 1.0; A[end,end] = 1.0
    return A
end

function assemble_b(f::Function)
    b = f.(x)
    b[1] = 0.0
    b[end] = 0.0
    return b
end

# Parameterized source function
function f(x, mu)
    return sin(pi * x) + mu * cos(2pi * x)
end

# Offline stage
μ_train = range(0.1, 1.0, length=20)
snapshots = Float64[]
A = assemble_A(N, h)

X = zeros(N, length(μ_train))
for (j, μ) in enumerate(μ_train)
    b = assemble_b(x -> f(x, μ))
    u = A\b
    X[:, j] = u
end

# Export data
writeout(A, "../data/stiffness.csv")
writeout(X, "../data/snapshots.csv")

# Execute postprocessing and plotting scripts
include("../postprocessing/fig:poisson_pod.jl")
include("../plotting/fig:poisson_pod.jl")
