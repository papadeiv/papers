include("../../../../inc/IO.jl")
include("../../../../inc/LatticeDynamics.jl")

# Construct the lattice
n = 100
m = 100 
N = n*m
L = Lattice(n, m)

# Model's parameters
k = 14.1::Float64 # (airway) smooth muscle mass
A = 0.63::Float64 # inter-airway coupling
Pi = 0.96::Float64 # inflection point of pressure-radius interdependence
Pb = 0.0::Float64 # breathing pressure
Pb0 = 7.25::Float64 # breating pressure's IC

# IC on the bifurcation parameter (smooth muscle activation)
K0 = 0.5::Float64

# IC of the lattice solution (augmented by the scalar, ramped parameter K) 
u0 = 0.93::Float64.*ones(Float64, length(L.grid))
push!(u0, K0)

# Define the noise level and timescale separation constants... 
σ = 0.080::Float64
ε = 0.010::Float64

# ... and associated functions
g(u, j, lattice) = σ
h(t) = ε

# Definition of the lattice deterministic dynamics with first-order nearest-neighbours coupling (discretisation of the laplacian)
f(u, j, lattice) = 1.0::Float64/(1.0::Float64 + exp(-(Pb0*N)/sum(u[1:(end-1)].^4) + u[end]*(k/u[j]) - (Pb0*N)/sum(u[1:(end-1)].^4)*A*(u[j] + sum([u[k] for k in get_neighbours(lattice,j)].^4))*(1.0::Float64 - u[j] + 1.5::Float64*(1-u[j])^2) + Pi)) - u[j]

# Solve the lattice SDE
t, K, u = evolve_lattice(L, f, g, h, u0, δt=1e-1, μf=1.5::Float64)

# Assemble data matrix of the lattice solution for export
u_mat = Matrix{Float64}(undef, N, length(u))
for n in 1:length(u)
        u_mat[:,n] = u[n] 
end

# Export the data
writeout(t, "../data/slide_08/time.csv")
writeout(K, "../data/slide_08/K.csv")
writeout(u_mat, "../data/slide_08/u.csv")

# Execute the postprocessing and plotting scripts
include("../postprocessing/slide_08_postprocessing.jl")
include("../plotting/slide_08_plotting.jl")

#=
# Definition of the independent vector variable u
u = zeros(Float64, length(L.grid))
u[6] = 1
u[10] = 1

# Definition of the ramped scalar bifurcation parameter scorporated from u
μ = 0.5::Float64

# Most abstract definition of f (taken as composition with sigmoid and r functions)
r(u, μ, j, lattice) = -(Pb0*N)/sum(u.^4) + μ*(k/u[j]) - (Pb0*N)/sum(u.^4)*A*(u[j] + sum([u[k] for k in get_neighbours(lattice,j)].^4))*(1.0::Float64 - u[j] + 1.5::Float64*(1-u[j])^2) + Pi
sigmoid(x) = 1.0::Float64/(1.0::Float64 + exp(-x))
f(u, μ, j, lattice) = sigmoid(-r(u, μ, j, lattice)) - u[j] 

# Definition of f without the definition of sigmoid
f(u, μ, j, lattice) = 1.0::Float64/(1.0::Float64 + exp(r(u, μ, j, lattice))) - u[j]

# Definition of f without the definition of sigmoid and r but with the parameter μ scorporated by the vector variable u
f(u, μ, j, lattice) = 1.0::Float64/(1.0::Float64 + exp(-(Pb0*N)/sum(u.^4) + μ*(k/u[j]) - (Pb0*N)/sum(u.^4)*A*(u[j] + sum([u[k] for k in get_neighbours(lattice,j)].^4))*(1.0::Float64 - u[j] + 1.5::Float64*(1-u[j])^2) + Pi)) - u[j]
=#
