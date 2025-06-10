include("../../../../inc/IO.jl")
include("../../../../inc/SystemAnalysis.jl")

# Define the noise level and timescale separation 
σ = 0.250::Float64
ε = 0.005::Float64

# Define the fast and slow variable vector field
f(x, μ) = -μ - 2*x + 3*(x^2) - (4/5)*(x^3)
g(t) = ε
η(x) = σ

# Define the IC
x0 = 2.789::Float64
μ0 = 0.400::Float64
u0 = [x0, μ0]

# Solve the fast-slow SDE
printstyled("Solving the fast-slow saddle-node SDE\n"; bold=true, underline=true, color=:light_blue)
t, μ, u = evolve_ramped_1d(f, g, η, u0, δt=1e-2, μf=1.400)

# Compute the stable equilibrium at each parameter value
X0 = [sort(get_equilibria(f, p, domain=[-20,20])[1], rev=true)[1] for p in μ]

# Get number of timesteps in the ramped solution
Nμ = length(μ)

# Define number of timesteps for the OUP equivalent
Nt = convert(Int64,1e3)

# Define matrix to store the parameter sweep solution of the OUP equivalent
U = Matrix{Float64}(undef, Nt+1, Nμ)

# Define vector to store the linearisations of the saddle-node
ϴ = Vector{Float64}(undef, Nμ)

# Loop over the parameter values
printstyled("Solving the OUP equivalent parameter sweep using $(Threads.nthreads()) threads\n"; bold=true, underline=true, color=:light_blue)
@showprogress Threads.@threads for n in 1:Nμ
        # Compute the linearised saddle-node system around the stable equilibrium
        F(x) = f(x, μ[n])
        ϴ[n] = ForwardDiff.derivative(F, X0[n])

        # Define a topologically equivalent OUP
        h(x, y) = ϴ[n]*(x - y)

        # Solve the OUP
        t_OUP, u_OUP = evolve_forward_1d(h, η, X0[n], IC=X0[n], Nt=Nt)

        # Export the solution
        U[:,n] = u_OUP[:,1]
end

# Export the data
writeout(hcat(μ, ϴ, X0), "../data/slides_10_to_13/parameters.csv")
writeout(hcat(t, μ, u), "../data/slides_10_to_13/solution.csv")
writeout(U, "../data/slides_10_to_13/OUP.csv")

# Execute the postprocessing and plotting scripts:
# Slide 10
include("../postprocessing/slide_10_postprocessing.jl")
include("../plotting/slide_10_plotting.jl")
# Slide 11
include("../postprocessing/slide_11_postprocessing.jl")
include("../plotting/slide_11_plotting.jl")
# Slide 12
include("../postprocessing/slide_12_postprocessing.jl")
include("../plotting/slide_12_plotting.jl")
# Slide 13
include("../postprocessing/slide_13_postprocessing.jl")
include("../plotting/slide_13_plotting.jl")
