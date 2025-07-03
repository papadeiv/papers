include("../../../../inc/PotentialLearning.jl")
include("../../../../inc/SystemAnalysis.jl")
include("../../../../inc/IO.jl")

# Define the saddle-node SDE 
f(x, μ) = -μ - x^2
σ = 0.200::Float64 
η(x) = σ

# Define the parameter range for the sweep
Nμ = convert(Int64,4e2)
μ_min = -(1.000::Float64)
μ_max = -(0.001::Float64)
μ0 = LinRange(μ_min, μ_max, Nμ)

# Compute the stable equilibrium at each parameter value
x0 = [sort(get_equilibria(f, μ, domain=[-20,20])[1], rev=true)[1] for μ in μ0]

# Solve the parameter sweep of the deterministic case
Nt = convert(Int64,1e3)
t, μ, u = evolve_fixed_1d(f, η, collect(μ0), x0, Nt=Nt)

# Define matrix to store the parameter sweep solution of the OUP equivalent
U = Matrix{Float64}(undef, Nt+1, Nμ)

# Define vector to store the linearisations of the saddle-node
ϴ = Vector{Float64}(undef, Nμ)

# Loop over the parameter values
printstyled("Solving the OUP equivalent parameter sweep using $(Threads.nthreads()) threads\n"; bold=true, underline=true, color=:light_blue)
@showprogress Threads.@threads for n in 1:Nμ
        # Compute the linearised saddle-node system around the stabloe equilibrium
        F(x) = f(x, μ0[n])
        ϴ[n] = ForwardDiff.derivative(F, x0[n])

        # Define a topologically equivalent OUP
        g(x, y) = ϴ[n]*(x - y)

        # Solve the OUP
        t, u = evolve_forward_1d(g, η, x0[n], IC=x0[n], Nt=Nt)

        # Export the solution
        U[:,n] = u[:,1]
end

# Export the data
writeout(hcat(t, μ, ϴ, x0), "../data/slide_07/parameters.csv")
writeout(u, "../data/slide_07/saddle_node.csv")
writeout(U, "../data/slide_07/OUP.csv")

# Execute the postprocessing and plotting scripts
include("../postprocessing/slide_07_postprocessing.jl")
include("../plotting/slide_07_plotting.jl")
