include("../../../../inc/SystemAnalysis.jl")
include("../../../../inc/IO.jl")

# Number of timesteps of each particle
Nt = convert(Int64,1e5)
# Number of (fixed) parameter values
Nμ = convert(Int64,1e1)

# Define the normal form of the monostable saddle-node
f(x, μ) = -μ - x^2 
# Define the normal form of the topologically equivalent OUP
h(x, ϴ) = -ϴ*x

# Specify the (additive) noise level
σ = 0.100
g(x) = σ

# Specify the parameter range and export it 
μ = LinRange(-2.000,-0.200,Nμ)
writeout(μ, "../data/parameter.csv")

# Loop over the parameter range 
printstyled("Simulating the SDE (monostable saddle node)\n"; bold=true, underline=true, color=:light_blue)
@showprogress for n in 1:Nμ
        # Propagate the saddle-node normal form forward in time  
        t, u = evolve_forward_1d(f, g, μ[n], Nt=Nt)

        # Export the data 
        writeout(hcat(t, u), "../data/solutions/monostable_$n.csv")

        # Compute the parameter of the topologically-equivalent OUP 
        eq = get_equilibria(f, μ[n])
        local ϴ = abs(-2*eq[1][1])

        # Propagate the OUP forward in time
        t, u = evolve_forward_1d(h, g, ϴ, Nt=Nt)

        # Export the data
        writeout(hcat(t, u), "../data/solutions/OUP_$n.csv")
end

# Execute the additional scripts for this simulation
include("../postprocessing/01_postprocessing.jl")
include("../plotting/01_plotting.jl")
