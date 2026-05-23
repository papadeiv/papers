"""
    Main script

Run this file to execute the simulation, analyse and plot the results.
"""

# Import the necessary packages and local modules
include("inc.jl")

# Import the simulation's scripts
include("./scripts/sim.jl")
include("./scripts/figs.jl")
include("./scripts/proc.jl")

# Define the main algorithm
function main(iteration::Integer, case::Integer)
        # Solve the ensemble problem
        ensemble = evolve(f, η, Λ, [a, α[case]], stepsize=dt, steps=Nt, particles=Ne)

        # Loop over the solutions in the ensemble
        solutions = Matrix{Float64}(undef, convert(Integer, Ne), 3)
        @showprogress for (index, sample) in enumerate(ensemble.state)
                # Suppress the output from Python
                old_stdout = pyimport("sys").stdout
                old_stderr = pyimport("sys").stderr
                devnull = pyimport("io").StringIO()
                pyimport("sys").stdout = devnull
                pyimport("sys").stderr = devnull

                # Solve the MLE problem for the given sample
                @pyinclude("./scripts/mle.py")
                solutions[index,:] = py"infer"(sample, dt)

                pyimport("sys").stdout = old_stdout
                pyimport("sys").stderr = old_stderr
        end

        # Export the solutions
        writeout(solutions, "α=$(α[case])/$iteration.csv")

        # Empty the ensemble solutions
        solutions = nothing
        ensemble = nothing
end

# Loop over the number of cases
for m in 1:length(α)
        # Plot the potential
        domain = LinRange(-0.5,0.5,1000)
        lines!(ax_p[m], domain, [(α[m]/2)*x^2 for x in domain], color = :black, linewidth = 3.0)

        # Loop over the number of executions
        for n in 1:N_exec
                # Execute the main
                printstyled("Case α = $(α[m]): iteration n.$n of $N_exec\n"; bold=true, underline=true, color=:light_blue)
                main(n, m)

                # Force garbage collection
                GC.gc()
        end

        # Analyse and plot the results
        θ = [α[m], a, σ]
        analyse(m, θ)
end

# Export the figure
savefig("oup_mle.pdf", fig)
