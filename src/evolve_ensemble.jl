include("./get_equilibria.jl")

# Solve the ensemble problem
using DifferentialEquations
function evolve_ensemble(f::Function, g::Function, μ::Float64; IC=nothing, δt=1e-2, Nt=1000::Int64, Ne=1000::Int64)
        # Compute the total time 
        T = δt*Nt

        # Deterministic dynamics 
        function drift!(dudt, u, μ, t)
                dudt[1] = f(u[1], μ[1]) 
        end
        # Stochastic dynamics
        function diffusion!(dudt, u, μ, t)
                dudt[1] = g(u[1]) 
        end

        # Specify the IC depending on the user input
        u0 = [0.0]
        if IC == nothing
                # Perturbed from a stable equilibrium 
                equilibria = get_equilibria(f, μ)
                stable = equilibria[1]
                u0 = [stable[1] + σ/20]
        else
                # Input from the user
                u0 = [IC]
        end
       
        # Define the differential equation
        dynamics = SDEProblem(drift!, diffusion!, u0, (0.0, T), μ)

        # Convert it into an ensemble problem
        ensemble = EnsembleProblem(dynamics)

        # Solve the problem forward in time
        sol = solve(ensemble, EM(), dt=δt, verbose=false, EnsembleDistributed(), trajectories=Ne)

        # Extract the timepoints for plotting purposes
        time = sol[1].t
        Nt = length(time)

        # Extract the timeseries
        states = Matrix{Float64}(undef, Ne, Nt) 
        for n in 1:Ne
                states[n,:] = [(sol[n].u)[t][1] for t in 1:Nt]
        end

        # Return the solutions and their timestamps
        return time, states
end
