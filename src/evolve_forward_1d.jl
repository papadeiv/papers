include("./get_equilibria.jl")

# Solve the SDE forward in time for a fixed parameter value 
using DifferentialEquations
function evolve_forward_1d(f::Function, g::Function, μ::Float64; IC=nothing, δt=5e-2, Nt=1000::Int64, saveat=δt)
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

        # Solve the SDE forward in time
        sol = solve(dynamics, EM(), dt=δt, verbose=false, saveat=saveat)

        # Extract the timestamps of the (saved) states
        time = sol.t
        Nt = length(time)

        # Extract the solution and timestamps
        states = [(sol.u)[t][1] for t in 1:Nt]

        # Return the solutions and their timestamps
        return time, states
end
