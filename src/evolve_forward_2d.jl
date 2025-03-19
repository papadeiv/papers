# Solve the 2-dimensional SDE forward in time for a fixed parameter value 
using DifferentialEquations
function evolve_forward_2d(f1::Function, f2::Function, g1::Function, g2::Function, μ::Float64, IC::Vector{Float64}; δt=5e-2, Nt=1000::Int64)
        # Compute the total time 
        T = δt*Nt

        # Deterministic dynamics 
        function drift!(dudt, u, μ, t)
                dudt[1] = f1(u[1], u[2], μ[1])
                dudt[2] = f2(u[1], u[2], μ[1])
        end
        # Stochastic dynamics
        function diffusion!(dudt, u, μ, t)
                dudt[1] = g1(u[1], u[2])
                dudt[2] = g2(u[1], u[2])
        end
        
        # Define the differential equation
        dynamics = SDEProblem(drift!, diffusion!, IC, (0.0, T), μ)

        # Solve the SDE forward in time
        sol = solve(dynamics, EM(), dt=δt, verbose=false)

        # Extract the solution and timestamps
        u1 = [(sol.u)[t][1] for t in 1:(Nt+1)]
        u2 = [(sol.u)[t][2] for t in 1:(Nt+1)]
        states = hcat(u1, u2)
        time = sol.t

        # Return the solutions and their timestamps
        return time, states
end
