"""
Dynamics utility functions associated to propagating a solution of a system of differential equations forward in time. 

Author: Davide Papapicco
Affil: U. of Auckland
Date: 21-08-2025
"""

#----------------------#
#                      # 
#   evolve_system.jl   #                     
#                      #
#----------------------#

"""
    evolve_ensemble()

Solves the SDE with drift `f` and diffusion `g` for an ensemble of Ne independent particles starting (deterministically) from `x0`.
"""
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
      
        # Define the differential equation
        dynamics = SDEProblem(drift!, diffusion!, IC, (0.0, T), μ)

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
