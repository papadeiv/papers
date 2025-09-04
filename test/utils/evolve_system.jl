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

"""
    evolve_ramped_1d()

Solves an ensemble of slow-fast SDEs with drift `f`, ramp `g` and diffusion `η`.
"""
function evolve_ramped_1d(f::Function, g::Function, η::Function, u0; δt=1e-2, Nt=1000::Int64, saveat=δt, μf=nothing, Ne=1::Int64)
        # Define endtime value
        T = 0.0::Float64

        # Compute the temporal quantities according to the user's choice
        if μf==nothing
                # Specify the simulation endtime according to the number of steps
                T = δt*Nt
        else
                # Specify the simulation endtime according to the parameter range
                T = (μf-u0[2])/g(0)
                Nt = T/δt
        end

        # Deterministic dynamics 
        function drift!(dudt, u, μ, t)
                dudt[1] = f(u[1], u[2])
                dudt[2] = g(u[2])
        end
        # Stochastic dynamics
        function diffusion!(dudt, u, μ, t)
                dudt[1] = η(u[1])
                dudt[2] = 0 
        end
        
        # Define the differential equation
        dynamics = SDEProblem(drift!, diffusion!, u0, (0.0, T))

        # Convert it into an ensemble problem
        ensemble = EnsembleProblem(dynamics)

        # Solve the problem forward in time
        printstyled("Solving ", Ne, " slow-fast SDEs\n"; bold=true, underline=true, color=:light_blue)
        sol = solve(ensemble, EM(), dt=δt, verbose=false, EnsembleDistributed(), trajectories=Ne)

        # Extract the timestamps and parameter values
        time = sol[1].t
        Nt = length(time)
        param = [(sol[1].u)[t][2] for t in 1:Nt]

        # Extract the timeseries
        states = Vector{Vector{Float64}}() 
        for n in 1:Ne
                push!(states, [(sol[n].u)[t][1] for t in 1:Nt])
        end

        # Return the solutions and their timestamps
        return (
                time = time,
                parameter = param,
                states = states
               )
end
