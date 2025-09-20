"""
Dynamics utility functions associated to propagating a solution of a system of differential equations forward in time. 

Author: Davide Papapicco
Affil: U. of Auckland
Date: 21-08-2025
"""

#--------------------------#
#                          # 
#   evolve_1d_systems.jl   #                     
#                          #
#--------------------------#

"""
$(TYPEDSIGNATURES)

Solve an ensemble of stationary IVPs in 1-dimension.

The ensemble problem is made of i.i.d. SDEs propagating forward in time from a unique, deterministic initial condition `u0`.
The SDEs are characterised by a parametric drift `f` (evaluated at fixed parameter value `μ`) and diffusion `η`.

## Keyword arguments
* `δt=1e-2`: timestep of the solver
* `saveat=δt`: timestep at which solutions are exported
* `Nt=1000::Int64`: total number of timesteps
* `Ne=1::Int64`: total number of particles in the ensemble

## Output
`solutions::Tuple`
* `solutions.time::Vector{Float64}`: timestamps of the trajectories
* `solutions.states::Vector{Vector{Float64}}`: trajectories of the ensemble 

## Example
"""
function evolve_1d(f::Function, η::Function, μ, u0::Vector; δt=1e-2, saveat=δt, Nt=1000::Int64, Ne=1::Int64)
        # Compute the total time 
        T = δt*Nt

        # Deterministic dynamics 
        function drift!(dudt, u, μ, t)
                dudt[1] = f(u[1], μ[1]) 
        end
        # Stochastic dynamics
        function diffusion!(dudt, u, μ, t)
                dudt[1] = η(u[1]) 
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
        states = Vector{Vector{Float64}}()
        for n in 1:Ne
                push!(states, [(sol[n].u)[t][1] for t in 1:Nt])
        end

        # Return the solutions and their timestamps
        return (
                time = time, 
                states = states
               )
end


"""
$(TYPEDSIGNATURES)

Solve an ensemble of non-autonomous IVPs in 1-dimensions.

The system of SDEs is non-autonomous in the sense that the drift `f` has now explicit time dependence.
Such time-dependence is encoded into the parameter shift `Λ` which turns the problem into a 2-dimensional one.
As such the user needs to provide `u0 = [x0,λ0]` as an initial condition.

The problem can either be solved on a time interval `T` or on a parameter interval `[u0[2],λf]`, in which case the user provides `λf::Float64` instead of `T::Vector{Float64}`.

## Keyword arguments
* `saveat=δt`: timestep at which solutions are exported
* `Nt=1000::Int64`: total number of timesteps
* `Ne=1::Int64`: total number of particles in the ensemble

## Output
`solutions::Tuple`
* `solutions.time::Vector{Float64}`: timestamps of the trajectories
* `solutions.param::Vector{Float64}`: solution of the parameter shift 
* `solutions.states::Vector{Vector{Float64}}`: trajectories of the ensemble 

## Example
"""
function evolve_shifted_1d(f::Function, Λ::Function, η::Function, u0::Vector, T::Vector; Nt=1000::Int64, saveat=nothing)
        # Compute the timestep
        δt = (T[end]-T[1])/Nt

        # Compute the export timestep
        Δt = 0.00::Float64
        if saveat == nothing
                Δt = δt
        else
                Δt = saveat 
        end

        # Deterministic dynamics 
        function drift!(dudt, u, λ, t)
                dudt[1] = f(u[1], u[2])
                dudt[2] = λ(t)
        end
        # Stochastic dynamics
        function diffusion!(dudt, u, λ, t)
                dudt[1] = η(u[1])
                dudt[2] = 0 
        end
        
        # Define the differential equation
        dynamics = SDEProblem(drift!, diffusion!, u0, (T[1], T[end]), Λ)

        # Solve the SDE forward in time
        sol = solve(dynamics, EM(), dt=δt, verbose=false, saveat=Δt)

        # Extract the timestamps
        time = sol.t
        Nt = length(time)

        # Extract the solution and the time evolution
        states = [(sol.u)[t][1] for t in 1:Nt]
        param = [(sol.u)[t][2] for t in 1:Nt]

        # Return the solutions and their timestamps
        return (
                time = time, 
                param = param, 
                states = states
               )
end

function evolve_shifted_1d(f::Function, Λ::Function, η::Function, u0::Vector, λf::Float64; δt=1e-2, saveat=δt, Ne=1::Int64)
        # Specify the simulation endtime and number of timesteps according to the parameter range
        T = (μf-u0[2])/Λ(0)
        Nt = T/δt

        # Deterministic dynamics 
        function drift!(dudt, u, λ, t)
                dudt[1] = f(u[1], u[2])
                dudt[2] = Λ(u[2])
        end
        # Stochastic dynamics
        function diffusion!(dudt, u, λ, t)
                dudt[1] = η(u[1])
                dudt[2] = 0 
        end
        
        # Define the differential equation
        dynamics = SDEProblem(drift!, diffusion!, u0, (0.0, T))

        # Convert it into an ensemble problem
        ensemble = EnsembleProblem(dynamics)

        # Solve the problem forward in time
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
