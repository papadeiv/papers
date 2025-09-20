"""
Dynamics utility functions associated to propagating a solution of a system of differential equations forward in time. 

Author: Davide Papapicco
Affil: U. of Auckland
Date: 21-08-2025
"""

#--------------------------#
#                          # 
#   evolve_2d_systems.jl   #                     
#                          #
#--------------------------#

"""
$(TYPEDSIGNATURES)

Solve an ensemble of stationary IVPs in 2-dimensions.

The ensmeble problem is made of i.i.d. SDEs propagating forward in time from a unique, deterministic initial condition `u0`.
The SDEs are characterised by a parametric drift `f1`,`f2`, corresponding to the x and y component respectively, and evaluated at fixed parameter value `μ`.
The diffusion is `η1`,`η2` for the x and y component respectively.

## Keyword arguments
* `δt=1e-2`: timestep of the solver
* `Nt=1000::Int64`: total number of timesteps
* `Ne=1::Int64`: total number of particles in the ensemble

## Output
`solutions::Tuple`
* `solutions.time::Vector{Float64}`: timestamps of the trajectories
* `solutions.states::Vector{Matrix{Float64}}`: trajectories of the ensemble arranged in a `Nx2` matrix where column 1 is the timeseries of the x component and column 2 is the trajectory of the y component 

## Example
"""
function evolve_2d(f1::Function, f2::Function, η1::Function, η2::Function, μ::Float64, u0::Vector{Float64}; δt=1e-2, Nt=1000::Int64)
        # Compute the total time 
        T = δt*Nt

        # Deterministic dynamics 
        function drift!(dudt, u, μ, t)
                dudt[1] = f1(u[1], u[2], μ[1])
                dudt[2] = f2(u[1], u[2], μ[1])
        end
        # Stochastic dynamics
        function diffusion!(dudt, u, μ, t)
                dudt[1] = η1(u[1], u[2])
                dudt[2] = η2(u[1], u[2])
        end
        
        # Define the differential equation
        dynamics = SDEProblem(drift!, diffusion!, u0, (0.0, T), μ)

        # Solve the SDE forward in time
        sol = solve(dynamics, EM(), dt=δt, verbose=false)

        # Extract the solution and timestamps
        time = sol.t
        Nt = length(time)
        u1 = [(sol.u)[t][1] for t in 1:Nt]
        u2 = [(sol.u)[t][2] for t in 1:Nt]
        states = hcat(u1, u2)

        # Return the solutions and their timestamps
        return (
                time = time, 
                states = states
               )
end

"""
$(TYPEDSIGNATURES)

Solve an ensemble of non-autonomous IVPs in 2-dimensions.

The system of SDEs is non-autonomous in the sense that the drift `f1`,`f2`, for the x and y component respectively, has explicit time dependence through the parameter shift `Λ`.
As such the user needs to provide `u0 = [x0,y0,λ0]` as an initial condition as the problem is solved on a time interval `T::Vector{Float64}`.
The diffusion is `η1`,`η2` for the x and y component respectively.

## Keyword arguments
* `Nt=1000::Int64`: total number of timesteps
* `Ne=1::Int64`: total number of particles in the ensemble

## Output
`solutions::Tuple`
* `solutions.time::Vector{Float64}`: timestamps of the trajectories
* `solutions.param::Vector{Float64}`: solution of the parameter shift 
* `solutions.states::Vector{Matrix{Float64}}`: trajectories of the ensemble arranged in a `Nx2` matrix where column 1 is the timeseries of the x component and column 2 is the trajectory of the y component 

## Example
"""
function evolve_shifted_2d(f1::Function, f2::Function, Λ::Function, η1::Function, η2::Function, u0, T; Nt=1000::Int64, saveat=nothing)
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
                dudt[1] = f1(u[1], u[2], u[3])
                dudt[2] = f2(u[1], u[2], u[3])
                dudt[3] = λ(t)
        end
        # Stochastic dynamics
        function diffusion!(dudt, u, λ, t)
                dudt[1] = η1(u[1], u[1])
                dudt[2] = η2(u[1], u[1])
                dudt[3] = 0 
        end
        
        # Define the differential equation
        dynamics = SDEProblem(drift!, diffusion!, u0, (T[1], T[end]), Λ)

        # Solve the SDE forward in time
        sol = solve(dynamics, EM(), dt=δt, verbose=false, saveat=Δt)

        # Extract the timestamps
        time = sol.t
        Nt = length(time)

        # Extract the solution and the time evolution
        param = [(sol.u)[t][3] for t in 1:Nt]
        u1 = [(sol.u)[t][1] for t in 1:Nt]
        u2 = [(sol.u)[t][2] for t in 1:Nt]
        states = hcat(u1, u2)

        # Return the solutions and their timestamps
        return (
                time = time, 
                param = param, 
                states = states
               )
end
