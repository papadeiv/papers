# Solve the 1-dimensional non-autonomous SDE in compactified time for a bounded parameter shift 

function evolve_shifted_1d(f::Function, Λ::Function, η::Function, u, T; Nt=1000::Int64, saveat=nothing)
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
        dynamics = SDEProblem(drift!, diffusion!, u, (T[1], T[end]), Λ)

        # Solve the SDE forward in time
        printstyled("Solving the non-autonomous SDE\n"; bold=true, underline=true, color=:light_blue)
        sol = solve(dynamics, EM(), dt=δt, verbose=false, saveat=Δt)

        # Extract the timestamps
        time = sol.t
        Nt = length(time)

        # Extract the solution and the time evolution
        states = [(sol.u)[t][1] for t in 1:Nt]
        param = [(sol.u)[t][2] for t in 1:Nt]

        # Return the solutions and their timestamps
        return time, param, states
end
