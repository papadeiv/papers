# Solve the 1-dimensional non-autonomous SDE forward in time for a linearly ramped parameter

function evolve_ramped_1d(f::Function, g::Function, η::Function, u0; δt=1e-2, Nt=1000::Int64, saveat=δt, μf=nothing)
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

        # Solve the SDE forward in time
        printstyled("Solving the non-autonomous SDE\n"; bold=true, underline=true, color=:light_blue)
        sol = solve(dynamics, EM(), dt=δt, verbose=false, saveat=saveat)

        # Extract the timestamps
        time = sol.t
        Nt = length(time)

        # Extract the solution and the time evolution
        states = [(sol.u)[t][1] for t in 1:Nt]
        param = [(sol.u)[t][2] for t in 1:Nt]

        # Return the solutions and their timestamps
        return time, param, states
end
