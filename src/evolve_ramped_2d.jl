# Solve the 2-dimensional non-autonomous SDE forward in time for a linearly ramped parameter

function evolve_ramped_2d(f1::Function, f2::Function, g1::Function, g2::Function, h::Function, u0; δt=5e-2, Nt=1000::Int64, saveat=δt, μf=nothing)
        # Compute the temporal quantities according to the user's choice
        if μf==nothing
                # Specify the simulation endtime according to the number of steps
                global T = δt*Nt
        else
                # Specify the simulation endtime according to the parameter range
                global T = (μf-u0[3])/h(0)
        end

        # Deterministic dynamics 
        function drift!(dudt, u, μ, t)
                dudt[1] = f1(u[1], u[2], u[3])
                dudt[2] = f2(u[1], u[2], u[3])
                dudt[3] = h(u[3])
        end
        # Stochastic dynamics
        function diffusion!(dudt, u, μ, t)
                dudt[1] = g1(u[1], u[2])
                dudt[2] = g2(u[1], u[2])
                dudt[3] = 0 
        end
        
        # Define the differential equation
        dynamics = SDEProblem(drift!, diffusion!, u0, (0.0, T))

        # Solve the SDE forward in time
        sol = solve(dynamics, EM(), dt=δt, verbose=false, saveat=saveat)

        # Extract the timestamps
        t = sol.t
        Nt = length(t)

        # Extract the solution and the time evolution
        x1 = [(sol.u)[t][1] for t in 1:Nt]
        x2 = [(sol.u)[t][2] for t in 1:Nt]
        μ = [(sol.u)[t][3] for t in 1:Nt]

        # Return the solutions and their timestamps
        return t, μ, hcat(x1, x2) 
end
