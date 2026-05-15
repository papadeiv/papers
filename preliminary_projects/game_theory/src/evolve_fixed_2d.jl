# Solve the parameter sweep of a 2-dimensional autonomous SDE

function evolve_fixed_2d(f1::Function, f2::Function, g1::Function, g2::Function, μ::Vector{Float64}, IC::Vector{Float64}; δt=5e-2, Nt=1000::Int64, Nμ=nothing)
        # Compute the total time 
        T = δt*Nt

        # Define the vector of fixed parameter values according to the user's choice
        if Nμ == nothing
                # The parameter values are already specified by the user
                global parameters = μ
        else
                # Generate uniformly spaced parameter values within a specified range
                global parameters = LinRange(μ[1], μ[end], Nμ)
        end

        # Define matrix to store the solutions and their timestamps at each parameter value
        time = Vector{Float64}(undef, Nt+1)
        u1 = Matrix{Float64}(undef, Nt+1, length(parameters))
        u2 = Matrix{Float64}(undef, Nt+1, length(parameters))

        # Loop over the parameter values
        for n in 1:length(parameters)
                # Solve the SDE for the current parameter value
                t, u = evolve_forward_2d(f1, f2, g1, g2, parameters[n], IC, δt=δt, Nt=Nt)

                # Export the timestamps only at the first parameter value
                n==1 && (time = t)

                # Export the solutions
                u1[:,n] = u[:,1]
                u2[:,n] = u[:,2]
        end

        # Concatenate the solutions in the 2 scalar variables
        states = [u1, u2]

        # Return the solutions, their timestamps and parameter values
        return time, parameters, states
end
