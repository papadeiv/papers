using ProgressMeter
using DifferentialEquations

# Solve the 1-dimensional SDE forward in time for fixed parameter values in a range 
function evolve_fixed_1d(f::Function, g::Function, μ::Vector{Float64}, IC::Vector{Float64}; δt=5e-2, Nt=1000::Int64, Nμ=nothing)
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
        time = Vector{Int64}(undef, length(parameters))
        states = Matrix{Float64}(undef, Nt+1, length(parameters))

        # Loop over the parameter values
        printstyled("Solving the SDE parameter sweep using $(Threads.nthreads()) threads\n"; bold=true, underline=true, color=:light_blue)
        @showprogress Threads.@threads for n in 1:length(parameters)
                # Solve the SDE for the current parameter value
                t, u = evolve_forward_1d(f, g, parameters[n], IC=IC[n], δt=δt, Nt=Nt)

                # Get the total number of realisations for the current parameter
                time[n] = length(u[:,1])

                # Export the solutions
                states[1:time[n],n] = u[:,1]
        end

        # Return the solutions, their timestamps and parameter values
        return time, parameters, states 
end
