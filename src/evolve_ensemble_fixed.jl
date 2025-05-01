include("./evolve_ensemble.jl")

# Solve the ensemble problem for a range of parameter values
using ProgressMeter
function evolve_ensemble_fixed(f::Function, g::Function, μ::Vector{Float64}, IC::Vector{Float64}; δt=1e-2, Nt=1000::Int64, Nμ=nothing, Ne=1000::Int64)
        # Define the vector of fixed parameter values according to the user's choice
        if Nμ == nothing
                # The parameter values are already specified by the user
                global parameters = μ
        else
                # Generate uniformly spaced parameter values within a specified range
                global parameters = LinRange(μ[1], μ[end], Nμ)
        end

        # Get the number of parameter values in the range
        Nμ = length(parameters)

        # Define arrays to store the ensemble solutions and their timestamps
        t_ensemble = Vector{Float64}(undef, Nt+1)
        u_ensemble = Vector{Matrix{Float64}}(undef, Nμ)

        # Loop over the parameter values
        printstyled("Solving the ensemble problems in the prescribed parameter range\n"; bold=true, underline=true, color=:light_blue)
        @showprogress for n in 1:length(parameters)
                # Solve the ensemble problem for the current parameter value
                t, u = evolve_ensemble(f, g, parameters[n], IC=IC[n], δt=δt, Nt=Nt, Ne=Ne)

                # Export the timestamps only at the first parameter value
                n==1 && (t_ensemble = t)

                # Export the ensemble solution
                u_ensemble[n] = u
        end

        # Return the results
        return t_ensemble, parameters, u_ensemble
end
