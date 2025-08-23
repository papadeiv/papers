# Assemble the escape time distribution

function ensemble_escapes(time, ensemble, equilibria)
        # Get the number of realizations
        Nt = length(ensemble[1,:])
        # Get the number of sample paths
        Ne = length(ensemble[:,1])

        # Loop over the ensemble
        hitting = Float64[]
        for n in 1:Ne
                # Get the n-th trajectory
                ut = ensemble[n,:] 
                # Initialise a placeholder variable for the escape time
                τ = 0
                # Initialise a boolean variable for the escape event
                escaped = false
                # Loop over its states sequence
                for t in 1:Nt
                        # Enforce absorbing BCs
                        if ut[t][1] > equilibria[2] 
                                τ = time[t]
                                escaped = true
                                break
                        end
                end
                # Update the hitting time distribution
                if escaped
                        push!(hitting, τ)
                end
        end

        # Assemble the cumulative distribution of escapes at each time
        distribution = sort(hitting) 
        return distribution 
end
