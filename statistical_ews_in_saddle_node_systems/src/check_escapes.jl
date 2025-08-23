# Assemble the escape time distribution

function check_escapes(u, boundaries)
        # Get the number of realizations
        Nt = length(u)

        # Initialise the time index of first escape
        Nt_exit = Nt

        # Initialise the boolean variable for escapes
        check = false

        # Loop over the timesteps
        for n in 1:Nt
                # Check if the trajectory has crossed the boundaries of the basin at least once (strict inequalities here are important)
                if u[n] < boundaries[1] || u[n] > boundaries[2]
                        check = true
                        Nt_exit = n-1
                        break
                end
        end

        # Return the boolean variable and the time index of first escape
        return check, Nt_exit 
end
