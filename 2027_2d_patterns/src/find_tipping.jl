# Identify the tipping point of a critical trasition in a timeseries

function find_tipping(ut::Vector{Float64}, width, check)
        # Get the length of the timeseries
        Nt = length(ut)
        
        # Assemble the sliding window
        window = get_window_parameters(Nt, width)
        Nw = window[1]
        Ns = window[2]

        # Define the split point in the sliding window to check the tipping point
        split = get_window_parameters(Nw, check)
        Np = split[1]

        #println("Nw = ", Nw, ", Np = ", Np)

        # Initialise the tipping point flags
        tip_chk = false 
        tip_idx = Nt

        # Loop over the strides of the sliding window
        printstyled("Searching for tipping points\n"; bold=true, underline=true, color=:light_blue)
        @showprogress for n in 1:Ns
                # Define the index of the timestep at the end and split end of the sliding window
                n_end = Nw + n - 1::Int64
                n_split = n_end - Np 

                # Extract the subseries in the window
                u_presplit = ut[n:n_split]
                u_postsplit = ut[(n_split+1):n_end]

                # Compute the mean and std in the presplit subseries
                u_mean = mean(u_presplit)
                u_std = std(u_presplit)

                # Check whether the realisations postsplit have all tipped
                u_inf = u_mean - 3*u_std
                u_sup = u_mean + 3*u_std
                if !any(u -> u_inf <= u <= u_sup, u_postsplit)
                        # Store the tip flags and stop the loop
                        tip_chk = true
                        tip_idx = n_split
                        break
                end
        end

        if tip_chk
                # Print this message if tipping has been bound
                printstyled("Tipping point found at timestep ", tip_idx, " of ", Nt,"\n"; bold=true, underline=true, color=:green)
          
        else
                # Print this message if not tipping has been bound
                printstyled("No tipping point found\n"; bold=true, underline=true, color=:red)
        end

        # Return the boolean acheck and the tip index
        return tip_chk, tip_idx
end
