# Derive the distribution of the stochastic realizations on the potential
function invert_equilibrium_distribution(bins, distribution, noise::Float64; N = nothing)
        # Compute the diffusion coefficient
        D = (noise^2)/2

        # Filter out the 0-valued entries in the distribution
        idx = findall(x -> x > 0.0, distribution)
        ys = [distribution[n] for n in idx]
        xs = [bins[n] for n in idx]

        # Define the normalisation constant based on user input
        if N == nothing
                # Assumption of a stationary OUP
                N=1/sqrt(4*pi*D)
        end

        # Compute the distribution on the potential using the stationarity assumption
        Vs = -D.*log.(ys./N)

        # Return the filtered datapoints and their location
        return xs, Vs
end
