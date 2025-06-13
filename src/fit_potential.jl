# Find optimal coefficients of the potential to fit the observed distribution 

function fit_potential(xn, yn, σ; initial_guess=nothing)
        # Define the stochastic diffusion
        D = (σ^2)/2.0::Float64

        # Compute a shift for the potential {c0} that sets V(xs)=0 to avoid numerical cancellation
        xs(μ) = (1/(3*μ[3]))*(sqrt((μ[2])^2 - 3*μ[1]*μ[3]) - μ[2])
        c0(μ) = - μ[1]*xs(μ) - μ[2]*(xs(μ))^2 - μ[3]*(xs(μ))^3

        # Define an arbitrary cubic with the the above constraint on {c0}
        V(x, μ) = c0(μ) + μ[1]*x + μ[2]*(x^2) + μ[3]*(x^3)
        # Define the unnormalised pdf as an exponential of the abritrary cubic
        f(x, μ) = exp(-(1.0::Float64/D)*(V(x, μ)))

        # Define the normalisation constant as a function of the 4 parameters
        N(μ) = get_normalisation_constant(f, (-(1/(3*μ[3]))*(sqrt((μ[2])^2 - 3*μ[1]*μ[3]) + μ[2]), Inf::Float64), parameters=μ)

        # Define the target of the optimisation problem: normalised pdf
        p(x, μ) = N(μ)*exp.(-(1.0::Float64/D).*(c0(μ) .+ μ[1]*x .+ μ[2]*(x.^2) .+ μ[3]*(x.^3)))

        # Initial guess for the optimal parameters
        if initial_guess == nothing
                xs, Vs = invert_equilibrium_distribution(xn, yn, σ)
                initial_guess = approximate_potential(xs, Vs, degree=3).coeffs[2:end] 
        end

        # Define the lower and upper bounds for the coefficients
        lower = [-45.0, -25.0, 0.01] 
        upper = [45.0, 25.0, 4.0]

        # Get the non-linear least-squares solution
        c = curve_fit(p, xn, yn, initial_guess, lower=lower, upper=upper).param

        # Return the solution
        return c 
end
