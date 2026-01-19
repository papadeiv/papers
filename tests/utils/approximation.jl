"""
Approximation methods to reconstruct the scalar potential of a conservative system from the empirical distribution of a sample path.

Author: Davide Papapicco
Affil: U. of Auckland
Date: 26-09-2025
"""

"""
$(TYPEDSIGNATURES)

Computes the approximate normalisation constant of a pdf by Gauss-Kronrod quadrature.

The input function `f::Function` depends on some `parameters::Vector{Float64}` from which the approximate normalisation constant is computed over a specified `domain::Tuple`.
If no `domain` is given, then `f::Function` is assumed to depend on a cubic potential of coefficients `{c1,c2,c3}` provided in `parameters`.

## Keyword Arguments
* `accuracy=1e-8`: relative accuracy criterion for the quadrature iterations.

## Output
* `normalisation_constant::Float64`: normalisation constant to make `f` a pdf 

## Example
"""
function normalise(f::Function, parameters, domain; accuracy=1e-8)
        # Define the integral problem over the domain
        if parameters == nothing
                global integral = IntegralProblem(f, domain)
        else
                global integral = IntegralProblem(f, domain, parameters)
        end

        # Solve the definite integral by using adaptive Gauss-Kronrod quadrature
        quadrature = solve(integral, QuadGKJL(); maxiters=10000, reltol=accuracy, abstol=accuracy)

        # Return the approximation
        return 1.0::Float64/(quadrature.u)
end

function normalise(f::Function, parameters; accuracy=1e-12)
        # Compute the location of the local minima and maxima
        μ = parameters
        xs = +(1/(3*μ[3]))*(sqrt((μ[2])^2 - 3*μ[1]*μ[3]) - μ[2])
        xu = -(1/(3*μ[3]))*(sqrt((μ[2])^2 - 3*μ[1]*μ[3]) + μ[2])

        # Define the integration interval
        I = (-Inf,Inf)

        # Initialise the interval
        if xs > xu
                # Avoid numerical cancellation in the quadrature
                if abs(xu-xs) < 1e1
                        I = (xu, +Inf) 
                else
                        # Find a small enough value of the integrand
                        subinterval = collect(LinRange(xs,xu,100000))
                        idx = findfirst(p -> p < 1e-16, [f(x, μ) for x in subinterval])
                        xu = subinterval[idx]
                        I = (xu, +Inf)
                end
        else
                # Avoid numerical cancellation in the quadrature
                if abs(xu-xs) < 1e1
                        I = (-Inf, xu)
                else
                        # Find a small enough value of the integrand
                        subinterval = collect(LinRange(xs,xu,100000))
                        idx = findfirst(p -> p < 1e-16, [f(x, μ) for x in subinterval])
                        xu = subinterval[idx]
                        I = (-Inf, xu)
                end
        end

        # Define the integral problem over the domain
        if parameters == nothing
                global integral = IntegralProblem(f, I)
        else
                global integral = IntegralProblem(f, I, parameters)
        end

        # Solve the definite integral by using adaptive Gauss-Kronrod quadrature
        quadrature = solve(integral, QuadGKJL(; order=1000); maxiters=100000, reltol=accuracy, abstol=accuracy)

        #display(1.0::Float64/(quadrature.u))
        
        # Return the approximation
        return 1.0::Float64/(quadrature.u)
end

"""
$(TYPEDSIGNATURES)

Fits a polynomial potential of a certain `degree` by solving a linear least-squares problem on its probability density function.

The LLS above is defined by `bins` and `distribution` corresponding to the bins'centerpoints and the heights of an histogram creating an empirical distribution. 
It uses the asymptotic limit of the Fokker-Plank equation of diffusion level `noise` so that an explicit equilibrium (stationary) distribution explicit depends on the potential of the Langevin dynamics. 

If instead of `bins` and `distribution` one only provides a `timeseries`, then the potential, which is assumed to be a cubic, is approximated by solving a nonlinear least-square problem.

## Keyword Arguments
* `N=nothing`: normalisation constant of the equilibrium distribution. If unknown (default) then the function uses the OUP assumption to compute it
* `n_coeff=3::Int`: degree of the polynomial function to be fitted (defaults to `3`)
* `n_bins=nothing`: number of bins to approximate the histogram of the stationary distribution as computed by `fit_distribution` (defaults to the `floor` integer of 2% of the total number of timesteps in `timeseries`)
* `noise=nothing`: additive noise of the scalar SDE generating the `timeseries` (defaults to `std(timeseries)`)
* `initial_guess=nothing`: initial value for the vector of coefficients to be used by the optimiser (if nothing is passed it defaults to the solution of `approximate_solution` which uses the OUP assumption)
* `optimiser=1e-2`: standard deviation of the noise injected in the `initial_guess` (defaults to `1e-2`)
* `attempts=1000::Int`: number of tries on the initial guess for the optimiser in case of failure or numerical instability (defaults to `1000`)
* `verbose=false`: print on screen whether the method converged or not 

## Output
`potential::Tuple`
* `potential.points::Vector{Float64}`:
* `potential.potential::Vector{Float64}`:
* `potential.fit::Vector{Float64}`:

## Example
"""
function fit_potential(bins, distribution, degree, noise::Float64; N = nothing)
        # Compute the diffusion coefficient
        D = (noise^2)/2

        # Filter out the 0-valued entries in the distribution
        idx = findall(x -> x > 0.0, distribution)
        y = [distribution[n] for n in idx]
        x = [bins[n] for n in idx]

        # Define the normalisation constant based on user input
        if N == nothing
                # Assumption of a stationary OUP
                N=1/sqrt(2*pi*D)
        end

        # Compute the distribution on the potential using the stationarity assumption
        U = -D.*log.(y./N)

        # Solve the linear least-square problem to fit U(x) to V
        V = Polynomials.fit(x, U, degree)

        return (
                points = x,
                potential = U,
                fit = V.coeffs
               )
end

function fit_potential(timeseries; n_bins=nothing, noise=nothing, initial_guess=nothing, optimiser=1e-2, attempts=1000, verbose = false)
        # Number of bins for the histogram
        if n_bins == nothing
                n_bins = convert(Int64, floor(0.02*length(timeseries)))
        end
        Nb = n_bins

        # Additive noise in the SDE
        if noise == nothing
                noise = std(timeseries)
        end
        σ = noise

        # Noise in the optimiser's steps
        β = optimiser

        # Fit an empirical distribution to the timeseries data
        bins, hist = fit_distribution(timeseries, n_bins=Nb+1)

        # Initial guess for the non-linear 0-problem 
        if initial_guess == nothing
                initial_guess = (fit_potential(bins, hist, 3, σ)).fit[2:4]
        end

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
        N(μ) = normalise(f, μ)

        # Define the target of the optimisation problem: normalised pdf
        p(x, μ) = N(μ)*exp.(-(1.0::Float64/D).*(c0(μ) .+ μ[1]*x .+ μ[2]*(x.^2) .+ μ[3]*(x.^3)))

        # Define the lower and upper bounds for the coefficients
        lower = [-45.0, -25.0, -40.0] 
        upper = [45.0, 25.0, 40.0]
        
        # First attempt to solve the non-linear least-squares problem
        try
                solution = curve_fit(p, bins, hist, initial_guess, lower=lower, upper=upper, show_trace=true).param
                return (
                        points = bins,
                        guess = initial_guess,
                        fit = solution
                       )
        catch e
                if isa(e, ArgumentError) && (
                                             occursin("matrix contains Infs or NaNs", e.msg) ||
                                             occursin("Initial guess must be within bounds", e.msg)
                                            )
                        # No print
                elseif isa(e, DomainError)
                        # No print
                else
                        rethrow(e)
                end
        end

        # Perturbed attempts
        tries = 0
        while tries < attempts
                try
                        # Perturb the initial guess
                        perturbed_guess = initial_guess + β.*randn(3)
                        # Attempt to solve the nonlinear problem
                        solution = curve_fit(p, bins, hist, perturbed_guess, lower=lower, upper=upper).param
                        return (
                                points = bins,
                                guess = initial_guess,
                                fit = solution
                               )
                catch e
                        if isa(e, ArgumentError) && (
                                                     occursin("matrix contains Infs or NaNs", e.msg) ||
                                                     occursin("Initial guess must be within bounds", e.msg)
                                                    )
                                tries += 1
                        elseif isa(e, DomainError)
                                tries += 1
                        else
                                rethrow(e)
                        end
                end
        end

        # Return the linear solution
        if verbose
                debug("Curve fitting failed after $(tries) attempts.")
        end

        return (
                points = bins,
                guess = initial_guess,
                fit = initial_guess
               )
end
