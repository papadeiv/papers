"""
Methods to approximate the scalar potential of a conservative vector field from the realisations of a stochastic process whose drift is the gradient of said vector field.

Author: Davide Papapicco
Affil: U. of Auckland
Date: 21-08-2025
"""

#---------------------------------#
#                                 # 
#   potential_reconstruction.jl   #                     
#                                 #
#---------------------------------#

struct Potential
        bins::Vector
        hist::Vector
        coeff::Vector
        guess::Vector
        norm::Float64
end

"""
$(TYPEDSIGNATURES)

Constructs a histogram approximating the probability distribution of a timeseries.

To get meaningful histograms the timeseries `u` should be detrended first.

## Keyword Arguments
* `interval=nothing`: The domain `(a, b)` over which the histogram is constructed (defaults to `(minimum(u),maximum(u))`).
* `n_bins=200::Int64`: Number of bins used for the histogram (defaults to `200`).

## Output
`bins, pdf, normalisation_constant`
* `bins::Vector{Float64}`: center points of the bins of the histogram 
* `pdf::Vector{Float64}`: heights of the bars of the histogram 
* `normalisation_constant::Float64`: `LinearAlgebra.norm(hist)` 

## Example
"""
function fit_distribution(u; interval = nothing, n_bins = 200::Int64)
        # Get the range of values of the distribution
        if interval == nothing
                global u_min = u[argmin(u)]
                global u_max = u[argmax(u)]
                global range = u_max - u_min
        else
                global u_min = interval[1]
                global u_max = interval[end]
                global range = interval[end] - interval[1] 
        end

        # Define the edges of the bins of the histogram
        x = LinRange(u_min - range*0.05, u_max + range*0.05, n_bins)

        # Derive the center-points for the locations of the bins in the plot
        bins = [(x[n+1]+x[n])/2 for n in 1:(length(x)-1)]

        # Fit the histogram through the defined bins
        hist = StatsBase.fit(Histogram, u, x)

        # Normalise the histrogram to get an empirical pdf
        pdf = (LinearAlgebra.normalize(hist, mode = :pdf)).weights
                
        # Export the data
        return bins, pdf, LinearAlgebra.norm(hist)
end

"""
$(TYPEDSIGNATURES)

Uses the stationary asymptotic limit of the OUP to derive the scalar potential given the histogram `bins` and `distribution` values of a probability density and a scalar-values `noise`.

## Keyword Arguments
* `N=nothing`: the normalisation constant of the pdf (defaults to `1/sqrt(pi*noise^2)`). 

## Output
`xs, Vs`
* `xs::Vector{Float64}`: discrete coordinates of the inverted potential 
* `Vs::Vector{Float64}`: values of the inverted potential 

## Example
"""
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
                N=1/sqrt(2*pi*D)
        end

        # Compute the distribution on the potential using the stationarity assumption
        Vs = -D.*log.(ys./N)

        # Return the filtered datapoints and their location
        return xs, Vs
end

"""
$(TYPEDSIGNATURES)

Solves the linear least-squares problem defined by `data_x` and `data_y` to fit a polynomial scalar potential.

## Keyword Arguments
* `degree=3::Int`: degree of the fitted polynomial

## Output
* `interpolant::Function`: result of the least-square fit

## Example
"""
function approximate_potential(data_x, data_f; degree=3::Int64)
        interpolant = Polynomials.fit(data_x, data_f, degree)
        return interpolant
        #=
        # Equivalent self implementation of polynomial least-squares
        # Get number of observations in the data
        N = length(data_x)
        # Initialise and construct the polynomial model matrix
        A = Matrix{Float64}(undef, N, degree+1)
        for n in 1:N
                for d in 0:degree
                        A[n,d+1] = data_x[n]^d 
                end
        end
        # Compute the covariance matrix
        Σ = inv(transpose(A)*A)
        # Compute the pseudoinverse of A
        M = Σ*transpose(A)
        # Initialise and construct the observation vector
        b = data_f
        # Solve the normal equations
        c = M*b
        # Compute the residual
        r = b - A*c
        # Define the polynomial
        U = Polynomial(c)
        # Export the polynomial fit and the residual
        return A, Σ, c, r
        =#
end

"""
$(TYPEDSIGNATURES)

Computes the approximate normalisation constant of a pdf by Gauss-Kronrod quadrature.

Here `f::Function` is assumed to depend on a cubic potential parametrised by `parameters` and defined over a finite domain.

## Keyword Arguments
* `accuracy=1e-8`: relative accuracy criterion for the quadrature iterations (defaults to `1e-8`).

## Output
* `normalisation_constant::Float64`: normalisation constant to make `f` a pdf 

## Example
"""
function get_normalisation_constant(f::Function, parameters; accuracy=1e-8)
        # Compute the location of the local minima and maxima
        μ = parameters
        xs = +(1/(3*μ[3]))*(sqrt((μ[2])^2 - 3*μ[1]*μ[3]) - μ[2])
        xu = -(1/(3*μ[3]))*(sqrt((μ[2])^2 - 3*μ[1]*μ[3]) + μ[2])

        # Define the integration interval
        I = (-Inf,Inf)
        if xs > xu
                I = (xu, +Inf) 
        else
                I = (-Inf, xu) 
        end

        # Define the integral problem over the domain
        if parameters == nothing
                global integral = IntegralProblem(f, I)
        else
                global integral = IntegralProblem(f, I, parameters)
        end
        # Solve the definite integral by using adaptive Gauss-Kronrod quadrature
        quadrature = solve(integral, QuadGKJL(); maxiters=10000, reltol=accuracy, abstol=accuracy)
        # Return the approximation
        return 1.0::Float64/(quadrature.u)
end

"""
$(TYPEDSIGNATURES)

Sets up and solves the non-linear least-square optimisation problem to fit a scalar potential to a `timeseries` using the stationarity assumption.

The potential is assumbed to be cubic.

## Keyword Arguments
* `n_coeff=3::Int`: degree of the polynomial function to be fitted (defaults to `3`)
* `n_bins=nothing`: number of bins to approximate the histogram of the stationary distribution as computed by `fit_distribution` (defaults to the `floor` integer of 2% of the total number of timesteps in `timeseries`)
* `noise=nothing`: additive noise of the scalar SDE generating the `timeseries` (defaults to `std(timeseries)`)
* `initial_guess=nothing`: initial value for the vector of coefficients to be used by the optimiser (if nothing is passed it defaults to the solution of `approximate_solution` which uses the OUP assumption)
* `optimiser=1e-2`: standard deviation of the noise injected in the `initial_guess` (defaults to `1e-2`)
* `attempts=1000::Int`: number of tries on the initial guess for the optimiser in case of failure or numerical instability (defaults to `1000`)
* `verbose=false`: print on screen whether the method converged or not 

## Output
* `coefficients::Vector{Float64}`: 3 polynomial coefficients of the cubic solving the nonlinear problem (c0 is assumed to be 0 and does not impact the potential's shape)

## Example
"""
function fit_potential(timeseries; n_coeff=3, n_bins=nothing, noise=nothing, initial_guess=nothing, optimiser=1e-2, attempts=1000, verbose = false)
        # Number of coefficients for the non-linear least-squares problem
        Nc = n_coeff

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

        # Initial guess for non-linear 0-problem 
        if initial_guess == nothing
                xs, Vs = invert_equilibrium_distribution(bins, hist, σ)
                initial_guess = (approximate_potential(xs, Vs, degree=Nc).coeffs[2:(Nc+1)])
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
        N(μ) = get_normalisation_constant(f, μ)

        # Define the target of the optimisation problem: normalised pdf
        p(x, μ) = N(μ)*exp.(-(1.0::Float64/D).*(c0(μ) .+ μ[1]*x .+ μ[2]*(x.^2) .+ μ[3]*(x.^3)))

        # Define the lower and upper bounds for the coefficients
        lower = [-45.0, -25.0, -40.0] 
        upper = [45.0, 25.0, 40.0]

        
        # First attempt to solve the non-linear least-squares problem
        try
                solution = curve_fit(p, bins, hist, initial_guess, lower=lower, upper=upper).param
                return solution
        catch e
                if isa(e, ArgumentError) && (
                                             occursin("matrix contains Infs or NaNs", e.msg) ||
                                             occursin("Initial guess must be within bounds", e.msg)
                                            )
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
                        perturbed_guess = initial_guess + β.*randn(Nc)
                        # Attempt to solve the nonlinear problem
                        solution = curve_fit(p, bins, hist, perturbed_guess, lower=lower, upper=upper).param
                        return solution
                catch e
                        if isa(e, ArgumentError) && (
                                                     occursin("matrix contains Infs or NaNs", e.msg) ||
                                                     occursin("Initial guess must be within bounds", e.msg)
                                                    )
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
        return initial_guess
end

"""
$(TYPEDSIGNATURES)

Compute the vertical and horizontal shifts of a polynomial scalar potential generated by `c` by centering it on the local minima `x0` of a target potential `U` (usually the ground truth) evaluated at parameter value `μ`.

## Output
* `xs, Vs`: locations `xs` of the shifted potential `Vs` 
"""
function shift_potential(U::Function, x0, μ, c)
        # Compute the stable equilibrium (center of the shift)
        xs = +(1/(3*c[3]))*(sqrt((c[2])^2 - 3*c[1]*c[3]) - c[2])

        # Compute the shifts
        δx = x0 - xs 
        δy = U(x0, μ) - (Polynomial([0.0; c]))(xs)

        # Define the shifted potential
        Vs(x) = δy + c[1]*(x - δx) + c[2]*(x - δx)^2 + c[3]*(x - δx)^3

        return xs, Vs
end
