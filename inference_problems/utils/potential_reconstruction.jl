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
