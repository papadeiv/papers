"""
    Postprocessing script

In here we define the quantities related to the computation of EWSs from raw data.
"""

# Parameters of the scalar potential method
width = 0.250::Float64                      # Relative size of the sliding window
Nc = convert(Int64, 3e0)                    # Solution space dim. of the method 
Na =  convert(Int64, 1e4)                   # Number of attempts per guess 
β = 1e-2                                    # Std of the guess perturbation 

# Scalar potential of the conservative system 
U(x, μ) = μ*x + x^2 - x^3 + (1/5)*(x^4)     # Potential (ground truth)
Ux(x, μ) = -f(x, μ)                         # First derivative ( == vector field) 
Uxx(x, μ) = 2 - 6*x + (12/5)*(x^2)          # Second derivative ( == Jacobian) 
U3x(x, μ) = - 6 + (24/5)*x                  # Third derivative ( == Hessian) 

# Reconstructed dynamics 
V(x,c) = c[1]*x + c[2]*(x^2) + c[3]*(x^3)   # Potential
Vxx(x,c) = 2*c[2] + 6*c[3]*x                # Second derivative

# Stationary probability distribution
p(x,c) = exp(-(1.0::Float64/D)*(V(x,c)))

# Data structures for storing the results of the analysis
solutions = Vector{Vector{Float64}}()        # Solutions of the method 
escapes = Vector{Vector{Float64}}()          # Estimated escape EWS
parameters = Vector{Vector{Float64}}()       # R.V.s under analysis 
error = Vector{Float64}()                    # Numerical error of the reconstrucion

# Converts the non-stationary timeseries into an ensemble of subseries associated to the strides of a sliding window 
function preprocess_solution(timestamps, timeseries, width)
        # Find the tipping point
        tipping = find_tipping(timeseries)
        idx = tipping.index

        # Extract the subseries up to the tipping
        t = timestamps[1:idx] 
        u = timeseries[1:idx]
        Nt = length(u)

        # Get the sliding window parameters
        window = get_window_parameters(Nt, width)
        Nw = window.size 
        Ns = window.strides

        # Convert the sliding window subseries into ensemble timeseries
        printstyled("Converting the truncated sample path to an ensemble of ", Ns," trajectories of ", Nw, " steps\n"; bold=true, underline=true, color=:light_blue)
        timesteps = [t[n:(n+Nw-1)] for n in 1:Ns] 
        ensemble = [u[n:(n+Nw-1)] for n in 1:Ns] 

        # Export the parameters of the ensemble problem 
        return (
                tipping_point = idx,
                timesteps = timesteps,
                trajectories = ensemble 
               ) 
end

# Location of the local minima and local maxima of the ground truth and reconstructed potential 
function get_bounds(μ, coeffs)
        # Compute the extrema of the ground truth 
        true_equilibria = get_equilibria(f, μ, domain=[-10,10])

        # Compute the extrema of the (shifted) reconstruction
        xs = +(1/(3*coeffs[3]))*(sqrt((coeffs[2])^2 - 3*coeffs[1]*coeffs[3]) - coeffs[2])
        xu = -(1/(3*coeffs[3]))*(sqrt((coeffs[2])^2 - 3*coeffs[1]*coeffs[3]) + coeffs[2])

        return (
                true_min = true_equilibria.stable[2],
                true_max = true_equilibria.unstable[1],
                approx_min = xs,
                approx_max = xu 
               ) 
end

# Numerical error of the potential reconstruction (trapezoid rule on L2-norm)
function get_error(μ, Vs; Nh=1000)
        # Compute the equilibria of the vector field (extrema of the potential)
        equilibria = get_equilibria(f, μ, domain=[-10,10])

        # Create uniform partition of the domain of integration
        domain = LinRange(equilibria.unstable[1], equilibria.stable[2], Nh)
        dx = domain[2] - domain[1] 

        # Define the integrand
        E = [(U(x, μ) - Vs(x))^2 for x in domain]

        return sqrt(dx*(sum(E) - 0.5*(E[1]+E[end])))
end
