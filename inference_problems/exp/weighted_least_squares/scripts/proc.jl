"""
    Postprocessing script

In here we define the quantities related to the computation of EWSs from raw data.
"""

# Parameters of the scalar potential method
Nb = convert(Int64, floor(0.010*Nt))            # Number of bins in the histogram
Nc = convert(Int64, 3e0)                        # Solution space dim. of the method 

# Scalar and reconstructed potential 
U(x, μ) =  + μ*x + (1.0/3.0)*x^3
V(x, c) = c[1]*x + c[2]*(x^2) + c[3]*(x^3)

# Ensemble properties 
error = Matrix{Float64}(undef, convert(Int64, Ne), 6)
ensemble_skewness = Vector{Float64}(undef, convert(Int64, Ne))

# Shift the reconstructed potential to match the local minimum of the ground truth
function shift_potential(U::Function, x0, c)
        # Compute the stable equilibrium (center of the shift)
        xs = +(1/(3*c[3]))*(sqrt((c[2])^2 - 3*c[1]*c[3]) - c[2])

        # Compute the shifts
        δx = x0 - xs 
        δy = U(x0, μ) - (Polynomial([0.0; c]))(xs)

        # Define the shifted potential
        Vs(x) = δy + c[1]*(x - δx) + c[2]*(x - δx)^2 + c[3]*(x - δx)^3

        return Vs
end

# Numerical error of the potential reconstruction (trapezoid rule on L2-norm)
function get_error(Vs; Nh=1000)
        # Create uniform partition of the domain of integration
        domain = LinRange(-sqrt(-μ), sqrt(-μ), Nh)
        dx = domain[2] - domain[1] 

        # Define the integrand
        E = [(U(x, μ) - Vs(x))^2 for x in domain]

        return sqrt(dx*(sum(E) - 0.5*(E[1]+E[end])))
end


