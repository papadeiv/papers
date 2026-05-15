"""
    Postprocessing script

Collection of quantities and functions used to postprocess and analyse the results of a simulation.
"""

# Scalar potential of the conservative system 
U(x, μ) = 1 + μ*x - 2*(x^2) + x^4               # Ground truth...
V(x, c) = c[1]*x + c[2]*(x^2) + c[3]*(x^3)      # Arbitrary
ρ(x, c) = exp(-V(x, c)/D)
xs(C) = +(1/(3*C[3]))*(sqrt((C[2])^2 - 3*C[1]*C[3]) - C[2])
xu(C) = -(1/(3*C[3]))*(sqrt((C[2])^2 - 3*C[1]*C[3]) + C[2])

# Error array 
error = Matrix{Float64}(undef, convert(Integer, Ne), 4)

# Parameter index
global μ_idx = 1::Integer

# Shift the reconstructed potential to match the local minimum of the ground truth
function shift_potential(x0, c, μ)
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
function get_error(Vs, equilibria, μ; Nh=5000)
        # Create uniform partition of the domain of integration
        domain = LinRange(equilibria.unstable[1], equilibria.stable[2], Nh)
        dx = domain[2] - domain[1] 

        # Define the integrands
        E_num = [(U(x, μ) - Vs(x))^2 for x in domain]
        E_den = [(U(x, μ))^2 for x in domain]

        # Compute the ratio
        numerator = sqrt(dx*(sum(E_num) - 0.5*(E_num[1]+E_num[end])))
        denominator = sqrt(dx*(sum(E_den) - 0.5*(E_den[1]+E_den[end])))
        relative_error = numerator/denominator

        # Return the relative error in L-2
        return relative_error 
end
