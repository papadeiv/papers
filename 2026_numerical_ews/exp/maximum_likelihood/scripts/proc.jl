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

# Global index for exporting data
global glb_idx_good = 1::Integer
global glb_idx_bad = 1::Integer
 
# Shift the reconstructed potential to match the local minimum of the ground truth
function shift_potential(x0, c)
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
function get_error(Vs, equilibria; Nh=1000)
        # Create uniform partition of the domain of integration
        domain = LinRange(equilibria.unstable[1], equilibria.stable[2], Nh)
        dx = domain[2] - domain[1] 

        # Define the integrand
        E = [(U(x, μ) - Vs(x))^2 for x in domain]

        return sqrt(dx*(sum(E) - 0.5*(E[1]+E[end])))
end
