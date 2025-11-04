"""
    Postprocessing script

In here we define the quantities related to the computation of EWSs from raw data.
"""

# Parameters of the scalar potential method
Nb = convert(Int64, floor(0.010*Nt))            # Number of bins in the histogram
Nc = convert(Int64, 3e0)                        # Solution space dim. of the method 
β = 1e-3                                        # Std of the guess perturbation 

# Scalar and reconstructed potential 
U(x, μ) =  + μ*x + (1.0/3.0)*x^3
V(x, c) = c[1]*x + c[2]*(x^2) + c[3]*(x^3)

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
