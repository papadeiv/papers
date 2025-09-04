"""
    Postprocessing script

In here we define the quantities related to the computation of EWSs from raw data.
"""

# Parameters of the scalar potential method
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
 
# Data structures for storing the results of the analysis
c = Matrix{Float64}(undef, Ne, Nc)          # Solutions of the method 
escape = Matrix{Float64}(undef, Ne, 2)      # Estimated escape EWS
parameters = Matrix{Float64}(undef, Ne, 4)  # R.V.s under analysis 
