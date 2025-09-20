"""
    Postprocessing script

In here we define the quantities related to the computation of EWSs from raw data.
"""

# Parameters of the scalar potential method
Nμ = length(μ0)                               # Number of parameter ranges 
Nb = max(convert(Int64, floor(0.01*Ne)), 10)  # Number of bins in the plotting histogram 
Nc = convert(Int64, 3e0)                      # Solution space dim. of the method 
Na = convert(Int64, 1e4)                      # Number of attempts per guess 
β = 1e-2                                      # Std of the guess perturbation 

# Scalar potential of the conservative system 
U(x, μ) = μ*x + x^2 - x^3 + (1/5)*(x^4)       # Potential (ground truth)
Uxx(x, μ) = 2 - 6*x + (12/5)*(x^2)            # Second derivative ( == Jacobian) 

# Reconstructed dynamics 
V(x,c) = c[1]*x + c[2]*(x^2) + c[3]*(x^3)     # Potential
Vxx(x,c) = 2*c[2] + 6*c[3]*x                  # Second derivative
 
# Data structures for storing the results of the analysis
coefficients = Matrix{Float64}(undef, Ne, Nc) # Solutions of the method 
escape = Matrix{Float64}(undef, Ne, Nμ)       # Estimated escape EWS
analysis = Vector{Vector{Float64}}()          # Transformations of the solutions 

# Perform statistical analysis on the NLLS solutions 
function analyse(solutions)
        # Estimated equilibria
        xs = +(1/(3*solutions[3]))*(sqrt((solutions[2])^2 - 3*solutions[1]*solutions[3]) - solutions[2])
        xu = -(1/(3*solutions[3]))*(sqrt((solutions[2])^2 - 3*solutions[1]*solutions[3]) + solutions[2])
 
        # Value of the potential at equilibria 
        Vs = V(xs, solutions)     
        Vu = V(xu, solutions)     

        # Local curvature of the potential at equilibria
        Vxxs = Vxx(xs, solutions)
        Vxxu = Vxx(xu, solutions)

        # Energy barrier
        ΔV = Vu - Vs

        # Large-deviation principles
        LDPs = exp(Vs)
        LDPu = exp(Vu)
        LDP = exp(-ΔV/D)

        return [xs, xu, Vs, Vu, Vxxs, Vxxu, ΔV, LDPs, LDPu, LDP]
end
