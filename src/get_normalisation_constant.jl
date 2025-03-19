using Integrals

# Numerically approximate the normalisation constant of a pdf via quadrature
function get_normalisation_constant(f::Function, I::Tuple{Float64, Float64}; accuracy=1e-6)
        # Define the integral problem over the domain
        integral = IntegralProblem(f, I)
        # Solve the definite integral by using adaptive Gauss-Kronrod quadrature
        quadrature = solve(integral, QuadGKJL(); reltol=accuracy)
        # Return the approximation
        return 1.0::Float64/(quadrature.u)
end
