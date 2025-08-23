# Numerically approximate the normalisation constant of a pdf via quadrature

function get_normalisation_constant(f::Function, I::Tuple{Float64, Float64}; parameters=nothing, accuracy=1e-8)
        # Define the integral problem over the domain and with the user's parameters 
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
