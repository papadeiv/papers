# Derive the stationary points of a polynomial 

function get_stationary_points(V::Polynomial; domain=[-Inf,Inf])
        # Differentiate the polynomial
        dV = derivative(V)

        # Find the roots of the first-derivative of the polynomial
        points = roots(dV)

        # Return the points sorted from smallest to largest 
        return sort(points) 
end
