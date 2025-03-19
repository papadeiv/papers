using LinearAlgebra, Polynomials

# Approximate the potential function from the data using polynomial least-squares regression
function fit_potential(x, y; degree=3::Int64)
        # Perform the least-squares regression
        A, Σ, c, r = get_coefficients(x, y, degree=degree)
        # Construct the Polynomial structure out of the coefficients
        U = Polynomial(c) 
        # Export the polynomial fit and its residual
        return U, r
end
