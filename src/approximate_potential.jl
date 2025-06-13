# Approximate the potential function from the data using polynomial least-squares regression 

function approximate_potential(data_x, data_f; degree=3::Int64)
        interpolant = Polynomials.fit(data_x, data_f, degree)
        return interpolant
        #=
        # Equivalent self implementation of polynomial least-squares
        # Get number of observations in the data
        N = length(data_x)
        # Initialise and construct the polynomial model matrix
        A = Matrix{Float64}(undef, N, degree+1)
        for n in 1:N
                for d in 0:degree
                        A[n,d+1] = data_x[n]^d 
                end
        end
        # Compute the covariance matrix
        Σ = inv(transpose(A)*A)
        # Compute the pseudoinverse of A
        M = Σ*transpose(A)
        # Initialise and construct the observation vector
        b = data_f
        # Solve the normal equations
        c = M*b
        # Compute the residual
        r = b - A*c
        # Define the polynomial
        U = Polynomial(c)
        # Export the polynomial fit and the residual
        return A, Σ, c, r
        =#
end
