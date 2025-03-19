using LinearAlgebra

function estimate_error(x, y, σ, d; estimate=1)
        # Number of samples
        n = length(x)
        # Number of (polynomial coefficients)
        m = d+1
        # Get the linear problem of the least-squares regression
        A, Σ = get_coefficients(x, y, degree=d)[1:2]
        # Assemble the pseudoinverse of A
        H = Σ*transpose(A)
        # Compute the upper bound according to the user choice
        if estimate==1
                global e = [3*σ*opnorm(H, Inf) for k in 1:m]
        elseif estimate==2
                global e = [3*σ*n*maximum(H) for k in 1:m]
        elseif estimate==3
                global e = [3*σ*sqrt(n)*norm(H[k,:]) for k in 1:m]
        elseif estimate==4
                global e = [3*σ*sqrt(n*m)/minimum(svd(A).S) for k in 1:m]
        elseif estimate==5
                M = (1/n).*(transpose(A)*A)
                global e = [σ*sqrt(2)/sqrt(pi*n*eigmin(M)) for k in 1:m]
        end
        # Return the error estimate
        return e
end
