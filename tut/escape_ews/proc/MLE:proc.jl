using QuadGK, Optim, Statistics

function estimate_mle(X::Vector{Float64})

    # Define the negative log-likelihood
    function negloglik(params)
        c1, c2, c3 = params

        sigma = 0.250

        xs = (1/(3*c3))*(sqrt((c2)^2 - 3*c1*c3) - c2)
        c0 = - c1*xs - c2*(xs)^2 - c3*(xs)^3

        V(x) = c0 + c1*x + c2*x^2 + c3*x^3

        # Compute unnormalized log-density at each data point
        V_vals = V.(X)
        nll_data = (2 / sigma^2) * sum(V_vals)

        # Choose integration bounds based on data range
        xu = -(1/(3*c3))*(sqrt((c2)^2 - 3*c1*c3) + c2)
        xmin, xmax = xu, +10.0 

        # Estimate normalizing constant Z(theta)
        integrand(x) = exp(-2 * V(x) / sigma^2)
        Z, err = quadgk(integrand, xmin, xmax, rtol=1e-6)

        if Z <= 0 || !isfinite(Z)
            return Inf
        end

        println("params = ", params, ", N = ", Z)

        # Compute the negative log-likelihood
        nll_norm = length(X) * log(Z)
        nll = nll_data + nll_norm 

        return nll
    end

    # Initial guess
    init = [-5e-6, 1.33, 1.08]

    result = optimize(negloglik, init, NelderMead(); autodiff = :forward)

    return Optim.minimizer(result)
end

include("../../../inc/IO.jl")

# Import the data from csv 
solution = readin("../data/MLE/solution.csv")
t = solution[:,1]
u = solution[:,2]

# Find the MLE
ϑ = estimate_mle(u) 
display(ϑ)
