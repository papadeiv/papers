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

# Jacobian approximation (central finite differences)
function approx_jacobian(f, x, c; h=eps()^(1.0/3.0))
        n = length(x)
        m = length(c)
        J = zeros(Float64, n, m)
        for j in 1:m
                c_plus = copy(c)
                c_minus = copy(c)
                c_plus[j] += h 
                c_minus[j] -= h
                f_plus = f.(x, Ref(c_plus)) 
                f_minus = f.(x, Ref(c_minus)) 
                J[:,j] .= (f_plus .- f_minus)./(2*h)
        end
        return J
end

function cost(x, y, c, model)
        r = model.(x, Ref(c)) .- y
        0.5 * sum(r.^2)
end

# Nonlinear regression
function fit_gradient_descent(xdata, ydata, model, guess; jac=nothing, maxiter=1000, α=1e-3, tol=1e-10, verbose=true)
        c = copy(guess)
        for k in 1:maxiter
                xs = +(1/(3*c[3]))*(sqrt((c[2])^2 - 3*c[1]*c[3]) - c[2])
                xu = -(1/(3*c[3]))*(sqrt((c[2])^2 - 3*c[1]*c[3]) + c[2])
                @printf("iter %5d | c = [%.8f, %.8f, %.8f] | xs = %.8f | xu = %.8f\n", k-1, c[1], c[2], c[3], xs, xu)

                J = zeros(Float64, length(xdata), length(c))
                if jac==nothing
                        J = approx_jacobian(model, xdata, c)
                else
                        J = hcat([jac(x, c)[1] for x in xdata], [jac(x, c)[2] for x in xdata]) 
                end
                r = model.(xdata, Ref(c)) .- ydata    # Residuals
                g = J'*r                              # Gradient cost

                #@printf("iter %5d | ||r||_2 = %.8f | ||J||_1 = %.8f\n", k-1, norm(r, 2), norm(J, 1))
                #df = DataFrame(dfdc1 = J[:,1], dfdc2 = J[:,2], dfdc3 = J[:,3])
                #CSV.write("../../res/data/jacobian/$k.csv", df)

                @printf("           | g = [%.8f, %.8f, %.8f]\n", g[1], g[2], g[3])
                @printf("\n")

                c_new = c .- α .* g
                if norm(c_new - c) < tol
                        verbose && @printf("Converged at iter %d\n", k)
                        return c_new
                end

                c = c_new
                if verbose && k % 10 == 0
                        # Print iteration info
                        #@printf("iter %5d | cost = %.6e | gradient = [%.6f, %.6f, %.6f] | c = [%.6f, %.6f, %.6f]\n", k, cost(xdata, ydata, c, model), g[1], g[2], g[3], c[1], c[2], c[3])

                        # Plot the Jacobian approximation 
                        fig = Figure(; size = (1500, 500))
                        ax1 = Axis(fig[1,1], xlabel=L"x", ylabel=L"\partial_{c_1}\;f")
                        lines!(ax1, xdata, J[:,1], color=:orange)
                        ax2 = Axis(fig[1,2], xlabel=L"x", ylabel=L"\partial_{c_2}\;f")
                        lines!(ax2, xdata, J[:,2], color=:orange)
                        ax3 = Axis(fig[1,3], xlabel=L"x", ylabel=L"\partial_{c_3}\;f")
                        lines!(ax3, xdata, J[:,3], color=:orange)
                        savefig("cubic/jacobian/bad/$k.png", fig)
                end
        end
        verbose && println("Reached max iterations")
        return c
end
