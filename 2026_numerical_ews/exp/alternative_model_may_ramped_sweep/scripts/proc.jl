"""
    Postprocessing script

Collection of quantities and functions used to postprocess and analyse the results of a simulation.
"""

# Scalar potential of the conservative system 
U(x, μ) = -(r/2.0)*x^2 + (r/(3.0*k))*x^3 + μ*(x - sqrt(h)*atan(x/sqrt(h)))    # Ground truth
V(x, c) = c[1]*x + c[2]*(x^2) + c[3]*(x^3)                                    # Arbitrary cubic potential

# Data structure to store the approximate early-warning signals
ews = [Matrix{Float64}(undef, convert(Integer, Ne), length(M0)) for _ in 1:3]

# Build the critical manifold (branch of bifurcation diagram) in the given parameter range 
function compute_crit_man(initial_condition, parameter_range)
        # Redefine the vector field for the API 
        F(x, μ) = @. f(x, μ[1])

        # Define the 0-problem
        zero_problem = BifurcationProblem(F, initial_condition, parameter_range[1], 1; record_from_solution = (x, p; k...) -> x = x[1])

        # Compute the branches of the diagram
        bifurcation_diagram = continuation(zero_problem, PALC(), ContinuationPar(p_min = parameter_range[1], p_max = parameter_range[end], dsmax = 0.01, max_steps=1000))

        # Store the bifurcation diagram in the appropriate data structure
        parameter_value = Vector{Float64}(undef, length(bifurcation_diagram.sol))
        equilibrium_drift = Vector{Float64}(undef, length(bifurcation_diagram.sol))
        for (index, point) in enumerate(bifurcation_diagram.sol)
                parameter_value[index] = point.p
                equilibrium_drift[index] = (point.x)[1]
        end

        # Interpolate the values across the bifurcation diagram to match the sampling points in the parameter range
        interpolation = linear_interpolation(parameter_value, equilibrium_drift)
        critical_manifold = interpolation.(parameter_range)

        return critical_manifold
end

# Solve the MLE on the EM approximation
function solve_EM_MLE(solution)
        Xn = solution[1:end-1]
        Y  = (solution[2:end] .- solution[1:end-1])./δt
        Φ = hcat(ones(length(Xn)), Xn, Xn.^2)
        β = Φ\Y
        c = [-β[1], -β[2]/2, -β[3]/3]
        return c
end

# Compute a parameter shift for the reconstructed potential
function shift_potential(x0, c, μ)
        # Compute the stable equilibrium (center of the shift)
        xs = +(1/(3*c[3]))*(sqrt((c[2])^2 - 3*c[1]*c[3]) - c[2])

        # Compute the shifts
        δx = x0 - xs 
        δy = U(x0, μ) - (Polynomials.Polynomial([0.0; c]))(xs)

        # Construct the shifted potential
        Vs(x) = δy + c[1]*(x - δx) + c[2]*(x - δx)^2 + c[3]*(x - δx)^3

        return Vs
end

# Compute the early-warning signal
function compute_ews(solutions)
        # Compute estimated stable and unstable equilibria of the cubic
        xs = +(1/(3*solutions[3]))*(sqrt((solutions[2])^2 - 3*solutions[1]*solutions[3]) - solutions[2])
        xu = -(1/(3*solutions[3]))*(sqrt((solutions[2])^2 - 3*solutions[1]*solutions[3]) + solutions[2])

        # Compute the ews
        ΔV = abs(V(xu, solutions) - V(xs, solutions))
        return exp(-ΔV)
end
