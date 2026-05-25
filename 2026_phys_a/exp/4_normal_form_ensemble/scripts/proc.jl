"""
    Postprocessing script

Collection of quantities and functions used to postprocess and analyse the results of a simulation.
"""

# Scalar potential of the conservative system 
U(x, μ) =  + μ*x + (1.0/3.0)*x^3                # Potential (ground truth)
V(x, c) = c[1]*x + c[2]*(x^2) + c[3]*(x^3)      # Potential

# SOlve the LLS problem
function solve_lls(solution)
        # Define the observation vectors 
        Xn = solution[1:end-1]
        Y  = (solution[2:end] .- solution[1:end-1])./dt

        # Assemble the model matrix
        A = hcat(ones(length(Xn)), Xn, Xn.^2)

        # Solve the (linear) least-squares problem
        β = A\Y

        # Compute the coefficients of the potential
        θ = [-β[1], -β[2]/2, -β[3]/3]
        return θ 
end
 
# Shift the reconstructed potential to match the local minimum of the ground truth
function shift_potential(U::Function, x0, c)
        # Compute the stable equilibrium (center of the shift)
        xs = +(1/(3*c[3]))*(sqrt((c[2])^2 - 3*c[1]*c[3]) - c[2])

        # Compute the shifts
        δx = x0 - xs 
        δy = U(x0, μ) - (Polynomial([0.0; c]))(xs)

        # Define the shifted potential
        Vs(x) = δy + c[1]*(x - δx) + c[2]*(x - δx)^2 + c[3]*(x - δx)^3

        return xs, Vs
end

#=
# Numerical error of the potential reconstruction (trapezoid rule on L2-norm)
function get_error(Vs; Nh=1000)
        # Create uniform partition of the domain of integration
        domain = LinRange(equilibria.unstable[1], equilibria.stable[1], Nh)
        dx = domain[2] - domain[1] 

        # Define the integrand
        E = [(U(x, μ) - Vs(x))^2 for x in domain]

        return sqrt(dx*(sum(E) - 0.5*(E[1]+E[end])))
end

function analyse(solutions)
        # Reconstruct a shifted potential to match the ground truth (for error and plotting purposes)
        xs, Vs = shift_potential(U, sqrt(-μ), solutions)

        # Compute estimated stable and unstable equilibria of the cubic
        xs = +(1/(3*solutions[3]))*(sqrt((solutions[2])^2 - 3*solutions[1]*solutions[3]) - solutions[2])
        xu = -(1/(3*solutions[3]))*(sqrt((solutions[2])^2 - 3*solutions[1]*solutions[3]) + solutions[2])

        # Compute the escape rate's components
        ΔV = abs(V(xu, solutions) - V(xs, solutions))
        prefactor = sqrt(abs(Vxx(xu, solutions))*Vxx(xs, solutions))
        decay = exp(-ΔV)

        # Define empty vector
        analysis = Float64[]

        # Compute the random variables at the estimated stable equilibrium 
        push!(analysis, xs)                     # Equilibrium
        push!(analysis, V(xs, solutions))       # Potential value
        push!(analysis, Vxx(xs, solutions))     # Curvature value
        push!(analysis, exp(V(xs,solutions)))   # Large deviation

        # Compute the random variables at the estimated unstable equilibrium 
        push!(analysis, xu)                     # Equilibrium
        push!(analysis, V(xu, solutions))       # Potential value
        push!(analysis, Vxx(xu, solutions))     # Curvature value
        push!(analysis, exp(V(xu,solutions)))   # Large deviation

        # Compute the random variables associated to the escape ews
        push!(analysis, ΔV)                     # Potential barrier 
        push!(analysis, prefactor)              # Prefactor of the escape rate 
        push!(analysis, decay)                  # Exponential decay of the escape rate 
        push!(analysis, get_error(Vs))          # Numerical error of the shifted potential reconstruction

        # Update the results vector
        push!(results, analysis)

        # Return the shifted potential (for plotting)
        return Vs 
end
=#
