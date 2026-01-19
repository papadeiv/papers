"""
    Postprocessing script

Collection of quantities and functions used to postprocess and analyse the results of a simulation.
"""

# Parameters of the scalar potential method
Nb = convert(Int64, floor(0.010*Nt))            # Number of bins in the histogram
Nc = convert(Int64, 3e0)                        # Solution space dim. of the method 
Na = convert(Int64, 1e4)                        # Number of attempts per guess 
Nx = 0::Integer                                 # Number of escapes prior to the end of the simulation
β = 1e-3                                        # Std of the guess perturbation 

# Scalar potential of the conservative system 
U(x, μ) =  + μ*x + (1.0/3.0)*x^3                # Potential (ground truth)
Uxx(x, μ) = + 2*x                               # Second derivative ( == Jacobian) 

# Reconstructed dynamics 
V(x, c) = c[1]*x + c[2]*(x^2) + c[3]*(x^3)      # Potential
Vxx(x, c) = 2*c[2] + 6*c[3]*x                   # Second derivative
 
# Data structures
solutions = Vector{Vector{Float64}}()           # Solutions of the inference method 
results = Vector{Vector{Float64}}()             # Statistical distribution over the ensemble

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

        # Update the results vector
        push!(results, analysis)

        # Return the shifted potential (for plotting)
        return Vs 
end

function characterise()
        # Define type of runs
        runs = ["bad", "good"]

        # Loop over the type
        for n in 1:2 
                # Define empty matrix
                characterisation = Matrix{Float64}(undef, 50, 6)

                # Loop over the runs
                for m in 1:50
                        # Import the run
                        run = Matrix(CSV.read("../../res/data/"*runs[n]*"/$m.csv", DataFrame))
                        t = run[:,1] 
                        u = run[:,2] 

                        # Compute its statistical characteristics
                        characterisation[m,1] = moment(u, 2, mean(u))                 # Centered variance
                        characterisation[m,2] = moment(u, 3, mean(u))/std(u)^3        # Centered skewness 
                        characterisation[m,3] = moment(u, 4, mean(u))/std(u)^4 - 3.0  # Centered (excess) kurtosis
                        characterisation[m,4] = moment(u, 5, mean(u))/std(u)^5
                        characterisation[m,5] = moment(u, 6, mean(u))/std(u)^6 - 15.0
                        characterisation[m,6] = moment(u, 7, mean(u))/std(u)^7
                end

                # Construct the dataframe with its header
                df = DataFrame(variance = characterisation[:,1], skewness = characterisation[:,2], kurtosis = characterisation[:,3])

                # Display disjoint union of only some columns
                println(runs[n]*" runs")
                @printf("2nd moment = %.6f | 3rd moment = %.6f | 4th moment = %.6f | 5th moment = %.6f | 6th moment = %.6f | 7th moment = %.6f\n\n", 
                        mean(characterisation[:,1]),mean(characterisation[:,2]),mean(characterisation[:,3]),mean(characterisation[:,4]),mean(characterisation[:,5]),mean(characterisation[:,6]))
        end
end
