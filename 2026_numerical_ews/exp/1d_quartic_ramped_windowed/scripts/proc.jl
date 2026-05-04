"""
    Postprocessing script

Collection of quantities and functions used to postprocess and analyse the results of a simulation.
"""

# Parameters of the scalar potential method
width = 0.333::Float64                      # Relative size of the sliding window
Na =  convert(Int64, 1e4)                   # Number of attempts per guess 
β = 1e-2                                    # Std of the guess perturbation 

# Scalar potential of the conservative system 
U(x, μ) = 1 + μ*x - 2*(x^2) + x^4               # Ground truth
V(x, c) = c[1]*x + c[2]*(x^2) + c[3]*(x^3)      # Arbitrary cubic potential

# Data structures
solutions = Vector{Vector{Float64}}()           # Solutions of the inference method 
results = Vector{Vector{Float64}}()             # Statistical distribution over the ensemble

# Converts the non-stationary timeseries into an ensemble of subseries associated to the strides of a sliding window 
function preprocess_solution(timestamps, timeseries, width)
        # Find the tipping point
        tipping = find_tipping(timeseries)
        idx = tipping.index

        # Extract the subseries up to the tipping
        t = timestamps[1:idx] 
        u = timeseries[1:idx]
        Nt = length(u)

        # Assemble the sliding window
        window = build_window(Nt, width)
        Nw = window.size 
        Ns = window.strides

        # Convert the sliding window subseries into an ensemble of timeseries
        printstyled("Converting the truncated sample path to an ensemble of ", Ns," trajectories of ", Nw, " steps\n"; bold=true, underline=true, color=:light_blue)
        timesteps = [@view t[n:(n+Nw-1)] for n in 1:Ns] 
        ensemble = [@view u[n:(n+Nw-1)] for n in 1:Ns] 

        # Export the parameters of the ensemble problem 
        return (
                tipping_point = idx,
                timesteps = timesteps,
                trajectories = ensemble 
               ) 
end

function shift_potential(U::Function, x0, μ, c)
        # Compute the stable equilibrium (center of the shift)
        xs = +(1/(3*c[3]))*(sqrt((c[2])^2 - 3*c[1]*c[3]) - c[2])

        # Compute the shifts
        δx = x0 - xs 
        δy = U(x0, μ) - (Polynomial([0.0; c]))(xs)

        # Define the shifted potential
        Vs(x) = δy + c[1]*(x - δx) + c[2]*(x - δx)^2 + c[3]*(x - δx)^3

        return xs, Vs
end

# Numerical error of the potential reconstruction (trapezoid rule on L2-norm)
function get_error(Vs, μ; Nh=1000)
        # Create uniform partition of the domain of integration
        domain = LinRange(equilibria.unstable[1], equilibria.stable[2], Nh)
        dx = domain[2] - domain[1] 

        # Define the integrand
        E = [(U(x, μ) - Vs(x))^2 for x in domain]

        return sqrt(dx*(sum(E) - 0.5*(E[1]+E[end])))
end

function analyse(solutions, parameter)
        # Reconstruct a shifted potential to match the ground truth (for error and plotting purposes)
        xs, Vs = shift_potential(U, x0[1], parameter, solutions)

        # Compute estimated stable and unstable equilibria of the cubic
        xs = +(1/(3*solutions[3]))*(sqrt((solutions[2])^2 - 3*solutions[1]*solutions[3]) - solutions[2])
        xu = -(1/(3*solutions[3]))*(sqrt((solutions[2])^2 - 3*solutions[1]*solutions[3]) + solutions[2])

        # Define empty vector
        analysis = Float64[]
        push!(analysis, parameter)

        # Compute the ews
        ΔV = abs(V(xu, solutions) - V(xs, solutions))
        push!(analysis, exp(-ΔV))

        # Compute the approximation error
        push!(analysis, get_error(Vs, parameter))

        # Update the results vector
        push!(results, analysis)

        # Return the shifted potential (for plotting)
        return Vs 
end

function export_data(idx)
        # Extract the coefficients 
        mat = transpose(reduce(hcat, solutions))
        μ, c1, c2, c3 = eachcol(mat)

        # Construct the dataframe with its header
        df = DataFrame(parameter = μ, c1 = c1, c2 = c2, c3 = c3)

        # Export the matrix in csv
        CSV.write("../../res/data/ramped/solutions_$idx.csv", df)
 
        # Extract the transformations 
        mat = transpose(reduce(hcat, results))
        μ, LDP, error = eachcol(mat)

        # Construct the dataframe with its header
        df = DataFrame(parameter = μ, ews = LDP, error = error)

        # Export the matrix in csv
        CSV.write("../../res/data/ramped/results_$idx.csv", df)
end

function import_data(idx)
        # Import the least-squares solution of the potential reconstruction
        df = CSV.read("../../res/data/ramped/solutions_$idx.csv", DataFrame)
        μ = df.parameter
        c1 = df.c1
        c2 = df.c2
        c3 = df.c3

        # Sort the arrays according to increasing values of the parameter
        permutation_idx = sortperm(μ)
        μc_sorted = μ[permutation_idx]
        c1_sorted = c1[permutation_idx]
        c2_sorted = c2[permutation_idx]
        c3_sorted = c3[permutation_idx]
       
        # Import the ews and error timeseries
        df = CSV.read("../../res/data/ramped/results_$idx.csv", DataFrame)
        μ = df.parameter
        ews = df.ews
        error = df.error

        # Sort the arrays according to increasing values of the parameter
        permutation_idx = sortperm(μ)
        μr_sorted = μ[permutation_idx]
        ews_sorted = ews[permutation_idx]
        error_sorted = error[permutation_idx]

        return (
                coefficients = [[μc_sorted[n], c1_sorted[n], c2_sorted[n], c3_sorted[n]] for n in eachindex(μc_sorted)],
                analysis = hcat(μr_sorted, ews_sorted, error_sorted) 
               )
end
