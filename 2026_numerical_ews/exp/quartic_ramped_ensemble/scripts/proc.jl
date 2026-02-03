"""
    Postprocessing script

Collection of quantities and functions used to postprocess and analyse the results of a simulation.
"""

# Parameters of the scalar potential method
Na = convert(Int64, 1e4)                        # Number of attempts per guess 
β = 1e-2                                        # Std of the guess perturbation 

# Scalar potential of the conservative system 
U(x, μ) = 1 + μ*x - 2*(x^2) + x^4               # Ground truth
V(x, c) = c[1]*x + c[2]*(x^2) + c[3]*(x^3)      # Arbitrary cubic potential

# Data structures
solutions = Vector{Vector{Float64}}()           # Solutions of the inference method 
results = Vector{Vector{Float64}}()             # Statistical distribution over the ensemble
glb_idx = 1::Integer                            # Global index to loop over the parameter values

function shift_potential(x0, μ, c)
        # Compute the stable equilibrium (center of the shift)
        xs = +(1/(3*c[3]))*(sqrt((c[2])^2 - 3*c[1]*c[3]) - c[2])

        # Compute the shifts
        δx = x0 - xs 
        δy = U(x0, μ) - (Polynomial([0.0; c]))(xs)

        # Define the shifted potential
        Vs(x) = δy + c[1]*(x - δx) + c[2]*(x - δx)^2 + c[3]*(x - δx)^3

        return xs, Vs
end

function analyse(solutions, parameter)
        # Compute estimated stable and unstable equilibria of the cubic
        xs = +(1/(3*solutions[3]))*(sqrt((solutions[2])^2 - 3*solutions[1]*solutions[3]) - solutions[2])
        xu = -(1/(3*solutions[3]))*(sqrt((solutions[2])^2 - 3*solutions[1]*solutions[3]) + solutions[2])

        # Define empty vector
        analysis = Float64[]
        push!(analysis, parameter)

        # Compute the ews
        ΔV = abs(V(xu, solutions) - V(xs, solutions))
        push!(analysis, exp(-ΔV))

        # Update the results vector
        push!(results, analysis)
end

function export_data()
        # Extract the coefficients 
        mat = transpose(reduce(hcat, solutions))
        xs, μ, c1, c2, c3 = eachcol(mat)

        # Construct the dataframe with its header
        df = DataFrame(parameter = μ, qse = xs, c1 = c1, c2 = c2, c3 = c3)

        # Export the matrix in csv
        CSV.write("../../res/data/ramped_ensemble/solutions/$glb_idx.csv", df)
 
        # Extract the transformations 
        mat = transpose(reduce(hcat, results))
        μ, LDP = eachcol(mat)

        # Construct the dataframe with its header
        df = DataFrame(parameter = μ, ews = LDP)

        # Export the matrix in csv
        CSV.write("../../res/data/ramped_ensemble/results/$glb_idx.csv", df)
end

function import_data()
        # Import the least-squares solution of the potential reconstruction
        df = CSV.read("../../res/data/ramped_ensemble/solutions/$glb_idx.csv", DataFrame)
        μc = df.parameter
        xs = df.qse
        c1 = df.c1
        c2 = df.c2
        c3 = df.c3
      
        # Import the ews and error timeseries
        df = CSV.read("../../res/data/ramped_ensemble/results/$glb_idx.csv", DataFrame)
        μr = df.parameter
        ews = df.ews
        error = df.error

        return (
                coefficients = [[μc[n], c1[n], c2[n], c3[n]] for n in eachindex(μc)],
                analysis = hcat(μr, ews, error) 
               )
end
