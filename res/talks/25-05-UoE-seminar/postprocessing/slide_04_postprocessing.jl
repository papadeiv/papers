include("../../../../inc/IO.jl")
include("../../../../inc/SystemAnalysis.jl")

# Define the normal form
f(x, μ) = -μ - 2*x + 3*(x^2) - (4/5)*(x^3)

# Import the data from csv 
u_B = readin("../data/slide_04/B-tipping.csv")
μ = u_B[:,2] 
u = u_B[:,3]

# Define arrays to store the values of the stable and unstable equilibria
x1 = Float64[]
μ1 = Float64[]
x2 = Float64[]
μ2 = Float64[]
x3 = Float64[]
μ3 = Float64[]

# Loop over the parameter values
for n in 1:length(μ)
        # Get the equilibria at the current parameter value
        eq = get_equilibria(f, μ[n], domain=[-1,4])
        # Get the stable equilibrium 
        stable = sort(eq[1], rev=true)
        # Avoid double exporting the same equilibrium
        if μ[n] < 1.6212
                push!(x1, stable[1])
                push!(μ1, μ[n])
        elseif stable[1] < 0.0
                push!(x2, stable[1])
                push!(μ2, μ[n])
        end
        # Check if you are in the region of bistability
        if length(stable)==2 
                push!(x2, stable[2])
                push!(μ2, μ[n])
        end
        # Check if the unstable equilibrium exists
        unstable = eq[2]
        if length(unstable) == 1
                push!(x3, unstable[1])
                push!(μ3, μ[n])
        end
end

# Export the data
writeout(hcat(μ1, x1), "../data/slide_04/stable_eq_1.csv")
writeout(hcat(μ2, x2), "../data/slide_04/stable_eq_2.csv")
writeout(hcat(μ3, x3), "../data/slide_04/unstable_eq.csv")
