include("../../../../inc/IO.jl")
include("../../../../inc/SystemAnalysis.jl")


# Define the dynamics of the frozen system 
a = -(0.25::Float64)
b = 1.20::Float64
c = -(0.40::Float64)
d = -(0.30::Float64)
e = 3.00::Float64
K = 2.00::Float64
f(x, λ) = -((x + a + b*λ)^2 + c*tanh(λ - d))*(x - K/(cosh(e*λ)))

# Import the data from csv 
solution = readin("../data/figure_01/solutions/1.csv")
t = solution[:,1] 
μ = solution[:,2] 
u = solution[:,3]

ε = readin("../data/figure_01/rates.csv")
Nε = length(ε)

# Define arrays to store the values of the stable and unstable equilibria
x1 = Float64[]
μ1 = Float64[]
x2 = Float64[]
μ2 = Float64[]
x3 = Float64[]
μ3 = Float64[]

# Loop over the parameter values
printstyled("Computing the bifurcation diagram of the frozen system\n"; bold=true, underline=true, color=:light_blue)
@showprogress for n in 1:length(μ)
        # Get the equilibria at the current parameter value
        eq = get_equilibria(f, μ[n], domain=[-10,+10])
        # Get the stable equilibrium 
        stable = sort(eq[1], rev=true)

        # Check if you are in the region of bistability
        if length(stable)==2 
                push!(x1, stable[1])
                push!(μ1, μ[n])
                push!(x2, stable[2])
                push!(μ2, μ[n])
        else
                push!(x1, stable[1])
                push!(μ1, μ[n])
        end
        # Check if the unstable equilibrium exists
        unstable = eq[2]
        if length(unstable) == 1
                push!(x3, unstable[1])
                push!(μ3, μ[n])
        end
end

# Define a neighbourhood of the Λ-connected endstate X+
X = x1[end]  
Y = x3[end]
δ = 1e-1
U = [Y, X+δ] 

# Empty vector to store the rates for which R-tipping is observed
εc = Float64[]
idx = Int64[]

# Loop over the rates
printstyled("Identifying the region of irreversible R-tipping\n"; bold=true, underline=true, color=:light_blue)
@showprogress for n in 1:Nε
        # Import the data from csv 
        local solution = readin("../data/figure_01/solutions/$n.csv")
        x = solution[:,3]

        # Check for R-tipping
        if x[end] > U[end] || x[end] < U[1]
                push!(idx, n)
                push!(εc, ε[n])
        end
end

# Export the data
writeout(hcat(μ1, x1), "../data/figure_01/stable_eq_1.csv")
writeout(hcat(μ2, x2), "../data/figure_01/stable_eq_2.csv")
writeout(hcat(μ3, x3), "../data/figure_01/unstable_eq.csv")
writeout(hcat(idx, εc), "../data/figure_01/critical_rates.csv")
