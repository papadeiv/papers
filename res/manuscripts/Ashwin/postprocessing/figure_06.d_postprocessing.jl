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

# Define the parameter shift
s(t) = tanh(t)
C = [0.0, 0.0, 0.0, 0.0, -1.6, 0.0, 2.0, 0.0, 0.8, 0.0, -0.5]
Λ(t) = C[1]*(s(t)) + C[2]*(s(t))^2 + C[3]*(s(t))^3 + C[4]*(s(t))^4 + C[5]*(s(t))^5 + C[6]*(s(t))^6 + C[7]*(s(t))^7 + C[8]*(s(t))^8 + C[9]*(s(t))^9 + C[10]*(s(t))^10 + C[11]*(s(t))^11

# Import the data from csv 
solution = readin("../data/figure_06/solution.csv")
t = solution[:,1] 
μ = solution[:,2] 
u = solution[:,3]

# Define arrays to store the values of the stable and unstable equilibria
x1 = Float64[]
t1 = Float64[]
x2 = Float64[]
t2 = Float64[]
x3 = Float64[]
t3 = Float64[]

# Loop over the parameter values
printstyled("Computing the bifurcation diagram of the frozen system\n"; bold=true, underline=true, color=:light_blue)
@showprogress for n in 1:length(t)
        # Get the equilibria at the current parameter value
        eq = get_equilibria(f, Λ(t[n]), domain=[-10,+10])
        # Get the stable equilibrium 
        stable = sort(eq[1], rev=true)

        # Check if you are in the region of bistability
        if length(stable)==2 
                push!(x1, stable[1])
                push!(t1, t[n])
                push!(x2, stable[2])
                push!(t2, t[n])
        else
                push!(x1, stable[1])
                push!(t1, t[n])
        end
        # Check if the unstable equilibrium exists
        unstable = eq[2]
        if length(unstable) == 1
                push!(x3, unstable[1])
                push!(t3, t[n])
        end
end

# Export the data
writeout(hcat(t1, x1), "../data/figure_06/stable_eq_1.csv")
writeout(hcat(t2, x2), "../data/figure_06/stable_eq_2.csv")
writeout(hcat(t3, x3), "../data/figure_06/unstable_eq.csv")
