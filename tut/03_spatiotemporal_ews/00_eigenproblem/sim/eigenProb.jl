using ApproxFun 

# Spatial domain of the operator
dx = -32.0*π..32.0*π
x = Fun(dx)
S = space(x)

# Boundary conditions
D = Dirichlet(S)
N = Neumann(S)
P = periodic(S,0)

# Define the spatial derivative
∂x = Derivative(dx)

# Number of grid nodes
n = 10#00

# Vector values of the parameter
Ω = LinRange(-1,1,2)

# Initialise array of leading eigenvalues
λmax = Float64[] 

# Solve the eigenproblem at different parameter values
using ProgressMeter
@showprogress for µ in Ω
        # Linearisation of the non-linear reaction term about different steady-states
        local Df = µ - 1.0
        # Construction of the linear operator
        local L =  -∂x^4 - 2.0*∂x^2 + Df
        # Solution of the eigenproblem
        local λ, v = ApproxFun.eigs(D, L, n, tolerance=1e-10);
        # Get the datatype of v
        display(typeof(v))
        # Extraction of the real part of the point spectrum
        local σ = sort(real(λ), rev=true)
        # Extract the leading eigenvalue from the point spectrum
        push!(λmax, σ[1])
end
println(λmax)

# Plot of the leading eigenvalue at different parameter values 
using CairoMakie, Makie.Colors
CairoMakie.activate!(; px_per_unit = 2)
fig = Figure(; size = (1200, 800), backgroundcolor = :transparent)
ax = Axis(fig[1, 1],
    ylabel = L"\lambda_{\text{max}}",
    backgroundcolor = :transparent,
    #title = latexstring("u_0(x) = \\text{Train of pulses}, \\quad t=", @sprintf("%.3f", t_h[i]), ", \\quad µ =", @sprintf("%.3f", µ_h[1,i])),
    xticksvisible = false,
    #xticklabelsvisible = false,
    xgridvisible = false,
    yticksize = 10,
    yticklabelsize = 24,
    xticksize = 10,
    xticklabelsize = 24,
    yticksvisible = false,
    #yticklabelsvisible = false,
    ygridvisible = false,
    #limits = (a, b, -2.0, 2.0),
    xlabel = L"\mu",
    xlabelsize = 32,
    xlabelcolor = :black,
    #xlabelpadding = -20.0,
    xticks = [-1,0,1],
    #ylabel = L"u",
    ylabelsize = 32,
    ylabelcolor = :black,
    #ylabelpadding = -20.0,
    #yticks = [0.0],
)
scatter!(ax, Ω, λmax, color = :blue, strokecolor = :black, strokewidth = 1.5, markersize = 30)
save("./test.png", fig)
