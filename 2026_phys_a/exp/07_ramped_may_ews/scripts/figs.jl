"""
    Figures layout

Generation of the layouts and formats of the figures.
"""

# Specify the figure dimensions and style
fig = Figure(; size = (1300, 650), figure_padding = (30,30,30,30))

# Set thickness and size of axes elements 
border = 2.0
labels = 20
ticks = 22

# Set boundaries of the potential range
upper = 4
lower = -8

# Variance EWS 
ax1 = Axis(fig[1,1],
           spinewidth = border,
           xgridvisible = false,
           ygridvisible = false,
           xlabelvisible = true,
           ylabelvisible = true,
           xticklabelsvisible = true,
           yticklabelsvisible = true,
           xtickalign = 1,
           ytickalign = 1,
           xtickformat = values -> ["$(trunc(value, digits=1))" for value in values],
           ytickformat = values -> ["$(trunc(value, digits=1))" for value in values],
           xtickwidth = border,
           ytickwidth = border,
           xlabel = L"\mu",
           ylabel = L"\text{Var}(x)",
           xlabelsize = labels,
           ylabelsize = labels,
           xticklabelsize = labels,
           yticklabelsize = labels,
          )

# Escape EWS 
ax2 = Axis(fig[1,2],
           spinewidth = border,
           xgridvisible = false,
           ygridvisible = false,
           xlabelvisible = true,
           ylabelvisible = true,
           xticklabelsvisible = true,
           yticklabelsvisible = true,
           xtickalign = 1,
           ytickalign = 1,
           xtickformat = values -> ["$(trunc(value, digits=1))" for value in values],
           ytickformat = values -> ["$(trunc(value, digits=1))" for value in values],
           xtickwidth = border,
           ytickwidth = border,
           xlabel = L"\mu",
           ylabel = L"\exp(-\Delta V)",
           xlabelsize = labels,
           ylabelsize = labels,
           xticklabelsize = labels,
           yticklabelsize = labels,
          )

# Compute and lot the analytical EWS (escape and variance)
domain = LinRange(μ0, μf, 1000)
escape_ews = Vector{Float64}(undef, length(domain))
variance_ews = Vector{Float64}(undef, length(domain)) 
for (index, μ) in enumerate(domain)
        # Compute the equilibria
        equilibria = get_equilibria(f, μ, domain=[0,10])
        a = maximum(equilibria.stable)
        b = maximum(equilibria.unstable)
        
        # Compute the variance of the stationary solution of the FPE
        I = (b, Inf)
        pdf_integral = IntegralProblem(ρ, I, μ)
        Z = (solve(pdf_integral, QuadGKJL(; order=20000); maxiters=10000)).u
        p(x, μ) = ρ(x, μ)/Z
        v(x, μ) = p(x, μ)*(x - a)^2 
        var_integral = IntegralProblem(v, I, μ)
        variance_ews[index] = (solve(var_integral, QuadGKJL(; order=20000); maxiters=10000)).u

        # Compute the modified escape rate
        ΔV = U(b, μ) - U(a, μ)
        escape_ews[index] = exp(-ΔV)
end
lines!(ax1, domain, variance_ews, color = :black, linewidth = 3.0)
lines!(ax2, domain, escape_ews, color = :black, linewidth = 3.0)
