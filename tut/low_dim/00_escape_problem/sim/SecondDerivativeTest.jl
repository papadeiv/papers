using LaTeXStrings, CairoMakie, Makie.Colors
include("../ApproximatePotential.jl")

x = LinRange(-2,2,300)
f = [y^2 for y in x]
fx = [2*y for y in x]
fxx = [2*(0*y+1) for y in x]

g = f + randn(300)./100000
dg = approx_derivative(x, g)
dgg = approx_derivative(x[1:(end-1)], dg)

# Create the figure 
CairoMakie.activate!(; px_per_unit = 2)
fig = Figure(; size = (1600, 1200))#, backgroundcolor = :transparent)

ax = Axis(fig[1,1],
                # Background
                #backgroundcolor = :transparent,
                xgridvisible = false,
                ygridvisible = false,
                limits = ((x[1],x[end]), nothing),
                # Title
                #title = L"\mu = %$parameter",
                titlevisible = true,
                titlesize = 22,
                titlealign = :center,
                titlegap = 4.0,
                # Axes labels
                xlabel = "time",
                ylabel = L"x",
                xlabelvisible = false,
                ylabelvisible = false,
                xlabelsize = 20,
                ylabelsize = 20,
                xlabelcolor = :black,
                ylabelcolor = :black,
                xlabelpadding = 0.0,
                ylabelpadding = 0.0,
                # Axes scale, position and direction
                xscale = identity, 
                yscale = identity, #log10,
                xreversed = false,
                yreversed = false,
                xaxisposition = :bottom,
                yaxisposition = :left,
                # Ticks
                #xticks = equilibria[:,1],
                #yticks = [time[1], time[end]],
                xticksvisible = true,
                yticksvisible = true,
                xticksize = 6,
                yticksize = 6,
                # Ticks labels
                xticklabelsvisible = true,
                yticklabelsvisible = true,
                xticklabelsize = 18,
                yticklabelsize = 18,
                xticklabelalign = (:right, :top),
                yticklabelalign = (:right, :center),
                xticklabelcolor = :black,
                yticklabelcolor = :black,
                xtickformat = "{:.2f}",
                ytickformat = "{:.3f}",
)
lines!(ax, x, f) 
lines!(ax, x, fx) 
lines!(ax, x, fxx) 
scatter!(ax, x[1:end], g, color = (:blue, 0.5)) 
scatter!(ax, x[1:(end-1)], dg, color = (:blue, 0.5)) 
scatter!(ax, x[1:(end-2)], dgg, color = (:red, 0.5)) 
save("./SecondDerivativeTest.png", fig)
