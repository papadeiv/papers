include("EscapeProblem.jl")

μ_crit = -0.370544
μ0 = -0.365 #μ_crit + 0.0005 # -0.365
idx = 890
μ_range=LinRange(μ0,-0.371231,1000)
σ=0.010

R = Float64[] 
X = Float64[]
ΔV = Float64[] 
EXP = Float64[] 
Vxx_a = Float64[]
Vxx_b = Float64[]
Prefactor = Float64[]

for p in 1:length(μ_range)
        μ = μ_range[p]
        equilibria = get_equilibria(μ)

        rate = Kr(equilibria[1,1], equilibria[2,1], μ, σ)
        push!(R, rate)

        E = V(equilibria[2,1], μ) - V(equilibria[1,1], μ)
        push!(ΔV, E)

        x = -2*E/σ^2
        push!(X, x)
        push!(EXP, exp(x))

        curv_in_a = Vxx(equilibria[1,1], μ)
        push!(Vxx_a, curv_in_a)

        curv_in_b = abs(Vxx(equilibria[2,1], μ))
        push!(Vxx_b, curv_in_b)

        C = (1/(2*pi))*sqrt(curv_in_a*curv_in_b)
        push!(Prefactor, C)

        println("$(p): μ=$(μ)")
        println("   R = $(rate)")
        println("   ΔV = $(E)")
        println("   X = $(x)")
        println("   EXP = $(exp(-2*E/σ^2))")
        println("   C = $(C)\n")
end

using LaTeXStrings, CairoMakie, Makie.Colors
CairoMakie.activate!(; px_per_unit = 4)
fig = Figure(; size = (1200, 1800))#, backgroundcolor = :transparent)

# Escape rate and exponential 
ax1L = Axis(fig[1,1:3],
          # Background
          #backgroundcolor = :transparent,
          xgridvisible = false,
          ygridvisible = false,
          limits = ((μ_range[end],μ_range[1]), nothing),
          # Title
          title = L"R := C\,e^{-2\,\frac{\Delta V}{\sigma^2}}",
          titlevisible = true,
          titlesize = 40,
          titlealign = :center,
          titlegap = -80.0,
          # Axes labels
          xlabel = "time",
          ylabel = L"R",
          xlabelvisible = false,
          ylabelvisible = true,
          xlabelsize = 20,
          ylabelsize = 40,
          xlabelcolor = :black,
          ylabelcolor = :blue,
          xlabelpadding = 0.0,
          ylabelpadding = 20.0,
          # Axes scale, position and direction
          xscale = identity, 
          yscale = identity, #log10,
          xreversed = true,
          yreversed = false,
          xaxisposition = :bottom,
          yaxisposition = :left,
          # Ticks
          xticks = [μ_range[1], μ_range[end]],
          #yticks = equilibria[:,1],
          xticksvisible = true,
          yticksvisible = true,
          xticksize = 6,
          yticksize = 6,
          # Ticks labels
          xticklabelsvisible = false,
          yticklabelsvisible = true,
          xticklabelsize = 18,
          yticklabelsize = 18,
          xticklabelalign = (:right, :top),
          yticklabelalign = (:right, :center),
          xticklabelcolor = :black,
          yticklabelcolor = :blue,
          xtickformat = "{:.3f}",
          ytickformat = "{:.4f}",
)
lines!(ax1L, μ_range, R, color = :blue)
lines!(ax1L, [μ_range[idx], μ_range[idx]], [R[argmin(R)], R[argmax(R)]], linestyle = :dash, color = :black, linewidth = 2.5)

ax1R = Axis(fig[1,1:3],
          # Background
          #backgroundcolor = :transparent,
          xgridvisible = false,
          ygridvisible = false,
          limits = ((μ_range[end],μ_range[1]), nothing),
          # Title
          title = "Title",
          titlevisible = false,
          titlesize = 22,
          titlealign = :center,
          titlegap = 4.0,
          # Axes labels
          xlabel = "time",
          ylabel = L"e^{x}",
          xlabelvisible = false,
          ylabelvisible = true,
          xlabelsize = 20,
          ylabelsize = 40,
          xlabelcolor = :black,
          ylabelcolor = :red,
          xlabelpadding = 0.0,
          ylabelpadding = 20.0,
          # Axes scale, position and direction
          xscale = identity, 
          yscale = identity, #log10,
          xreversed = true,
          yreversed = false,
          xaxisposition = :bottom,
          yaxisposition = :right,
          # Ticks
          xticks = [μ_range[1], μ_range[end]],
          #yticks = equilibria[:,1],
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
          yticklabelalign = (:left, :center),
          xticklabelcolor = :black,
          yticklabelcolor = :red,
          xtickformat = "{:.3f}",
          ytickformat = "{:.1f}",
)
hidespines!(ax1R)
hidexdecorations!(ax1R)
lines!(ax1R, μ_range, EXP, color = :red, linewidth = 1.5)

# Energy difference and exponential
ax2L = Axis(fig[2,1:3],
          # Background
          #backgroundcolor = :transparent,
          xgridvisible = false,
          ygridvisible = false,
          limits = ((μ_range[end],μ_range[1]), nothing),
          # Title
          title = "Title",
          titlevisible = false,
          titlesize = 22,
          titlealign = :center,
          titlegap = 4.0,
          # Axes labels
          xlabel = "time",
          ylabel = L"\Delta V",
          xlabelvisible = false,
          ylabelvisible = true,
          xlabelsize = 20,
          ylabelsize = 40,
          xlabelcolor = :black,
          ylabelcolor = :blue,
          xlabelpadding = 0.0,
          ylabelpadding = 20.0,
          # Axes scale, position and direction
          xscale = identity, 
          yscale = identity, #log10,
          xreversed = true,
          yreversed = false,
          xaxisposition = :bottom,
          yaxisposition = :left,
          # Ticks
          xticks = [μ_range[1], μ_range[end]],
          #yticks = equilibria[:,1],
          xticksvisible = true,
          yticksvisible = true,
          xticksize = 6,
          yticksize = 6,
          # Ticks labels
          xticklabelsvisible = false,
          yticklabelsvisible = true,
          xticklabelsize = 18,
          yticklabelsize = 18,
          xticklabelalign = (:right, :top),
          yticklabelalign = (:right, :center),
          xticklabelcolor = :black,
          yticklabelcolor = :blue,
          xtickformat = "{:.3f}",
          #ytickformat = "{:.4f}",
)
lines!(ax2L, μ_range, ΔV, color = :blue, linewidth = 1.5)
lines!(ax2L, [μ_range[idx], μ_range[idx]], [ΔV[argmin(ΔV)], ΔV[argmax(ΔV)]], linestyle = :dash, color = :black, linewidth = 2.5)

ax2R = Axis(fig[2,1:3],
          # Background
          #backgroundcolor = :transparent,
          xgridvisible = false,
          ygridvisible = false,
          limits = ((μ_range[end],μ_range[1]), nothing),
          # Title
          title = "Title",
          titlevisible = false,
          titlesize = 22,
          titlealign = :center,
          titlegap = 4.0,
          # Axes labels
          xlabel = "time",
          ylabel = L"x:=-2\,\frac{\Delta V}{\sigma^2}",
          xlabelvisible = false,
          ylabelvisible = true,
          xlabelsize = 20,
          ylabelsize = 40,
          xlabelcolor = :black,
          ylabelcolor = :red,
          xlabelpadding = 0.0,
          ylabelpadding = 20.0,
          # Axes scale, position and direction
          xscale = identity, 
          yscale = identity, #log10,
          xreversed = true,
          yreversed = false,
          xaxisposition = :bottom,
          yaxisposition = :right,
          # Ticks
          xticks = [μ_range[1], μ_range[end]],
          #yticks = equilibria[:,1],
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
          yticklabelalign = (:left, :center),
          xticklabelcolor = :black,
          yticklabelcolor = :red,
          xtickformat = "{:.3f}",
          #ytickformat = "{:.1f}",
)
hidespines!(ax2R)
hidexdecorations!(ax2R)
lines!(ax2R, μ_range, X, color = :red, linewidth = 1.5)

# Potential function local curvature 
ax3L = Axis(fig[3,1:3],
          # Background
          #backgroundcolor = :transparent,
          xgridvisible = false,
          ygridvisible = false,
          limits = ((μ_range[end],μ_range[1]), nothing),
          # Title
          title = "Title",
          titlevisible = false,
          titlesize = 22,
          titlealign = :center,
          titlegap = 4.0,
          # Axes labels
          xlabel = L"\mu",
          ylabel = L"V^{''}(x\in\{a,b\})",
          xlabelvisible = true,
          ylabelvisible = true,
          xlabelsize = 40,
          ylabelsize = 40,
          xlabelcolor = :black,
          ylabelcolor = :blue,
          xlabelpadding = 0.0,
          ylabelpadding = 20.0,
          # Axes scale, position and direction
          xscale = identity, 
          yscale = identity, #log10,
          xreversed = true,
          yreversed = false,
          xaxisposition = :bottom,
          yaxisposition = :left,
          # Ticks
          xticks = [μ_range[1], μ_range[end]],
          #yticks = equilibria[:,1],
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
          yticklabelcolor = :blue,
          xtickformat = "{:.6f}",
          #ytickformat = "{:.2f}",
)
lines!(ax3L, μ_range, Vxx_a, color = :blue)
lines!(ax3L, μ_range, Vxx_b, color = :blue, linestyle = :dash)
lines!(ax3L, [μ_range[idx], μ_range[idx]], [Vxx_a[argmin(Vxx_a)], Vxx_a[argmax(Vxx_a)]], linestyle = :dash, color = :black, linewidth = 2.5)

ax3R = Axis(fig[3,1:3],
          # Background
          #backgroundcolor = :transparent,
          xgridvisible = false,
          ygridvisible = false,
          limits = ((μ_range[end],μ_range[1]), nothing),
          # Title
          title = "Title",
          titlevisible = false,
          titlesize = 22,
          titlealign = :center,
          titlegap = 4.0,
          # Axes labels
          xlabel = L"\mu",
          ylabel = L"C:=\frac{1}{2\pi}\sqrt{V^{''}(a)V^{''}(b)}",
          xlabelvisible = true,
          ylabelvisible = true,
          xlabelsize = 40,
          ylabelsize = 40,
          xlabelcolor = :black,
          ylabelcolor = :red,
          xlabelpadding = 0.0,
          ylabelpadding = 20.0,
          # Axes scale, position and direction
          xscale = identity, 
          yscale = identity, #log10,
          xreversed = true,
          yreversed = false,
          xaxisposition = :bottom,
          yaxisposition = :right,
          # Ticks
          xticks = [μ_range[1], μ_range[end]],
          #yticks = equilibria[:,1],
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
          yticklabelalign = (:left, :center),
          xticklabelcolor = :black,
          yticklabelcolor = :red,
          xtickformat = "{:.3f}",
          #ytickformat = "{:.3f}",
)
hidespines!(ax3R)
hidexdecorations!(ax3R)
lines!(ax3R, μ_range, Prefactor, color = :red)

save("../../results/precursors/escape_problem/hitting_times/KramersFormula.png", fig)
