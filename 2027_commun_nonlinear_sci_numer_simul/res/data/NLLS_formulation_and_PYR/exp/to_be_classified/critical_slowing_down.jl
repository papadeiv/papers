using DynamicalSystems, DifferentialEquations
using LaTeXStrings, CairoMakie, Makie.Colors
using Statistics, ProgressMeter 

# Noise level and timescale separation
σ = 0.0
ε = 0.01
# Perturbation
δx = -0.25
# Timestep
δt = 1e-3
# Return-time array
CSD_t = Vector{Float64}()
CSD_μ = Vector{Float64}()

# Deterministic dynamics 
function iip_det!(f, x, y, t)
        f[1] = x[2]*x[1]-(x[1])^2
        f[2] = ε
        return nothing
end
# Stochastic dynamics
function iip_stoc!(f, x, y, t)
        f[1] = +σ
        f[2] = 0.0
        return nothing
end

# First IC
x0 = [0.0,-3.0]
T = 30.00
# Define and solve the problem 
normal_form = SDEProblem(iip_det!, iip_stoc!, x0, (0.0, T))
sol = solve(normal_form, EM(), dt=δt, verbose=false)
xt = sol[1,:]
μt = sol[2,:]

# Second IC (endstate of the first run)
x0 = [xt[end]+δx, μt[end]]
T = 50
# Define and solve the problem 
normal_form = SDEProblem(iip_det!, iip_stoc!, x0, (0.0, T))
sol = solve(normal_form, EM(), dt=δt, verbose=false)
x = sol[1,:]
μ = sol[2,:]
# Compute the return time
time_idx = 0
for n in 1:length(x)
        if abs(x[n]) > 1e-5
                global time_idx = time_idx + 1
        end
end
push!(CSD_t, time_idx*δt)
push!(CSD_μ, μ[1])
# Chain together the solutions
xt = [xt; x]
μt = [μt; μ]

# Third IC (endstate of the second run)
x0 = [xt[end]+δx, μt[end]]
# Define and solve the problem 
normal_form = SDEProblem(iip_det!, iip_stoc!, x0, (0.0, T))
sol = solve(normal_form, EM(), dt=δt, verbose=false)
x = sol[1,:]
μ = sol[2,:]
# Compute the return time
time_idx = 0
for n in 1:length(x)
        if abs(x[n]) > 1e-5
                global time_idx = time_idx + 1
        end
end
push!(CSD_t, time_idx*δt)
push!(CSD_μ, μ[1])
# Chain together the solutions
xt = [xt; x]
μt = [μt; μ]

# Fourth IC (endstate of the third run)
x0 = [xt[end]+δx, μt[end]]
T = 60
# Define and solve the problem 
normal_form = SDEProblem(iip_det!, iip_stoc!, x0, (0.0, T))
sol = solve(normal_form, EM(), dt=δt, verbose=false)
x = sol[1,:]
μ = sol[2,:]
# Compute the return time
time_idx = 0
for n in 1:length(x)
        if abs(x[n]) > 1e-5
                global time_idx = time_idx + 1
        end
end
push!(CSD_t, time_idx*δt)
push!(CSD_μ, μ[1])
# Chain together the solutions
xt = [xt; x]
μt = [μt; μ]

# Fifth IC (endstate of the fourth run)
x0 = [xt[end]+δx, μt[end]]
# Define and solve the problem 
normal_form = SDEProblem(iip_det!, iip_stoc!, x0, (0.0, T))
sol = solve(normal_form, EM(), dt=δt, verbose=false)
x = sol[1,:]
μ = sol[2,:]
# Compute the return time
time_idx = 0
for n in 1:length(x)
        if abs(x[n]) > 1e-5
                global time_idx = time_idx + 1
        end
end
push!(CSD_t, time_idx*δt)
push!(CSD_μ, μ[1])
# Chain together the solutions
xt = [xt; x]
μt = [μt; μ]

# Sixth IC (endstate of the fifth run)
x0 = [xt[end]+δx, μt[end]]
# Define and solve the problem 
normal_form = SDEProblem(iip_det!, iip_stoc!, x0, (0.0, T))
sol = solve(normal_form, EM(), dt=δt, verbose=false)
x = sol[1,:]
μ = sol[2,:]
# Compute the return time
time_idx = 0
for n in 1:length(x)
        if abs(x[n]) > 1e-5
                global time_idx = time_idx + 1
        end
end
push!(CSD_t, time_idx*δt)
push!(CSD_μ, μ[1])
# Chain together the solutions
xt = [xt; x]
μt = [μt; μ]

# Plot the solutions
CairoMakie.activate!(; px_per_unit = 3)
fig = Figure(; size = (1200, 800))#, backgroundcolor = :transparent)
ax1 = Axis(fig[1,1:4],
    # Background
    backgroundcolor = :transparent,
    xgridvisible = false,
    ygridvisible = false,
    limits = ((-3,0),(-0.5,0.5)),
    # Title
    title = "Title",
    titlevisible = false,
    titlesize = 25,
    titlealign = :center,
    titlegap = -38.0,
    # x-axis
    xlabel = L"time",
    xlabelvisible = false,
    xlabelsize = 30,
    xlabelcolor = :black,
    xlabelpadding = -32.0,
    xticks = [-3,-2,-1,0],
    xticksvisible = true,
    xticksize = 20,
    xticklabelsvisible = false,
    xticklabelsize = 28,
    #xtickformat = values -> [L"\sqrt{%$(round(value; digits = 1))}" for value in values],
    xscale = identity, #log10,
    xaxisposition = :bottom,
    # y-axis
    ylabel = L"x",
    ylabelvisible = true,
    ylabelsize = 30,
    ylabelcolor = :black,
    ylabelpadding = -38.0,
    yticks = [-0.5,0.5],
    yticksvisible = true,
    yticksize = 20,
    yticklabelsvisible = true,
    yticklabelsize = 23,
    ytickformat = "{:.1f}",
    yscale = identity,
    yaxisposition = :left,
)
lines!(ax1, [-3.0,0.0], [0.0,0.0], color = :black, linewidth = 1)
lines!(ax1, μt, xt, linewidth = 2.0, color = :red)

# Plot the return-time and the eigenvalue at the equilibrium
ax2 = Axis(fig[2,1:4],
    # Background
    backgroundcolor = :transparent,
    xgridvisible = false,
    ygridvisible = false,
    limits = ((-3,0), (nothing, nothing)),
    # Title
    title = "Title",
    titlevisible = false,
    titlesize = 25,
    titlealign = :center,
    titlegap = -38.0,
    # x-axis
    xlabel = L"\mu",
    xlabelvisible = true,
    xlabelsize = 30,
    xlabelcolor = :black,
    xlabelpadding = -32.0,
    xticks = [-3,-2,-1,0],
    xticksvisible = true,
    xticksize = 20,
    xticklabelsvisible = true,
    xticklabelsize = 23,
    #xtickformat = values -> [L"\sqrt{%$(round(value; digits = 1))}" for value in values],
    xscale = identity, #log10,
    xaxisposition = :bottom,
    # y-axis
    ylabel = L"\text{recovery time}",
    ylabelvisible = true,
    ylabelsize = 30,
    ylabelcolor = :black,
    ylabelpadding = -28.0,
    yticks = [4,32],
    yticksvisible = true,
    yticksize = 20,
    yticklabelsvisible = true,
    yticklabelsize = 23,
    ytickformat = "{:.0f}",
    yscale = identity,
    yaxisposition = :left,
)
scatter!(ax2, CSD_μ, CSD_t, color = :blue, strokecolor = :black, strokewidth = 1.5, markersize = 15)
# Export figures
save("../../results/precursors/fig1.2.png", fig)
