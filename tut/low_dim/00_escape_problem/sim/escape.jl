using DynamicalSystems, DifferentialEquations
using LaTeXStrings, CairoMakie, Makie.Colors
using Statistics, ProgressMeter 

# Noise level and timescale separation
σ = 0.10
ε = 0.01
# Timestep
δt = 1e-1
T = 195.00
# IC 
x0 = [2.95,-0.35]

# Deterministic dynamics 
function iip_det!(f, x, y, t)
        f[1] = -x[2] - 2*x[1] + 3*(x[1])^2 - 0.8*(x[1])^3
        f[2] = ε
        return nothing
end
# Stochastic dynamics
function iip_stoc!(f, x, y, t)
        f[1] = +σ
        f[2] = 0.0
        return nothing
end

# Define and solve the problem 
normal_form = SDEProblem(iip_det!, iip_stoc!, x0, (0.0, T))
sol = solve(normal_form, EM(), dt=δt, verbose=false)
xt = sol[1,:]
μt = sol[2,:]

# Number of trajectories in the ensemble 
Ne = 2000
# Number of escaped trajectories
N_esc = 0
counter_escaped = zeros(Float64, length(μt))
# Number of unescaped (bounded) trajectories
N_bnd = 0

# Define arrays to store the ensemble sample paths for escaped and unescaped trajectories 
xt_bnd = fill(Float64[], 0) 
μt_bnd = fill(Float64[], 0)
xt_esc = fill(Float64[], 0)
μt_esc = fill(Float64[], 0)

using Polynomials
println("Simulating the ensemble sample paths")
# Loop over the ensemble sample paths
@showprogress for j in 1:Ne
        # Define an escape boolean variable
        escaped = false
        # Solve for the j-th trajectory
        local sol = solve(normal_form, EM(), dt=δt, verbose=false)
        # Extract the solutions
        local xt = sol[1,:] 
        local μt = sol[2,:]
        # Loop over the timesteps of the trajectory
        for t in 1:length(xt)
                # Get the state at the current timestep
                local x = xt[t]
                # Get the parameter value at the current timestep
                local μ = μt[t]
                # Construct the potential function at the current parameter value and extract its roots
                local equilibria = roots(Polynomial([-μ,-2,+3,-0.8]))
                # Get the unstable equilibria
                local x_uns = real(equilibria[2])
                # Check wheter the trajectory has overshot past the unstable equilibria
                if x < x_uns
                        escaped = true
                        for tt in t:length(xt)
                                counter_escaped[tt] = counter_escaped[tt] + 1
                        end
                        break
                end
        end
        # Update the ensemble arrays based on the escape boolean variable
        if escaped
                push!(xt_esc, xt)
                push!(μt_esc, μt)
                global N_esc = N_esc + 1
        else
                push!(xt_bnd, xt)
                push!(μt_bnd, μt)
                global N_bnd = N_bnd + 1
        end
end
println((1/Ne)*counter_escaped[end])

# Plot the ensemble's sample paths
CairoMakie.activate!(; px_per_unit = 4)
fig1 = Figure(; size = (1200, 800))#, backgroundcolor = :transparent)
ax1 = Axis(fig1[1,1:2],
    # Background
    #backgroundcolor = :transparent,
    xgridvisible = false,
    ygridvisible = false,
    limits = ((μt[1],μt[end]), (-0.5,3.25)),
    # Title
    #title = L"r = r_c = %$R",
    titlevisible = false,
    titlesize = 22,
    titlealign = :center,
    titlegap = 4.0,
    # x-axis
    xlabel = L"\mu",
    xlabelvisible = false,
    xlabelsize = 20,
    xlabelcolor = :black,
    xlabelpadding = -20.0,
    #xticks = [-1,0],
    xticksvisible = true,
    xticksize = 6,
    xticklabelsvisible = false,
    xticklabelsize = 18,
    #xtickformat = "{:.0f}",
    xscale = identity, #log10,
    xaxisposition = :bottom,
    # y-axis
    #ylabel = L"x",
    ylabelvisible = true,
    ylabelsize = 20,
    ylabelcolor = :black,
    ylabelpadding = -25.0,
    yticks = [0,1,2,3],
    yticksvisible = true,
    yticksize = 6,
    yticklabelsvisible = true,
    yticklabelsize = 18,
    ytickformat = "{:.0f}",
    yscale = identity,
    yaxisposition = :left,

)
ax2 = Axis(fig1[2,1:2],
    # Background
    #backgroundcolor = :transparent,
    xgridvisible = false,
    ygridvisible = false,
    limits = ((μt[1],μt[end]), (-0.5,3.25)),
    # Title
    #title = L"r = r_c = %$R",
    titlevisible = false,
    titlesize = 22,
    titlealign = :center,
    titlegap = 4.0,
    # x-axis
    xlabel = L"\mu",
    xlabelvisible = true,
    xlabelsize = 20,
    xlabelcolor = :black,
    xlabelpadding = -20.0,
    #xticks = [-1,0],
    xticksvisible = true,
    xticksize = 6,
    xticklabelsvisible = true,
    xticklabelsize = 18,
    #xtickformat = "{:.0f}",
    xscale = identity, #log10,
    xaxisposition = :bottom,
    # y-axis
    #ylabel = L"x",
    ylabelvisible = true,
    ylabelsize = 20,
    ylabelcolor = :black,
    ylabelpadding = -25.0,
    yticks = [0,1,2,3],
    yticksvisible = true,
    yticksize = 6,
    yticklabelsvisible = true,
    yticklabelsize = 18,
    ytickformat = "{:.0f}",
    yscale = identity,
    yaxisposition = :left,
)
# Unescaped trajectories
for n in 1:N_bnd # == size(xt,1)
        lines!(ax1, μt_bnd[n], xt_bnd[n], linewidth = 0.5, color = (:gray,0.15))
end
# Escaped trajectories
for n in 1:N_esc # == size(xt_esc,1)
        lines!(ax2, μt_esc[n], xt_esc[n], linewidth = 0.5, color = (:red,0.25))
end
# Plot the bifurcation diagram
a = μt[1]
b = μt[end]
μ=LinRange(a,b,1000)
x1 = Float64[]
μ1 = Float64[]
x2 = Float64[]
μ2 = Float64[]
x3 = Float64[]
μ3 = Float64[]
equilibria = ComplexF64[3]
using Polynomials
for n in 1:length(μ) 
        global equilibria = roots(Polynomial([-μ[n],-2,+3,-0.8]))
        if imag(equilibria[1])≈0.0
                push!(x1, real(equilibria[1]))
                push!(μ1, μ[n])
        end
        if imag(equilibria[2])≈0.0
                push!(x2, real(equilibria[2]))
                push!(μ2, μ[n])
        end
        if imag(equilibria[3])≈0.0
                push!(x3, real(equilibria[3]))
                push!(μ3, μ[n])
        end
end
lines!(ax1, μ1, x1, color = :black, linewidth = 1.5)
lines!(ax1, μ2, x2, color = :black, linestyle = :dash, linewidth = 1.5)
lines!(ax1, μ3, x3, color = :black, linewidth = 1.5)
lines!(ax2, μ1, x1, color = :black, linewidth = 1.5)
lines!(ax2, μ2, x2, color = :black, linestyle = :dash, linewidth = 1.5)
lines!(ax2, μ3, x3, color = :black, linewidth = 1.5)
save("../../results/precursors/trajectories.png", fig1)

# Plot the escape time 
function V(x, μ)
        return μ*x + x^2 - x^3 +0.2*x^4 
end

function Vxx(x)
        return 2 - 6*x + (12/5)*x^2
end

function r(μ)
        equilibria = roots(Polynomial([-μ,-2,+3,-0.8]))
        a = real(equilibria[3])
        b = real(equilibria[2])
        den = 2*pi*sqrt(abs(Vxx(a))*abs(Vxx(b)))
        ΔV = V(b,μ)-V(a,μ)
        σ = 0.05
        return [(1/den)*exp(-ΔV/σ),den,ΔV]
end
escape = [r(μt[n])[1] for n in 1:length(μt)]
den = [r(μt[n])[2] for n in 1:length(μt)]
ΔV = [r(μt[n])[3]/(σ^2/2) for n in 1:length(μt)]

CairoMakie.activate!(; px_per_unit = 3)
fig3 = Figure(; size = (1200, 1200))#, backgroundcolor = :transparent)
ax3 = Axis(fig3[1,1:2],
    # Background
    #backgroundcolor = :transparent,
    xgridvisible = false,
    ygridvisible = false,
    limits = ((μt[1], μt[end]), (-0.05,0.4)),
    # Title
    title = "Title",
    titlevisible = false,
    titlesize = 25,
    titlealign = :center,
    titlegap = -38.0,
    # x-axis
    xlabel = L"\mu",
    xlabelvisible = false,
    xlabelsize = 30,
    xlabelcolor = :black,
    xlabelpadding = 10.0,
    #xticks = [-3,-2,-1,0],
    xticksvisible = true,
    xticksize = 20,
    xticklabelsvisible = false,
    xticklabelsize = 20,
    #xtickformat = values -> [L"\sqrt{%$(round(value; digits = 1))}" for value in values],
    xscale = identity, #log10,
    xaxisposition = :bottom,
    # y-axis
    ylabel = L"\text{escape rate}",
    ylabelvisible = true,
    ylabelsize = 30,
    ylabelcolor = :black,
    ylabelpadding = 10.0,
    yticks = [0,0.1,0.3,0.4],
    yticksvisible = true,
    yticksize = 20,
    yticklabelsvisible = true,
    yticklabelsize = 20,
    #ytickformat = "{:.0f}",
    yscale = identity,
    yaxisposition = :left,
)
lines!(ax3, μt, escape, color = :black, linewidth = 4)
# Plot percentage of escaped trajectories at different y values
lines!(ax3, μt, (1/Ne).*counter_escaped, color = :red, linewidth = 2)

#=
ax4 = Axis(fig3[2,1:2],
    # Background
    #backgroundcolor = :transparent,
    xgridvisible = false,
    ygridvisible = false,
    limits = ((μt[1], μt[end]), nothing),
    # Title
    title = "Title",
    titlevisible = false,
    titlesize = 25,
    titlealign = :center,
    titlegap = -38.0,
    # x-axis
    xlabel = L"\mu",
    xlabelvisible = false,
    xlabelsize = 30,
    xlabelcolor = :black,
    xlabelpadding = 10.0,
    #xticks = [-3,-2,-1,0],
    xticksvisible = true,
    xticksize = 20,
    xticklabelsvisible = false,
    xticklabelsize = 20,
    #xtickformat = values -> [L"\sqrt{%$(round(value; digits = 1))}" for value in values],
    xscale = identity, #log10,
    xaxisposition = :bottom,
    # y-axis
    ylabel = L"2\,\pi\,\sqrt{|V^{''}(b)||V^{''}(a)|}",
    ylabelvisible = true,
    ylabelsize = 30,
    ylabelcolor = :black,
    ylabelpadding = 10.0,
    #yticks = [-0.5,0.5],
    yticksvisible = true,
    yticksize = 20,
    yticklabelsvisible = true,
    yticklabelsize = 20,
    #ytickformat = "{:.0f}",
    yscale = identity,
    yaxisposition = :left,
)
lines!(ax4, μt, den, color = :red, linewidth = 2)
=#

ax5 = Axis(fig3[2,1:2],
    # Background
    #backgroundcolor = :transparent,
    xgridvisible = false,
    ygridvisible = false,
    limits = ((μt[1], μt[end]), nothing),
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
    xlabelpadding = 10.0,
    #xticks = [-3,-2,-1,0],
    xticksvisible = true,
    xticksize = 20,
    xticklabelsvisible = true,
    xticklabelsize = 20,
    #xtickformat = values -> [L"\sqrt{%$(round(value; digits = 1))}" for value in values],
    xscale = identity, #log10,
    xaxisposition = :bottom,
    # y-axis
    ylabel = L"\frac{\Delta V}{D}",
    ylabelvisible = true,
    ylabelsize = 30,
    ylabelcolor = :black,
    ylabelpadding = 10.0,
    #yticks = [-0.5,0.5],
    yticksvisible = true,
    yticksize = 20,
    yticklabelsvisible = true,
    yticklabelsize = 20,
    #ytickformat = "{:.0f}",
    yscale = identity,
    yaxisposition = :left,
)
lines!(ax5, μt, ΔV, color = :blue, linewidth = 2)
save("../../results/precursors/escape.png", fig3)
