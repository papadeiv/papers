using DynamicalSystems
using LaTeXStrings, CairoMakie

# Define the dynamics
function iip_sn!(f, x, p, t)
        f[1] = (1/p[1])*(x[2] + x[3] + x[1]*(x[1]-1))
        f[2] = -(x[1] + (x[1])^2 + (x[1])^3 + (x[1])^4 + (x[1])^5)
        f[3] = p[2] 
        return nothing
end

# Define the set of ICs 
x0 = [[0.00, -0.50, 0.10], [1.00, -1.00, 0.10]]
# Define the parameter values
ε = 0.02
R = 1/2 + 1/4 + 1/8 + 1/16 + 1/32
r = 1.2
p = [ε, r]
# Define the temporal parameters
T1 = 100.00
δt = 1e-2

# Evolve the dynamical system from the first initial condition
saddle_node = ContinuousDynamicalSystem(iip_sn!, x0[2], p)
Xt, t = trajectory(saddle_node, T1; Δt=δt)
x = Xt[:,1]
y = Xt[:,2]
μ = Xt[:,3]

# Plot the time trajectory for x(t)
CairoMakie.activate!(; px_per_unit = 3)
fig1 = Figure(; size = (600, 400), backgroundcolor = :transparent)
ax1 = Axis(fig1[1, 1],
    backgroundcolor = :transparent,
    xlabel = L"t",
    ylabel = L"x(t)",
    limits = ((t[1],t[end]), (nothing,nothing))
)
lines!(ax1, t, x, color = :green, linewidth = 1.5)
# Plot the time trajectory for y(t)
ax2 = Axis(fig1[2, 1],
    backgroundcolor = :transparent,
    xlabel = L"t",
    ylabel = L"y(t)",
    limits = ((t[1],t[end]), (nothing,nothing))
)
lines!(ax2, t, y, color = :green, linewidth = 1.5)

# figure setup 
fig2 = Figure(; size = (600, 400), backgroundcolor = :transparent)
ax3 = Axis(fig2[1, 1],
    # Background
    backgroundcolor = :transparent,
    xgridvisible = false,
    ygridvisible = false,
    limits = ((-1.50,2.50), (-3.00,1.00)),
    # Title
    title = L"r = %$r \; > \; r_c = %$R",
    titlevisible = true,
    titlesize = 22,
    titlealign = :center,
    titlegap = 4.0,
    # x-axis
    xlabel = L"x",
    xlabelvisible = true,
    xlabelsize = 24,
    xlabelcolor = :black,
    xlabelpadding = -20.0,
    xticks = [-1.5,-1,-0.5,0,1,1.5,2,2.5],
    xticksvisible = true,
    xticksize = 6,
    xticklabelsvisible = true,
    xticklabelsize = 14,
    xtickformat = "{:.1f}",
    xscale = identity, #log10,
    xaxisposition = :bottom,
    # y-axis
    ylabel = L"y",
    ylabelvisible = false,
    ylabelsize = 24,
    ylabelcolor = :black,
    ylabelpadding = -15.0,
    yticks = [-3,-2,0,1],
    yticksvisible = true,
    yticksize = 6,
    yticklabelsvisible = false,
    yticklabelsize = 14,
    ytickformat = "{:.0f}",
    yscale = identity,
    yaxisposition = :left,
)
scatter!(ax3, x0[2][1], x0[2][2], color = :blue, strokecolor = :black, strokewidth = 1.5, markersize = 9)
text!(-1.1, 0.5, text = L"\text{(d)}", align = [:right, :bottom], color = :black, fontsize = 26)

# Plot the critical manifold at different parameter values 
x_values = LinRange(-10.0,10.0,1000)
y1_values = -μ[1] .- x_values.*(x_values .- 1.0)
lines!(ax3, x_values, y1_values, color = :black, linewidth = 0.75)
y2_values = -μ[100] .- x_values.*(x_values .- 1.0)
lines!(ax3, x_values, y2_values, color = :black, linestyle = :dash, linewidth = 0.75)
y3_values = -μ[200] .- x_values.*(x_values .- 1.0)
lines!(ax3, x_values, y3_values, color = :black, linestyle = :dash, linewidth = 0.75)
y4_values = -μ[300] .- x_values.*(x_values .- 1.0)
lines!(ax3, x_values, y4_values, color = :black, linestyle = :dash, linewidth = 0.75)
# Compute and plot the rate-dependent QSE: the numerical value is the (unique real) root of x^5+x^4+x^3+x^2+x-r=0
using Polynomials
all_roots = roots(Polynomial([-r,1,1,1,1,1]))
qse = 0.0
for n in 1:length(all_roots)
        if imag(all_roots[n])==0.0
                global qse = real(all_roots[n])
        end
end
invariant = qse.*ones(Float64, 1000)
lines!(ax3, invariant, x_values, color = :green, linewidth = 2.0)
# Plot the fold of the manifold
fold = 0.5.*ones(Float64, 1000)
lines!(ax3, fold, x_values, color = :red, linewidth = 1.0)
# Plot the trajectory in state space
lines!(ax3, x, y, color = :blue, linewidth = 1.5)

# Export the figures
save("../../../results/fast_slow/CompBombTraj.png", fig1)
save("../../../results/fast_slow/fig2.2.4.png", fig2)

#=
# Create animation
fig = Figure(; size = (600,400))
ax = Axis(fig[1,1],
          title = L"r = 1.2 > r_c = %$R",
          xlabel = L"x_1",
          ylabel = L"x_2",
          limits = ((-2.00,2.00), (-2.00,2.00))
)
scatter!(ax, x0[1][1], x0[1][2], color = :blue, strokecolor = :black, strokewidth = 1.5, markersize = 9)
scatter!(ax, x0[2][1], x0[2][2], color = :blue, strokecolor = :black, strokewidth = 1.5, markersize = 9)
lines!(ax, x1_values, x2_values, color = :black, linewidth = 0.75, linestyle = :dash)
lines!(ax, invariant, x1_values, color = :red, linewidth = 1.0)
lines!(ax, qse, x1_values, color = :black, linestyle = :dash, linewidth = 1.0)

frames = 1:size(x1)[1] 
framerate = 60

points1 = Observable(Point2f[(x0[1][1],x0[1][2])])
points2 = Observable(Point2f[(x0[2][1],x0[2][2])])
time = Observable(x0[1][3])
x2_values = @lift(-$time .- x1_values.*(x1_values .- 1.0))
lines!(ax, x1_values, x2_values, color = :black, linewidth = 0.75)
lambda = @lift(-$time)
scatter!(ax, 0.0, lambda, color = :black, markersize = 9)
text!(ax, 0.95, 0.95, text = @lift("λ = " .* string.(round($time + 0.0, digits = 3))))

record(fig, "../../../results/CompBombState.gif", frames; framerate=framerate) do frame 
       new_point1 = Point2f(x1[frame],x2[frame])
       new_point2 = Point2f(x4[frame],x5[frame])
       points1[] = push!(points1[], new_point1)
       points2[] = push!(points2[], new_point2)
       time[] = x3[frame]
       lines!(ax, points1, color = :green, linewidth = 1.5)
       lines!(ax, points2, color = :orange, linewidth = 1.5)
end
=#
