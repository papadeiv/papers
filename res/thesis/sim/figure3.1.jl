using DynamicalSystems
using LaTeXStrings, CairoMakie, Makie.Colors

# FAST-SLOW EQUATION

# Define the dynamics
function iip_bfold!(f, x, y, t)
        f[1] = -x[2]-(x[1])^2
        f[2] = 0.01
        return nothing
end

# Define the initial state and final time
x0 = [-0.50, -1.00]
T = 110.80
δt = 1e-1
y_values = LinRange(-1.125,0.00,200)

# Evolve the dynamical system from different initial conditions 
bfold = ContinuousDynamicalSystem(iip_bfold!, x0, nothing)
Xt, t = trajectory(bfold, T; Δt=δt)
X0 = Xt[:,1]
Y0 = Xt[:,2]

# Plot the (deterministic) critical manifold
CairoMakie.activate!(; px_per_unit = 3)
fig1 = Figure(; size = (600, 400), backgroundcolor = :transparent)
ax = Axis(fig1[1, 1],
    # Background
    backgroundcolor = :transparent,
    xgridvisible = false,
    ygridvisible = false,
    limits = ((y_values[1],y_values[end]+0.15), (-1.25,+1.25)),
    # Title
    title = "Title",
    titlevisible = false,
    titlesize = 25,
    titlealign = :center,
    titlegap = -38.0,
    # x-axis
    xlabel = L"\mu",
    xlabelvisible = true,
    xlabelsize = 32,
    xlabelcolor = :black,
    xlabelpadding = -35.0,
    xticks = [-1,0],
    xticksvisible = true,
    xticksize = 10,
    xticklabelsvisible = true,
    xticklabelsize = 25,
    xtickformat = "{:.0f}",
    xscale = identity, #log10,
    xaxisposition = :bottom,
    # y-axis
    ylabel = L"x",
    ylabelvisible = true,
    ylabelsize = 32,
    ylabelcolor = :black,
    ylabelpadding = -40.0,
    yticks = [-1,-0.5,0.5,1],
    yticksvisible = true,
    yticksize = 10,
    yticklabelsvisible = true,
    yticklabelsize = 25,
    ytickformat = "{:.1f}",
    yscale = identity,
    yaxisposition = :left,
)
stable = sqrt.(-y_values)
unstable = -sqrt.(-y_values)
lines!(ax, y_values, stable, color = :black, linewidth = 1.5)
lines!(ax, y_values, unstable, color = :black, linewidth = 1.5, linestyle = :dash)
text!(0, 0.9, text = L"\text{(a)}", align = [:left, :bottom], color = :black, fontsize = 30)

# SINGULAR LIMIT - Trajectory 1

# Define the dynamics
function iip_fs1!(f, x, y, t)
        f[1] = -x[2]-(x[1])^2
        f[2] = 0.0 
        return nothing
end

# Define the initial state and final time
x0 = [-0.50,-1.00]
T = 10.0
δt = 1e-1

# Evolve the dynamical system from different initial conditions 
bfold = ContinuousDynamicalSystem(iip_fs1!, x0, nothing)
Xt, t = trajectory(bfold, T; Δt=δt)
X1 = Xt[:,1]
Y1 = Xt[:,2]

# SINGULAR LIMIT - Trajectory 2

# Define the dynamics
function iip_fs2!(f, x, y, t)
        f[1] = -x[2]-(x[1])^2
        f[2] = 0.0001 
        return nothing
end

# Define the initial state and final time
x0 = [+1.00,-1.00]
T = 10000.0
δt = 100e0

# Evolve the dynamical system from different initial conditions 
bfold = ContinuousDynamicalSystem(iip_fs2!, x0, nothing)
Xt, t = trajectory(bfold, T; Δt=δt)
X2 = Xt[:,1]
Y2 = Xt[:,2]

# SINGULAR LIMIT - Trajectory 3

# Define the dynamics
function iip_fs3!(f, x, y, t)
        f[1] = -x[2]-(x[1])^2
        f[2] = 0.0 
        return nothing
end

# Define the initial state and final time
x0 = [-0.01,0.0]
T = 100.0
δt = 1e-1

# Evolve the dynamical system from different initial conditions 
bfold = ContinuousDynamicalSystem(iip_fs3!, x0, nothing)
Xt, t = trajectory(bfold, T; Δt=δt)
X3 = Xt[:,1]
Y3 = Xt[:,2]

# Plot singular limit trajectory 1 
scatter!(ax, Y1, X1, color = :blue, strokecolor = :black, strokewidth = 1.5, markersize = 9)
# Plot singular limit trajectory 2 
scatter!(ax, Y2, X2, color = :red, strokecolor = :black, strokewidth = 1.5, markersize = 9)
# Plot singular limit trajectory 3 
scatter!(ax, Y3, X3, color = :blue, strokecolor = :black, strokewidth = 1.5, markersize = 9)
# Plot the fast-slow trajectory 
lines!(ax, Y0, X0, color = :orange, linewidth = 2.5)
# Export the results
save("../../../results/fast_slow/fig2.1.1.png", fig1)
