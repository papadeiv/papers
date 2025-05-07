include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")
include("../../../../inc/SystemAnalysis.jl")

# Import the data from csv
u = readin("../data/figure_02/u.csv")
μ = readin("../data/figure_02/μ.csv")
ut = readin("../data/figure_02/ut.csv")
eq1 = readin("../data/figure_02/stable_eq_1.csv")
eq2 = readin("../data/figure_02/stable_eq_2.csv")
eq3 = readin("../data/figure_02/unstable_eq.csv")

# Create and customise the augmented timeseries figure
fig, ax = mkfig(size = [1200,600],
                pad = (30,60,30,30), # Order is: left, right, bottom, top 
                limits = ((-3,4), (-1,4)),
                lab = [L"\mathbf{\mu}", L"\mathbf{x}"],
                lab_size = [40,40],
                lab_pad = [-50.0,-50.0],
                ticks_lab_size = [40,40],
                x_ticks = [-3,4],
                y_ticks = [-1,4],
                ticks_lab_trunc = [0,0]
               )
# Plot the equilibria in the augmented phase space (i.e. bifurcation diagram)
lines!(ax, eq1[:,1], eq1[:,2], color = :darkgoldenrod2, linewidth = 4)
lines!(ax, eq2[:,1], eq2[:,2], color = :darkgoldenrod2, linewidth = 4)
lines!(ax, eq3[:,1], eq3[:,2], color = :darkgoldenrod2, linewidth = 4, linestyle = :dash)

# Plot the augmented timeseries
#lines!(ax, ut[:,1], ut[:,2], color = (:teal, 0.5), linewidth = 1)
lines!(ax, μ, u, color = (:teal, 0.5), linewidth = 1)

# Export the augmented timeseries plot 
save("../fig/fig2.png", fig)

# Define the scalar potential function 
V(x, p) = p*x + x^2 - x^3 + (1/5)*(x^4)

# Define the parameter value at which to plot the potential function
idx1 = 22000
μ1 = trunc(μ[idx1], digits=3)
idx2 = 36001
μ2 = trunc(μ[idx2], digits=3)
idx3 = 50001
μ3 = trunc(μ[idx3], digits=3)

########################
#       μ1=-0.80       #
########################

# Define the limits for the scalar potential plot
x_inf = -(2.0::Float64)
x_sup = 5.0::Float64
y_inf = -(6.0::Float64)
y_sup = 2.0::Float64 

# Define the domain of the function
domain = LinRange(x_inf, x_sup, 1000)

# Create and customise the scalar potential function figure
fig, ax = mkfig(size = [1000,1000],
                border = 15.0,
                limits = ((x_inf,x_sup), (y_inf,y_sup)),
                title = L"\mathbf{\mu=}\text{%$(μ1)}",
                toggle_title = true,
                title_size = 120,
                title_gap = 14.0,
                toggle_lab = [false,false],
                toggle_ticks = [false,false],
                toggle_ticks_lab = [false,false]
               )
# Plot the scalar potential (landscape)
lines!(ax, domain, [V(x,μ1) for x in domain], color = :darkgoldenrod2, linewidth = 10)
# Plot the state (ball in the landscape)
scatter!(ax, u[idx1], V(u[idx1],μ1), color = (:teal, 1.0), markersize = 75, strokewidth = 3)

# Export the scalar potential function plot 
save("../fig/fig2.1.png", fig)

########################
#       μ2=0.600       #
########################

# Define the limits for the scalar potential plot
x_inf = -(2.0::Float64)
x_sup = 4.0::Float64
y_inf = -(2.0::Float64)
y_sup = 5.0::Float64 

# Define the domain of the function
domain = LinRange(x_inf, x_sup, 1000)

# Create and customise the scalar potential function figure
fig, ax = mkfig(size = [1000,1000],
                border = 15.0,
                limits = ((x_inf,x_sup), (y_inf,y_sup)),
                title = L"\mathbf{\mu=}\text{%$(μ2)}",
                toggle_title = true,
                title_size = 120,
                title_gap = 14.0,
                toggle_lab = [false,false],
                toggle_ticks = [false,false],
                toggle_ticks_lab = [false,false]
               )
# Plot the scalar potential (landscape)
lines!(ax, domain, [V(x,μ2) for x in domain], color = :darkgoldenrod2, linewidth = 10)
# Plot the state (ball in the landscape)
scatter!(ax, u[idx2], V(u[idx2],μ2), color = (:teal, 1.0), markersize = 75, strokewidth = 3)

# Export the scalar potential function plot 
save("../fig/fig2.2.png", fig)

########################
#       μ3=2.000       #
########################

# Define the limits for the scalar potential plot
x_inf = -(2.0::Float64)
x_sup = 4.0::Float64
y_inf = -(2.0::Float64)
y_sup = 4.5::Float64 

# Define the domain of the function
domain = LinRange(x_inf, x_sup, 1000)

# Create and customise the scalar potential function figure
fig, ax = mkfig(size = [1000,1000],
                border = 15.0,
                limits = ((x_inf,x_sup), (y_inf,y_sup)),
                title = L"\mathbf{\mu=}\text{%$(μ3)}",
                toggle_title = true,
                title_size = 120,
                title_gap = 14.0,
                toggle_lab = [false,false],
                toggle_ticks = [false,false],
                toggle_ticks_lab = [false,false]
               )
# Plot the scalar potential (landscape)
lines!(ax, domain, [V(x,μ3) for x in domain], color = :darkgoldenrod2, linewidth = 10)
# Plot the state (ball in the landscape)
scatter!(ax, u[idx3], V(u[idx3],μ3), color = (:teal, 1.0), markersize = 75, strokewidth = 3)

# Export the scalar potential function plot 
save("../fig/fig2.3.png", fig)
