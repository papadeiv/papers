include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")
include("../../../../inc/SystemAnalysis.jl")

# Import the data from csv
u = readin("../data/slide_09/u.csv")
μ = readin("../data/slide_09/μ.csv")
ut = readin("../data/slide_09/ut.csv")
eq1 = readin("../data/slide_09/stable_eq_1.csv")
eq2 = readin("../data/slide_09/stable_eq_2.csv")
eq3 = readin("../data/slide_09/unstable_eq.csv")

# Create and customise the augmented timeseries figure
fig, ax = mkfig(size = [1200,800],
                bg_out = "#eeeeeeff",
                limits = ((-3,4), (-1,4)),
                lab = [L"\mathbf{\mu}", L"\mathbf{x}"],
                lab_pad = [-60.0,-60.0],
                x_ticks = [-3,4],
                y_ticks = [-1,4],
                ticks_lab_trunc = [0,0]
               )
# Plot the equilibria in the augmented phase space (i.e. bifurcation diagram)
lines!(ax, eq1[:,1], eq1[:,2], color = :black, linewidth = 4)
lines!(ax, eq2[:,1], eq2[:,2], color = :black, linewidth = 4)
lines!(ax, eq3[:,1], eq3[:,2], color = :black, linewidth = 4, linestyle = :dash)

# Plot the augmented timeseries
lines!(ax, ut[:,1], ut[:,2], color = (:red, 0.25), linewidth = 1)

# Export the augmented timeseries plot 
save("../fig/slide_09/timeseries.png", fig)

# Define the scalar potential function 
V(x, p) = p*x + x^2 - x^3 + (1/5)*(x^4)

# Define the limits for the scalar potential plot
x_inf = -(2.0::Float64)
x_sup = 4.0::Float64
y_inf = -(1.0::Float64)
y_sup = 4.0::Float64 

# Define the domain of the function
domain = LinRange(x_inf, x_sup, 1000)

# Define index for the potential snapshot
idx = 39500

# Create and customise the scalar potential function figure
fig, ax = mkfig(size = [1000,1000],
                bg_out = "#eeeeeeff",
                limits = ((x_inf,x_sup), (y_inf,y_sup)),
                lab = [L"\mathbf{x}", L"\mathbf{V(x)}"],
                lab_pad = [-60.0,-60.0],
                x_ticks = [x_inf,x_sup],
                y_ticks = [y_inf,y_sup],
                ticks_lab_trunc = [0,0]
               )
# Plot the scalar potential (landscape)
lines!(ax, domain, [V(x,μ[idx]) for x in domain], color = :black, linewidth = 4.5)
# Plot the state (ball in the landscape)
scatter!(ax, u[idx], V(u[idx],μ[idx]), color = (:red, 1.0), markersize = 45, strokewidth = 3)

# Export the scalar potential function plot 
save("../fig/slide_09/potential.png", fig)
