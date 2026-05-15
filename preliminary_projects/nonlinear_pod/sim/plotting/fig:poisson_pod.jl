include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")

# Define analytical solution 
function u(x)
        mu = 0.73
        return sin(pi*x)/(pi^2) + (mu*cos(2*pi*x))/(4*(pi^2)) - mu/(4*(pi^2))
end

# Domain and FEM discretization
const N = 200                      # Number of grid points
const x = LinRange(-1, 1, N)
const h = step(x)
μ_train = range(0.1, 1.0, length=20)

# Import data from csv
X = readin("../data/snapshots.csv")
S = readin("../data/singular_values.csv")
V = readin("../data/singular_vector.csv")
Vr = readin("../data/pod_basis.csv")
coefficients = readin("../data/coefficients.csv")
solutions = readin("../data/solutions.csv")
u_fom = solutions[:,1]
u_rom = solutions[:,2]

#################
#   Snapshots   #
#################

# Create and customise the figure 
fig, ax = mkfig(size = [1200,800],
                box_position = [1:2,1],
                bg_out = :white,
                pad = (35,30,30,35), # Order is: left, right, bottom, top 
                limits = ((x[1],x[end]), (-0.15,0.15)),
                lab = [L"\mathbf{x}", L"\mathbf{u}_{h}"],
                lab_pad = [-40.0,-60.0],
                lab_size = [30,30],
                x_ticks = [x[1],x[end]],
                y_ticks = [-0.15,0.15],
                ticks_lab_size = [30,30],
                ticks_lab_trunc = [0,2],
               )
# Loop over the columns of the snapshot matrix
for j in 1:size(X,2)
        lines!(ax, x, X[:, j], linewidth = 5, color = j, colormap = (:thermal, 0.25), colorrange = (1,size(X,2)))
end
# Plot the colorbar
Colorbar(fig[1:2,2],
         size = 45, 
         spinewidth = 5.0,
         limits = (μ_train[1], μ_train[end]), 
         colormap = (:thermal, 0.5),
         ticks = [μ_train[1], μ_train[end]],
         ticksize = 22,
         tickwidth = 5.0,
         tickformat = "{:.$(1)f}",
         ticklabelsize = 30,
         label = L"\mathbf{\mu}",
         labelpadding = -40.0,
         labelsize = 30
        )

# Export the figure
save("../fig/fig:poisson_fom.png", fig)

#############################
#   Singular values decay   #
#############################

# Create and customise the figure 
fig, ax = mkfig(size = [1200,800],
                bg_out = :white,
                pad = (30,60,10,35), # Order is: left, right, bottom, top 
                limits = ((0.75,5.25), (-0.35,5)),
                lab = [L"\text{index}(j)", L"\mathbf{\sigma}_{j}"],
                lab_pad = [-60.0,-40.0],
                x_ticks = [1,5],
                y_ticks = [0,5],
                ticks_lab_trunc = [0,0],
               )
# Plot the singular values decay
scatter!(ax, 1:length(S), S, markersize = 35, color = :red, strokewidth = 4, strokecolor = :black)
scatter!(ax, 1:2, S[1:2], markersize = 35, color = :green, strokewidth = 4, strokecolor = :black)

# Export the figure
save("../fig/fig:poisson_decay.png", fig)

##############################
#   Reduced order solution   #
##############################

# Create and customise the comparison figure 
fig, ax = mkfig(size = [2400,1600],
                box_position = [1:2,1:2],
                bg_out = :white,
                pad = (30,60,10,35), # Order is: left, right, bottom, top 
                limits = ((x[1],x[end]), (-0.15,0.15)),
                lab = [L"\mathbf{x}", L"\mathbf{u}_{h}"],
                lab_pad = [-60.0,-60.0],
                x_ticks = [x[1],x[end]],
                y_ticks = [-0.15,0.15],
                ticks_lab_trunc = [0,2],
               )
# Plot the FOM solution
lines!(ax, x, u_fom, linewidth = 7, color = (:red,1))
# Plot the ROM solution
lines!(ax, x, u_rom, linewidth = 10, linestyle = :dash, color = (:blue,1))

# Loop over the POD-basis vectors
for n in 1:2
        # Extract the coefficient of the ROM
        u_bar = trunc(coefficients[n], digits=2)
        # Create and customise the POD-basis figure 
        nullfig, ax = mkfig(fig=fig,
                            box_position = [3,n],
                            limits = ((x[1],x[end]), (-0.15,0.15)),
                            title = L"\bar{u}_{%$(n)} = %$(u_bar)",
                            toggle_title = true,
                            title_size = 50,
                            title_color = :black,
                            title_gap = 4.0,
                            lab = [L"\mathbf{x}", L"\mathbf{v}_{%$(n)}"],
                            lab_pad = [-60.0,-60.0],
                            x_ticks = [x[1],x[end]],
                            y_ticks = [-0.15,0.15],
                            ticks_lab_trunc = [0,2],
                           )
        # Plot the basis function
        lines!(ax, x, Vr[:,n], linewidth = 5, color = (:red,1))
        # Plot the weighted basis function
        lines!(ax, x, u_bar.*Vr[:,n], linewidth = 7, linestyle = :dash, color = (:blue,1))
end

# Export the figure
save("../fig/fig:poisson_rom.png", fig)
