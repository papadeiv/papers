include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")

# Import the timeseries data 
time = readin("../data/slide_08/time.csv")
μ = readin("../data/slide_08/K.csv")
u = readin("../data/slide_08/u.csv")
u_mean = readin("../data/slide_08/mean_field.csv")
ut_mean = readin("../data/slide_08/detrended_mean_field.csv")
data = readin("../data/slide_08/ews_variance.csv")
μ_var = data[:,1]
u_var = data[:,2]

# Get number of lattice sites
Ns = size(u)[1]
N = convert(Int64, sqrt(Ns))

# Get number of timesteps
Nt = size(u)[2]

##############
# Mean field #
##############

# Create and customise the mean field figure 
fig, ax = mkfig(size = [1200,600],
                bg_out = "#eeeeeeff",
                limits = ((μ[1], μ[end]), (0.5, 1)),
                lab = [L"\mathbf{\mu}", L"\textbf{mean field}"],
                lab_pad = [-60.0,-60.0],
                x_ticks = [μ[1], μ[end]],
                y_ticks = [0.5, 1],
               )
# Plot the mean field timeseries
lines!(ax, μ, u_mean, linewidth = 3.5, color = (:black, 1.0))

# Export the figure
save("../fig/slide_08/mean_field.png", fig)

########################
# Detrended mean field #
########################

# Create and customise the detrended mean field figure 
fig, ax = mkfig(size = [1200,600],
                bg_out = "#eeeeeeff",
                limits = ((μ[1], μ[end]), (-1, 2.75)),
                lab = [L"\mathbf{\mu}", L"\textbf{residuals}"],
                toggle_lab = [false, true],
                lab_pad = [-60.0,-60.0],
                x_ticks = [μ[1], μ[end]],
                y_ticks = [-1, 2.75],
                toggle_ticks_lab = [false, true],
               )
# Plot the detrended mean field timeseries
lines!(ax, μ, ut_mean.*1e3, linewidth = 1.5, color = (:black, 1.0))
#text!((1.1,1.15),text=L"\times 10^{-5}",fontsize=36)

# Export the figure
save("../fig/slide_08/residuals.png", fig)

#####################
# Temporal variance #
#####################

# Create and customise the temporal variance of the detrended mean field figure 
fig, ax = mkfig(size = [1200,600],
                bg_out = "#eeeeeeff",
                limits = ((μ[1], μ[end]), (0.0, 1.0)),
                lab = [L"\mathbf{\mu}", L"\textbf{variance}"],
                lab_pad = [-60.0,-60.0],
                x_ticks = [μ[1], μ[end]],
                y_ticks = [0.0, 1.0],
                ticks_lab_trunc = [1,1]
               )
# Plot the variance timeseries
lines!(ax, μ_var, u_var.*1e6, linewidth = 3.5, color = (:brown2, 1.0))
#text!((1.1,1.15),text=L"\times 10^{-5}",fontsize=36)

# Export the figure
save("../fig/slide_08/variance.png", fig)

#####################
# Lattice snapshots #
#####################

# Loop over the timesteps of the solution
printstyled("Generating snapshots of the lattice\n"; bold=true, underline=true, color=:light_blue)
@showprogress for n in 1:Nt
        # Get the solution at the current timestep
        solution = u[:,n]
        # Reshape the solution into a square snapshot
        solution = reshape(solution, (N, N))
        # Reshape the snapshot solution for the heatmap
        snapshot = Array{Float64}(undef, N, N)
        for m in 0:(N-1)
                snapshot[m+1,:] = solution[N-m,:]
        end
        # Create and customise the lattice snapshot figure 
        local fig, ax = mkfig(size = [1000,1000],
                pad = 30,
                bg_out = "#eeeeeeff",
                toggle_lab = [false,false],
                toggle_ticks = [false,false],
                toggle_ticks_lab = [false,false],
               )
        # Plot the snapshot as a heatmap field 
        heatmap!(ax, snapshot', colorrange = (0,1))

        # Export the figure
        save("../fig/slide_08/snapshots/$n.png", fig)
end
