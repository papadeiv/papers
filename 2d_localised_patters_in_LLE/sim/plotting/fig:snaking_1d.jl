include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")

using MAT

###################################
#     Bifurcation diagram (1)     #
###################################

# Import data from mat files
file = matread("../data/cont_rho/branch.mat")
key = file["branch"]
λ = key[:,2]
ρ = key[:,3]
u_2 = key[:,5]

# Define number of columns for the functions subplots
N_plt = 4::Int64

# Define indices for the 2-dimensional continuation 
ν_cont_plt = [21,321]
ν_idx_plt = ρ[ν_cont_plt]

# Extract the indices of the saddle-node bifurcations 
idx_bif = findall(diff(λ) .!=0)

# Check whether there have been bifurcations
if length(idx_bif) > 0
        # Separate the stable solutions from the unstable ones
        global idx_stb = Vector{Float64}(undef, length(idx_bif))

        # Update the stability array
        for n in 1:length(idx_bif)
                if λ[idx_bif[n]-1] > 0
                        idx_stb[n] = +1
                else
                        idx_stb[n] = -1
                end
        end

        # Update the bifurcation and stability indices to include the endpoints
        idx_bif = [1; idx_bif; length(λ)]
        idx_stb = [idx_stb[2]; idx_stb; idx_stb[end]]

        # Define a global variable for the number of bifurcations
        global N_bif = length(idx_bif)-1
else
        # Create a bifurcation index vector that spans the whole branch  
        global idx_bif = [1, length(λ)]
        
        # Create a single-entry stability vector
        if λ[1] > 0
                global idx_stb = [+1]
        else
                global idx_stb = [-1]
        end

        # Define a global variable for the number of bifurcations
        global N_bif = 1
end

# Create and customise the figure 
fig, ax = mkfig(size = [1200,1200],
                bg_out = :white,
                box_position = [1:3,1:4],
                pad = (15,40,25,25), # Order is: left, right, bottom, top 
                limits = ((1.111, 1.118), (0.0025, 0.335)),
                lab = [L"\mathbf{\rho}", L"\mathbf{||u - u_0||}_{2}"],
                lab_pad = [0.0,0.0],
                lab_size = [30,30],
                x_ticks = [1.111,1.11208,1.11467,1.118],
                y_ticks = [minimum(u_2),maximum(u_2)],
                ticks_lab_size = [30,30],
                ticks_lab_trunc = [3,2],
               )
# Plot vertical lines for the folding points
lines!(ax, [1.11208,1.11208], [0.0025,0.335], linewidth = 1.0, color = :black)
lines!(ax, [1.11467,1.11467], [0.0025,0.335], linewidth = 1.0, color = :black)
# Loop over the saddle-node bifurcation indices
for n in 1:(length(idx_bif)-1)
        # Check the stability of the subbranch
        if idx_stb[n] < 0
                # Plot the branch of stable solutions
                lines!(ax, ρ[idx_bif[n]:idx_bif[n+1]], u_2[idx_bif[n]:idx_bif[n+1]], linewidth = 3.0, color = :magenta4)
        else
                # Plot the branch of unstable solutions
                lines!(ax, ρ[idx_bif[n]:idx_bif[n+1]], u_2[idx_bif[n]:idx_bif[n+1]], linewidth = 3.0, color = :magenta4, linestyle = :dash)
        end
end

#####################
#     Solutions     #
#####################

# Define variables and parameters
L = 80.0::Float64
x = collect(LinRange(-L,L,1000)) 
Nh = 500::Int64
u_re = Matrix{Float64}(undef, 2*Nh, 2*N_plt)
u_im = Matrix{Float64}(undef, 2*Nh, 2*N_plt)
ρ_plt = Vector{Float64}(undef, 2*N_plt)
idx_plt = [1,51,101,201,321,481,801,861] 

# Loop over the plotting indices
for n in 1:(2*N_plt)
        # Read in solution from formatted index
        filename = @sprintf("../data/cont_rho/solution_%07d.mat", idx_plt[n]-1)
        file = matread(filename)
        
        # Extract the solution fields
        u = file["u"]
        u_re[:,n] = [reverse(u[1:Nh]); u[1:Nh]]
        u_im[:,n] = [reverse(u[Nh+1:end]); u[Nh+1:end]]

        # Extract the continuation parameter value
        ρ = file["p"]
        ρ_plt[n] = ρ[2]
end

# Locate the solutions you're about to plot in the bifurcation diagram
scatter!(ax, ρ_plt, u_2[idx_plt], markersize = 30, strokewidth = 5, strokecolor = :black, color = :lightslateblue)

# Array of boolean variables for plotting labels and ticks
plot_bool = [true, false, false, false]

# First row of solutions
for n in 1:4
        local nullfig, ax = mkfig(fig = fig,
                        box_position = [4,n],
                        pad = (15,35,25,145), # Order is: left, right, bottom, top 
                        limits = ((-L,L), (-0.6,1.5)),
                        lab = [L"\mathbf{x}", L"\mathbf{u}"],
                        toggle_lab = [false,plot_bool[n]],
                        lab_pad = [0.0,0.0],
                        lab_size = [30,30],
                        x_ticks = [-L,0,L],
                        y_ticks = [-0.6,1.5],
                        ticks_lab_size = [30,30],
                        toggle_ticks_lab = [false,plot_bool[n]],
                        ticks_lab_trunc = [0,1],
                       )
        # Plot the real field
        lines!(ax, x, u_re[:,n], linewidth = 3.0, color = :blue)
        # Plot the imaginary field
        lines!(ax, x, u_im[:,n], linewidth = 3.0, color = :red)
end

# Second row of solutions
for n in 1:4
        local nullfig, ax = mkfig(fig = fig,
                        box_position = [5,n],
                        pad = (15,35,25,145), # Order is: left, right, bottom, top 
                        limits = ((-L,L), (-0.6,1.5)),
                        lab = [L"\mathbf{x}", L"\mathbf{u}"],
                        toggle_lab = [true,plot_bool[n]],
                        lab_pad = [0.0,0.0],
                        lab_size = [30,30],
                        x_ticks = [-L,0,L],
                        y_ticks = [-0.6,1.5],
                        ticks_lab_size = [30,30],
                        toggle_ticks_lab = [true,plot_bool[n]],
                        ticks_lab_trunc = [0,1],
                       )
        # Plot the real field
        lines!(ax, x, u_re[:,n+4], linewidth = 3.0, color = :blue)
        # Plot the imaginary field
        lines!(ax, x, u_im[:,n+4], linewidth = 3.0, color = :red)
end

# Adjust whitespace between rows
rowgap!(fig.layout, 50)

# Export the figure 
save("../fig/fig:snaking_1d.png", fig)

# Locate the solutions you continued in nu (2-dimensional case) 
scatter!(ax, ν_idx_plt, u_2[ν_cont_plt], marker = :star5, markersize = 40, strokewidth = 3, strokecolor = :black, color = :yellow)

# Export the figure 
save("../fig/fig:snaking_1d_to_2d.png", fig)
