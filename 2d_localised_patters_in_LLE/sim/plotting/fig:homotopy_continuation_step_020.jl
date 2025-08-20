include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")

using MAT

###############################
#     Bifurcation diagram     #
###############################

# Import data from mat files
file = matread("../data/cont_nu_step_020/branch.mat")
key = file["branch"]
λ = key[:,2]
ν = key[:,3]
u_2 = key[:,5]

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
                box_position = [1:2,1:3],
                pad = (15,40,25,25), # Order is: left, right, bottom, top 
                limits = ((-0.025, 1.025), (0.025, 0.0525)),
                lab = [L"\mathbf{\nu}", L"\mathbf{||u - u_0||}_{2}"],
                lab_pad = [0.0,0.0],
                lab_size = [30,30],
                x_ticks = [0,1],
                y_ticks = [minimum(u_2),maximum(u_2)],
                ticks_lab_size = [30,30],
                ticks_lab_trunc = [0,3],
               )
# Loop over the saddle-node bifurcation indices
for n in 1:N_bif
        # Check the stability of the subbranch
        if idx_stb[n] < 0
                # Plot the branch of stable solutions
                lines!(ax, ν[idx_bif[n]:idx_bif[n+1]], u_2[idx_bif[n]:idx_bif[n+1]], linewidth = 5.0, color = :magenta4)
        else
                # Plot the branch of unstable solutions
                lines!(ax, ν[idx_bif[n]:idx_bif[n+1]], u_2[idx_bif[n]:idx_bif[n+1]], linewidth = 5.0, color = :magenta4#=, linestyle = :dash=#)
        end
end
# Plot the location of the 2 solutions displayed below
scatter!(ax, ν[1], u_2[1], markersize = 30, strokewidth = 5.0, color = :lightslateblue)
scatter!(ax, ν[end], u_2[end], markersize = 30, strokewidth = 5.0, color = :lightslateblue)

###########################
#     1-dim solutions     #
###########################

# Define variables and parameters
L = 80.0::Float64
r = collect(LinRange(-L,L,1000)) 
Nh = 500::Int64

# Import the 1-dimensional solution from mat files
file = matread("../data/cont_nu_step_020/solution_0000000.mat")
u = file["u"]
u_re = [reverse(u[1:Nh]); u[1:Nh]]
u_im = [reverse(u[Nh+1:end]); u[Nh+1:end]]
ν = file["p"]
ν_plt = trunc(Int64, ν[3])

# Create and customise the figure 
fig, ax = mkfig(fig = fig,
                bg_out = :white,
                box_position = [3,1],
                pad = (15,35,25,145), # Order is: left, right, bottom, top 
                title = L"\mathbf{\nu=%$ν_plt}",
                toggle_title = true,
                title_size = 30,
                limits = ((-L,L), (-0.6,1.5)),
                lab = [L"\mathbf{r}", L"\mathbf{u}"],
                lab_pad = [0.0,30.0],
                lab_size = [30,30],
                x_ticks = [-L,0,L],
                y_ticks = [-0.6,1.5],
                ticks_lab_size = [30,30],
                ticks_lab_trunc = [0,1],
               )
# Plot the real field
lines!(ax, r, u_re, linewidth = 3.0, color = :blue)
# Plot the imaginary field
lines!(ax, r, u_im, linewidth = 3.0, color = :red)

###########################
#     2-dim solutions     #
###########################

# Import the 2-dimensional solution from mat files
file = matread("../data/cont_nu_step_020/solution_0001320.mat")
u = file["u"]
u_re = u[1:Nh]
ν = file["p"]
ν_plt = trunc(Int64, ν[3])

# Define a discretised polar grid
r = range(0, L, length=Nh) 
θ = range(0, 2*pi, length=Nh)

# Reconstruct the 2-dimensional solution in polar coordinates
u_re_2d = [u_re[m] for n in θ, m in 1:length(r)]

# Define the size of the polar subplot
size = 450

# Create and customise the figure 
ax = PolarAxis(fig[3,2], 
               title = L"\mathbf{\nu=%$ν_plt}",
               titlesize = 30,
               titlegap = -10.0,
               rlimits = (:origin, L),
               rgridvisible = false,
               rticklabelsvisible = false,
               thetagridvisible = false,
               thetaticklabelsvisible = false,
               spinewidth = 5.0,
               width = size,
               height = size,
              )
# Plot the contour of the 2-dimensional field
s = surface!(ax, 0..2pi, 0..L, u_re_2d, color = u_re_2d, shading = NoShading, colormap = :viridis)
rowsize!(fig.layout, 3, Auto(1)) 
tightlimits!(ax)
# Plot the colorbar
Colorbar(fig[3,3],
         s, 
         size = 45,
         spinewidth = 5.0,
         ticks = [minimum(u_re_2d),maximum(u_re_2d)],
         label = L"\textbf{Re}\mathbf{(u)}",
         labelsize = 30,
         labelpadding = -30.0,
         ticklabelsize = 30,
         ticksize = 22, 
         tickwidth = 5.0,
         tickformat = "{:.1f}" 
        )

# Adjust whitespace between rows
colgap!(fig.layout, 1, 70)

# Export the figure 
save("../fig/fig:homotopy_continuation_step_020.png", fig)
