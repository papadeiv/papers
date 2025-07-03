include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")
include("../../../../inc/EarlyWarningSignals.jl")

# Import the data from csv
t1 = readin("../data/slide_03/t1.csv")
u1 = readin("../data/slide_03/u1.csv")

printstyled("Generating the figures\n"; bold=true, underline=true, color=:light_blue)

############################
# Deterministic autonomous #
############################

# Create and customise the first figure 
fig, ax = mkfig(size = [1800,600],
                border_color = RGBAf(0,0.236,0.236),
                limits = ((0,7), (-0.6,2.65)),
                lab = [L"\textbf{time}", L"\mathbf{x}"],
                lab_pad = [-50.0,-50.0],
                lab_color = [RGBAf(0,0.236,0.236),RGBAf(0,0.236,0.236)],
                x_ticks = [0,7],
                y_ticks = [-0.6,2.65],
                ticks_lab_trunc = [0,1]
               )
# Define vector of indices for plotting purposes
idx = collect(UnitRange(1::Int64, length(u1[1,:])))
randomised = collect(UnitRange(1::Int64, length(u1[1,:])))
StatsBase.direct_sample!(idx, randomised)
# Loop over the parameter values in the sweep
@showprogress for n in 1:length(u1[1,:])
        # Plot the solution at the current parameter values
        lines!(ax, t1, u1[:,n], color = randomised[n], colormap = [:teal, :turquoise], colorrange = (1,length(u1[1,:])), linewidth = 3.5)
end

# Export the first plot 
save("../fig/slide_03_2.1.png", fig)

################################
# Deterministic non-autonomous #
################################

# Import the data from csv
u2and3 = readin("../data/slide_03/u2&3.csv")
t = u2and3[:,1] 
μ = u2and3[:,2]
u2 = u2and3[:,3]
u3 = u2and3[:,4]

# Create and customise the second figure 
fig, ax = mkfig(size = [1800,600],
                border_color = RGBAf(0,0.236,0.236),
                limits = ((t[1],t[end]), (-0.6,2.65)),
                lab = [L"\textbf{time}", L"\mathbf{x}"],
                lab_pad = [-50.0,-50.0],
                lab_color = [RGBAf(0,0.236,0.236),RGBAf(0,0.236,0.236)],
                x_ticks = [t[1],t[end]],
                y_ticks = [-0.6,2.65],
                ticks_lab_trunc = [0,1]
               )
# Plot the non-stationary solution 
lines!(ax, t, u2, color = :teal, linewidth = 3.5)

# Export the second plot 
save("../fig/slide_03_2.2.png", fig)

#############################
# Stochastic non-autonomous #
#############################

# Create and customise the second figure 
fig, ax = mkfig(size = [1800,600],
                border_color = RGBAf(0,0.236,0.236),
                limits = ((t[1],t[end]), (-0.6,2.65)),
                lab = [L"\textbf{time}", L"\mathbf{x}"],
                lab_pad = [-50.0,-50.0],
                lab_color = [RGBAf(0,0.236,0.236),RGBAf(0,0.236,0.236)],
                x_ticks = [t[1],t[end]],
                y_ticks = [-0.6,2.65],
                ticks_lab_trunc = [0,1]
               )
# Plot the non-stationary solution 
lines!(ax, t, u3, color = :teal, linewidth = 3.5)

# Export the second plot 
save("../fig/slide_03_2.3.png", fig)
