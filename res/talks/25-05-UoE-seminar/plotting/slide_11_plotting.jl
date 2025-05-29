include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")
include("../../../../inc/EscapeProblem.jl")
include("../../../../inc/PotentialLearning.jl")
include("../../../../inc/TimeseriesAnalysis.jl")

# Define the stochastic diffusion
σ = 0.250::Float64
D = (σ^2)/2.0::Float64

# Compute a shift for the potential {c0} that sets V(xs)=0 to avoid numerical cancellation
xs(μ) = (1/(3*μ[3]))*(sqrt((μ[2])^2 - 3*μ[1]*μ[3]) - μ[2])
c0(μ) = - μ[1]*xs(μ) - μ[2]*(xs(μ))^2 - μ[3]*(xs(μ))^3

# Define an arbitrary cubic with the the above constraint on {c0}
V(x, μ) = c0(μ) + μ[1]*x + μ[2]*(x^2) + μ[3]*(x^3)
 
# Define the exact stationary probability distribution
f(x, μ) = exp(-(1.0::Float64/D)*(V(x, μ)))
N(μ) = get_normalisation_constant(f, (-(1/(3*μ[3]))*(sqrt((μ[2])^2 - 3*μ[1]*μ[3]) + μ[2]), Inf), parameters=μ)
p(x, μ) = N(μ)*f(x, μ)

# Import the data from csv
residuals = readin("../data/slide_10/residuals/1.csv")
Nw = length(residuals[:,1])

mean_linear_coefficients = readin("../data/slide_11/mean_linear_coefficients.csv")
mean_nonlinear_coefficients = readin("../data/slide_11/mean_nonlinear_coefficients.csv")
lin_linear_coefficients = readin("../data/slide_11/lin_linear_coefficients.csv")
lin_nonlinear_coefficients = readin("../data/slide_11/lin_nonlinear_coefficients.csv")
Ns = length(lin_nonlinear_coefficients[:,1])

solution = readin("../data/slides_10_to_13/solution.csv")
t = solution[:,1]
μ = solution[:,2]
Nt = length(t)

# Define a list of plotting indices
idx = collect(StepRange(1, Int64(10), Ns)) 

# Loop over the window's strides 
printstyled("Generating the figures for the equilibrium density using $(Threads.nthreads()) threads\n"; bold=true, underline=true, color=:light_blue)
@showprogress Threads.@threads for m in 1:length(idx) 
        # Get the current plotting index
        n = idx[m] 
        
        # Define the index of the timestep at the end of the sliding window
        n_end = Nw + n - 1::Int64

        # Import the data
        mean_hist = readin("../data/slide_10/histograms/SN_mean/$n.csv")
        mean_bins = mean_hist[:,1]
        mean_pdf = mean_hist[:,2]
        lin_hist = readin("../data/slide_10/histograms/SN_lin/$n.csv")
        lin_bins = lin_hist[:,1]
        lin_pdf = lin_hist[:,2]

        # Define domain limits for plotting purposes
        mean_bins_inf = mean_bins[1]
        mean_bins_sup = mean_bins[end]
        mean_domain = collect(LinRange(mean_bins_inf, mean_bins_sup, 1000))
        lin_bins_inf = lin_bins[1]
        lin_bins_sup = lin_bins[end]
        lin_domain = collect(LinRange(lin_bins_inf, lin_bins_sup, 1000))

        # Create and customise the mean detrended histogram figure 
        fig, ax = mkfig(size = [1200,800],
                        pad = (30,70,10,30), # Order is: left, right, bottom, top 
                        bg_out = colorant"rgb(147,237,213)",
                        box_position = [2,1],
                        border_color = RGBAf(0,0.236,0.236),
                        limits = ((mean_bins_inf,mean_bins_sup), (-0.1,6)),
                        lab = [L"\mathbf{x}", L"\mathbf{\rho}"],
                        lab_size = [40,40],
                        lab_color = [RGBAf(0,0.236,0.236), RGBAf(0,0.236,0.236)],
                        lab_pad = [-50.0,-40.0],
                        x_ticks = [mean_bins_inf,mean_bins_sup],
                        y_ticks = [0,6],
                        ticks_color = [RGBAf(0,0.236,0.236), RGBAf(0,0.236,0.236)],
                        ticks_lab_size = [40,40],
                        ticks_lab_xpos = [:right,:top],
                        ticks_lab_trunc = [1,0]
                       )
        # Plot the mean detrended histogram 
        lines!(ax,mean_bins, mean_pdf, linewidth = 5, color = (colorant"rgb(207,19,92)", 0.3))
        # Plot the least-squares solution for the mean detrended histogram
        lines!(ax, mean_domain, [p(x, mean_nonlinear_coefficients[n,:]) for x in mean_domain], linewidth = 3.5, color = colorant"rgb(207,19,92)")

        # Create and customise the linear detrended histogram figure 
        fig, ax = mkfig(fig = fig,
                        bg_out = colorant"rgb(147,237,213)",
                        box_position = [2,2],
                        border_color = RGBAf(0,0.236,0.236),
                        limits = ((lin_bins_inf,lin_bins_sup), (-0.1,6)),
                        toggle_lab = [true,false],
                        lab = [L"\mathbf{x}", L"\mathbf{\rho}"],
                        lab_size = [40,40],
                        lab_color = [RGBAf(0,0.236,0.236), RGBAf(0,0.236,0.236)],
                        lab_pad = [-50.0,-40.0],
                        x_ticks = [lin_bins_inf,lin_bins_sup],
                        y_ticks = [0,6],
                        ticks_color = [RGBAf(0,0.236,0.236), RGBAf(0,0.236,0.236)],
                        ticks_lab_size = [40,40],
                        ticks_lab_xpos = [:left,:top],
                        toggle_ticks_lab = [true,false],
                        ticks_lab_trunc = [1,0]
                       )
        # Plot the linear detrended histogram 
        lines!(ax,lin_bins, lin_pdf, linewidth = 5, color = (colorant"rgb(38,106,176)", 0.3))
        # Plot the least-squares solution for the linear detrended histogram
        lines!(ax, lin_domain, [p(x, lin_nonlinear_coefficients[n,:]) for x in lin_domain], linewidth = 3.5, color = colorant"rgb(38,106,176)")

        # Create and customise the least-squares solution coefficients timeseries figure 
        fig, ax = mkfig(fig = fig,
                        bg_out = colorant"rgb(147,237,213)",
                        box_position = [1,1:2],
                        border_color = RGBAf(0,0.236,0.236),
                        limits = ((μ[1],μ[end]), (0,4)),
                        lab = [L"\mathbf{\mu}", L"\mathbf{c_3}"],
                        lab_size = [40,40],
                        lab_color = [RGBAf(0,0.236,0.236), RGBAf(0,0.236,0.236)],
                        lab_pad = [-50.0,-40.0],
                        x_ticks = [μ[1],μ[Nw],μ[end]],
                        y_ticks = [0,4],
                        ticks_lab_size = [40,40],
                        ticks_color = [RGBAf(0,0.236,0.236), RGBAf(0,0.236,0.236)],
                        ticks_lab_trunc = [1,0]
                       )
        # Plot the solution coefficients from the mean detrended timeseries
        lines!(ax, μ[Nw:n_end], mean_nonlinear_coefficients[1:n,3], linewidth = 3.5, color = colorant"rgb(207,19,92)")
        scatter!(ax, μ[n_end], mean_nonlinear_coefficients[n,3], color = (colorant"rgb(207,19,92)",0.3), markersize = 30, strokewidth = 3)
        # Plot the solution coefficients from the linear detrended timeseries
        lines!(ax, μ[Nw:n_end], lin_nonlinear_coefficients[1:n,3], linewidth = 3.5, color = colorant"rgb(38,106,176)")
        scatter!(ax, μ[n_end], lin_nonlinear_coefficients[n,3], color = (colorant"rgb(38,106,176)",0.3), markersize = 30, strokewidth = 3)
 
        # Export the figure
        save("../fig/slide_11/$m.png", fig)
end
