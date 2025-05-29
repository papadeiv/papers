include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")
include("../../../../inc/EscapeProblem.jl")
include("../../../../inc/PotentialLearning.jl")

# Import the data from csv
residuals = readin("../data/slide_10/residuals/1.csv")
Nw = length(residuals[:,1])

mean_nonlinear_coefficients = readin("../data/slide_11/mean_nonlinear_coefficients.csv")
lin_nonlinear_coefficients = readin("../data/slide_11/lin_nonlinear_coefficients.csv")

shifts = readin("../data/slide_12/shifts.csv")
xs_mean = shifts[:,1] 
xs_lin = shifts[:,2]
ys_mean = shifts[:,3]
ys_lin = shifts[:,4]

errors = readin("../data/slide_12/errors.csv")
μ = errors[:,1] 
x_uns = errors[:,2] 
x_stb = errors[:,3] 
mean_error = errors[:,4] 
lin_error = errors[:,5] 
Ns = length(μ)

# Define a list of plotting indices
idx = collect(StepRange(1, Int64(10), Ns)) 

# Loop over the window's strides 
printstyled("Generating the figures for the scalar potential using $(Threads.nthreads()) threads\n"; bold=true, underline=true, color=:light_blue)
@showprogress Threads.@threads for m in 1:length(idx)
         # Get the current plotting index
         n = idx[m] 
 
         # Define the index of the timestep at the end of the sliding window
         n_end = Nw + n - 1::Int64

         # Extract the current solutions
         c_mean = mean_nonlinear_coefficients[n,:]
         c_lin = lin_nonlinear_coefficients[n,:]

         # Define domain limits for plotting purposes
         x_inf = -1.5 
         x_sup = 4.0
         domain = collect(LinRange(x_inf, x_sup, 2000))
         y_inf = -1.0
         y_sup = 3.0 

         # Define the analytical potential function
         U(x) = μ[n]*x + x^2 - x^3 + (1/5)*x^4

         # Define the shifted reconstructions
         Vs_mean(x) = ys_mean[n] + c_mean[1]*(x-xs_mean[n]) + c_mean[2]*(x-xs_mean[n])^2 + c_mean[3]*(x-xs_mean[n])^3
         Vs_lin(x) = ys_lin[n] + c_lin[1]*(x-xs_lin[n]) + c_lin[2]*(x-xs_lin[n])^2 + c_lin[3]*(x-xs_lin[n])^3

        # Create and customise the scalar potential figure 
        fig, ax = mkfig(size = [1200,1200],
                        bg_out = colorant"rgb(147,237,213)",
                        box_position = [2:3,1],
                        border_color = RGBAf(0,0.236,0.236),
                        limits = ((x_inf,x_sup), (y_inf,y_sup)),
                        lab = [L"\mathbf{x}", L"\mathbf{V(x)}"],
                        lab_size = [40,40],
                        lab_color = [RGBAf(0,0.236,0.236), RGBAf(0,0.236,0.236)],
                        lab_pad = [-50.0,-20.0],
                        x_ticks = [x_inf,x_sup],
                        y_ticks = [y_inf,y_sup],
                        ticks_color = [RGBAf(0,0.236,0.236), RGBAf(0,0.236,0.236)],
                        ticks_lab_size = [40,40],
                        ticks_lab_trunc = [0,0]
                       )
        # Plot the analytic potential 
        lines!(ax, domain, [U(x) for x in domain], linewidth = 5, color = :black)
        # Plot the reconstruction from the linear detrended timeseries 
        lines!(ax, domain, [Vs_lin(x) for x in domain], linewidth = 6.5, color = (colorant"rgb(38,106,176)",1.00))
        # Plot the reconstruction from the mean detrended timeseries
        lines!(ax, domain, [Vs_mean(x) for x in domain], linewidth = 2.5, color = (colorant"rgb(207,19,92)",0.9))

        # Create and customise the error decay figure 
        fig, ax = mkfig(fig = fig,
                        bg_out = colorant"rgb(147,237,213)",
                        box_position = [1,1],
                        border_color = RGBAf(0,0.236,0.236),
                        limits = ((μ[1],μ[end]), (0,160)),
                        lab = [L"\mathbf{\mu}", L"\mathbf{||V-V_{\textbf{loc}}||_2}"],
                        lab_size = [40,40],
                        lab_color = [RGBAf(0,0.236,0.236), RGBAf(0,0.236,0.236)],
                        lab_pad = [-50.0,-60.0],
                        x_ticks = [μ[1],μ[Nw],μ[end]],
                        y_ticks = [0,160],
                        ticks_color = [RGBAf(0,0.236,0.236), RGBAf(0,0.236,0.236)],
                        ticks_lab_size = [40,40],
                        ticks_lab_trunc = [1,0]
                       )
        # Plot the error from the mean detrended timeseries
        lines!(ax, μ[1:n], mean_error[1:n], linewidth = 3.5, color = colorant"rgb(207,19,92)")
        scatter!(ax, μ[n], mean_error[n], color = (colorant"rgb(207,19,92)",0.3), markersize = 30, strokewidth = 3)
        # Plot the error from the linear detrended timeseries
        lines!(ax, μ[1:n], lin_error[1:n], linewidth = 3.5, color = colorant"rgb(38,106,176)")
        scatter!(ax, μ[n], lin_error[n], color = (colorant"rgb(38,106,176)",0.3), markersize = 30, strokewidth = 3)
 
        # Export the figure
        save("../fig/slide_12/$m.png", fig)
end
