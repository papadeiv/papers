include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")
include("../../../../inc/SystemAnalysis.jl")

# Import the data from csv
time = readin("../data/slide_02/time.csv")
μ = readin("../data/slide_02/μ.csv")
x1 = readin("../data/slide_02/x1.csv")
x2 = readin("../data/slide_02/x2.csv")

# Define the limits for the phase space plot
x1_inf = -(1.0::Float64)
x1_sup = 3.5::Float64
x2_inf = -(1.5::Float64)
x2_sup = 1.5::Float64 

# Define a 2-dimensional bistable saddle-node normal form
f(x, p) = -2*p - 2*x + 3*(x^2) - (4/5)*(x^3)

# Loop over the parameter values
printstyled("Generating figures\n"; bold=true, underline=true, color=:light_blue)
@showprogress for n in 1:length(μ)
        # Get the equilibria of the scalar dynamics for x
        eq_x = get_equilibria(f, μ[n], domain=[-1,4])
        stable_eq_x = eq_x[1]
        unstable_eq_x = eq_x[2]
        # Get the equilibria of the scalar dynamics for y
        eq_y = -μ[n]

        # Create and customise the phase space figure
        fig, ax = mkfig(size = [1200,1200],
                    lab = [L"\mathbf{x_1}", L"\mathbf{x_2}"],
                    lab_pad = [-60.0,-60.0],
                    limits = ((x1_inf, x1_sup), (x2_inf, x2_sup)),
                    x_ticks = [x1_inf, x1_sup],
                    y_ticks = [x2_inf, x2_sup],
                    ticks_lab_trunc = [1,1]
                   )
        # Plot the vector field
        x1_axis = LinRange(x1_inf, x1_sup, 15)
        x2_axis = LinRange(x2_inf, x2_sup, 15)
        f1_field = [y-μ[n]-2*x+3*x^2-(4/5)*x^3 for x in x1_axis, y in x2_axis] 
        f2_field = [-y-μ[n] for x in x1_axis, y in x2_axis]
        arrows!(ax, x1_axis, x2_axis, f1_field, f2_field, arrowsize = 10, lengthscale = 0.035, arrowcolor = :black, linecolor = :black)
        # Plot the trajectories
        for j in 1:(n-1)
                lines!(ax, x1[:,j], x2[:,j], color = (:black,0.25), linewidth = 3)
        end
        lines!(ax, x1[:,n], x2[:,n], color = :black, linewidth = 3)
        # Plot the equilibria
        if n == 1
                scatter!(ax, unstable_eq_x[1], eq_y, color = :red, markersize = 30, strokewidth = 3)
                scatter!(ax, stable_eq_x[1], eq_y, color = :blue, markersize = 30, strokewidth = 3)
                scatter!(ax, stable_eq_x[2], eq_y, color = :blue, markersize = 30, strokewidth = 3)
        elseif n < length(μ) 
                scatter!(ax, stable_eq_x[1], eq_y, color = :blue, markersize = 30, strokewidth = 3)
                scatter!(ax, unstable_eq_x[1], eq_y, color = :red, markersize = 30, strokewidth = 3)
                scatter!(ax, stable_eq_x[2], eq_y, color = :blue, markersize = 30, strokewidth = 3)
        else
                scatter!(ax, stable_eq_x[1], eq_y, color = :blue, markersize = 30, strokewidth = 3)
        end
        # Plot the IC
        scatter!(ax, x1[1,1], x2[1,1], color = :yellow, markersize = 30, strokewidth = 3)
 
        # Export the phase plane figure 
        save("../fig/slide_02/phase_portrait/$n.png", fig)

        # Create and customise the timeseries figure
        fig, ax = mkfig(size = [1200,600],
                    lab = [L"\textbf{time}", L"\mathbf{x_1}"],
                    lab_pad = [-60.0,-60.0],
                    limits = ((0, time[end]), (-1, 3)),
                    x_ticks = [0, time[end]],
                    y_ticks = [-1, 3],
                    ticks_lab_trunc = [0,0]
                   )
        # Plot the trajectories
        for j in 1:(n-1)
                lines!(ax, LinRange(0.0, time[end], length(time)), x1[:,j], color = (:black,0.25), linewidth = 3)
        end
        lines!(ax, LinRange(0.0, time[end], length(time)), x1[:,n], color = :black, linewidth = 3)
        # Plot the equilibria
        if n == 1
                lines!(ax, [0.0,time[end]], [unstable_eq_x[1],unstable_eq_x[1]], color = (:red,0.25), linewidth = 3)
                lines!(ax, [0.0,time[end]], [stable_eq_x[1],stable_eq_x[1]], color = (:blue,0.25), linewidth = 3)
                lines!(ax, [0.0,time[end]], [stable_eq_x[2],stable_eq_x[2]], color = (:blue,0.25), linewidth = 3)
        elseif n < length(μ) 
                lines!(ax, [0.0,time[end]], [unstable_eq_x[1],unstable_eq_x[1]], color = (:red,0.25), linewidth = 3)
                lines!(ax, [0.0,time[end]], [stable_eq_x[1],stable_eq_x[1]], color = (:blue,0.25), linewidth = 3)
                lines!(ax, [0.0,time[end]], [stable_eq_x[2],stable_eq_x[2]], color = (:blue,0.25), linewidth = 3)
        else
                lines!(ax, [0.0,time[end]], [stable_eq_x[1],stable_eq_x[1]], color = (:blue,0.25), linewidth = 3)
        end
        # Plot the IC
        scatter!(ax, 0.0, x1[1,1], color = :yellow, markersize = 30, strokewidth = 3)

        # Export the timeseries figure 
        save("../fig/slide_02/timeseries/$n.png", fig)
end
