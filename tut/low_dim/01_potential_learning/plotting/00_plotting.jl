include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")
include("../../../../inc/PotentialLearning.jl")


# Import the parameter range
μ = readin("../data/parameter.csv")

# Get the number of parameter's values
Nμ = length(μ)

# Loop over the parameter values
printstyled("Creating the plots\n"; bold=true, underline=true, color=:light_blue)
@showprogress for n in 1:Nμ
        # Read the settings data 
        settings = readin("../data/settings/$n.csv")
        local σ = settings[1]
        local a = settings[2]
        local b = settings[3]
        local I = b - a
        local ϴ = settings[4]
        local D = (σ^2)/2
        # Read the OUP data
        OUP_distribution = readin("../data/distribution/OUP$n.csv")
        bins = OUP_distribution[:,1]
        pdf = OUP_distribution[:,2]
        OUP_potential = readin("../data/potential/OUP$n.csv")
        xs = OUP_potential[:,1] 
        Vs = OUP_potential[:,2]
        # Read the unnormalised gaussian samples
        histogram_distribution = readin("../data/distribution/histogram$n.csv")
        gaussian_bins = histogram_distribution[:,1]
        weights = histogram_distribution[:,2]
        histogram_potential = readin("../data/potential/histogram$n.csv")
        hist_xs = histogram_potential[:,1] 
        hist_Vs = histogram_potential[:,2]
        # Read the normalised gaussian samples
        gaussian_sample = readin("../data/distribution/gaussian$n.csv")
        gaussian_pdf = gaussian_sample[:,2]
        gaussian_potential = readin("../data/potential/gaussian$n.csv")
        gauss_xs = gaussian_potential[:,1] 
        gauss_Vs = gaussian_potential[:,2]

        # Get number of bins in the histograms 
        Nbins = length(bins)

        ##################################################
        #          Distribution of samples 
        ##################################################

        # Leftmost plot (unnormalised histogram)
        fig, ax = mkfig(size = [2000,1000],
                        bg_out = :white,
                        pad = (100,60,10,20),
                        box_position = [1,1],
                        lab = [L"\mathbf{x}", L"\mathbf{p(x)}"],
                        toggle_lab = [false, true],
                        lab_pad = [-60.0, -150.0],
                        toggle_ticks_lab = [false, true],
                        ticks_lab_trunc = [1,2] 
                       )
        # Customise the ticks of the subplot
        set_ticks(ax, gaussian_bins, weights)
        # Plot the points of the histogram
        scatter!(ax, gaussian_bins, weights, markersize = 20, color = :blue)
        # Plot vertical dashed lines from the x-axis to the points
        for m in 1:Nbins
                lines!(ax, [gaussian_bins[m], gaussian_bins[m]], [0.0, weights[m]], color = :blue, linestyle = :dash)
        end
        # Plot the edges of the bins
        scatter!(ax, LinRange(a - 0.05*I, b + 0.05*I, Nbins+1), zeros(Nbins+1), markersize = 5, color = :red)

        # Centre plot (normalised histogram from gaussian samples)
        nullfig, ax = mkfig(fig=fig,
                            box_position = [1,2],
                            toggle_lab = [false, false],
                            toggle_ticks_lab = [false, true],
                           )
        # Customise the ticks of the subplot
        set_ticks(ax, gaussian_bins, gaussian_pdf)
        # Plot the points of the histogram
        scatter!(ax, gaussian_bins, gaussian_pdf, markersize = 20, color = :blue)
        # Plot vertical dashed lines from the x-axis to the points
        for n in 1:Nbins
                lines!(ax, [gaussian_bins[n], gaussian_bins[n]], [0.0, gaussian_pdf[n]], color = :blue, linestyle = :dash)
        end
        # Plot the edges of the bins
        scatter!(ax, LinRange(a - 0.05*I, b + 0.05*I, Nbins+1), zeros(Nbins+1), markersize = 5, color = :red)
        # Actual Gaussian distribution
        X = LinRange(a - 0.05*I, b + 0.05*I, 1000)
        lines!(ax, X, [gaussian(x, μ[n], sqrt(D/ϴ)) for x in X], color = (:black, 0.5), linewidth = 3)

        # Rightmost plot (normalised histogram from OUP data)
        nullfig, ax = mkfig(fig=fig,
                            box_position = [1,3],
                            toggle_lab = [false, false],
                            toggle_ticks_lab = [false, true],
                           )
        # Customise the ticks of the subplot
        set_ticks(ax, bins, pdf)
        # Plot the points of the histogram
        scatter!(ax, bins, pdf, markersize = 20, color = :green)
        # Plot vertical dashed lines from the x-axis to the points
        for n in 1:Nbins
                lines!(ax, [bins[n], bins[n]], [0.0, pdf[n]], color = :green, linestyle = :dash)
        end
        # Plot the edges of the bins
        scatter!(ax, LinRange(a - 0.05*I, b + 0.05*I, Nbins+1), zeros(Nbins+1), markersize = 5, color = :red)
        # Actual Gaussian distribution
        lines!(ax, X, [gaussian(x, μ[n], sqrt(D/ϴ)) for x in X], color = (:black, 0.5), linewidth = 3)

        #################################################################
        #                       Potential function
        #################################################################

        # Leftmost plot (inverted unnormalised histogram)
        nullfig, ax = mkfig(fig = fig,
                            box_position = [2,1],
                            lab = [L"\mathbf{x}", L"\mathbf{V(x)}"],
                            toggle_lab = [true, true],
                            lab_pad = [-60.0, -60.0],
                            toggle_ticks_lab = [true, true],
                            ticks_lab_trunc = [1,2]
                           )
        # Customise the ticks of the subplot
        set_ticks(ax, hist_xs, hist_Vs)
        # Plot the points of the inverted histogram 
        scatter!(ax, hist_xs, hist_Vs, markersize = 20, color = :blue)
        # Plot vertical dashed lines from the x-axis to the points
        for m in 1:length(hist_xs)
                lines!(ax, [hist_xs[m], hist_xs[m]], [0.0, hist_Vs[m]], color = :blue, linestyle = :dash)
        end
        # Plot the edges of the bins
        scatter!(ax, LinRange(a - 0.05*I, b + 0.05*I, Nbins+1), zeros(Nbins+1), markersize = 5, color = :red)

        # Centre plot (inverted normalised histogram from gaussian samples)
        nullfig, ax = mkfig(fig = fig,
                            box_position = [2,2],
                            lab = [L"\mathbf{x}", L"\mathbf{V(x)}"],
                            toggle_lab = [true, false],
                            lab_pad = [-60.0, -60.0],
                            toggle_ticks_lab = [true, true],
                            ticks_lab_trunc = [1,2]
                           )
        # Customise the ticks of the subplot
        set_ticks(ax, gauss_xs, gauss_Vs)
        # Plot the points of the inverted histogram 
        scatter!(ax, gauss_xs, gauss_Vs, markersize = 20, color = :blue)
        # Plot vertical dashed lines from the x-axis to the points
        for m in 1:length(gauss_xs)
                lines!(ax, [gauss_xs[m], gauss_xs[m]], [0.0, gauss_Vs[m]], color = :blue, linestyle = :dash)
        end
        # Plot the edges of the bins
        scatter!(ax, LinRange(a - 0.05*I, b + 0.05*I, Nbins+1), zeros(Nbins+1), markersize = 5, color = :red)
        # Actual potential function
        U(x) = (x - μ[n])^2
        parabola = [U(x) for x in X]
        lines!(ax, X, parabola, color = (:black, 0.5), linewidth = 3)

        # Rightmost plot (inverted normalised histogram from the OUP data)
        nullfig, ax = mkfig(fig = fig,
                            box_position = [2,3],
                            lab = [L"\mathbf{x}", L"\mathbf{V(x)}"],
                            toggle_lab = [true, false],
                            lab_pad = [-60.0, -60.0],
                            toggle_ticks_lab = [true, true],
                            ticks_lab_trunc = [1,2]
                           )
        # Customise the ticks of the subplot
        set_ticks(ax, xs, Vs)
        # Plot the points of the inverted histogram 
        scatter!(ax, xs, Vs, markersize = 20, color = :green)
        # Plot vertical dashed lines from the x-axis to the points
        for m in 1:length(xs)
                lines!(ax, [xs[m], xs[m]], [0.0, Vs[m]], color = :green, linestyle = :dash)
        end
        # Plot the edges of the bins
        scatter!(ax, LinRange(a - 0.05*I, b + 0.05*I, Nbins+1), zeros(Nbins+1), markersize = 5, color = :red)
        # Actual potential function
        lines!(ax, X, parabola, color = (:black, 0.5), linewidth = 3)
 
        # Export figure
        save("../fig/distribution/$n.png", fig)
end
