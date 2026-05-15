"""
    Plotting script

Functions used to create the plots in each figure.
"""

# Plot the timeseries decomposition and their distributions
function plot_imfs(timestep, timeseries, imfs)
        #---------------------------------------------------#
        # Residuals by removing IMFS[11]                    #
        #---------------------------------------------------#
        
        trend = imfs[end,:]
        residual = timeseries .- trend

        fig = Figure()

        ax1 = Axis(fig[1,1:2], title="Original signal w/ trend (IMFS[11])")
        lines!(ax1, timestep, timeseries, color = :black, linewidth = 1.0)
        lines!(ax1, timestep, trend, color = :red, linewidth = 3.0)

        ax2 = Axis(fig[1,3], title="Empirical distribution")
        Nb = convert(Int64, ceil(abs(maximum(timeseries)-minimum(timeseries))/(3.49*std(timeseries)*(length(timeseries))^(-1.0/3.0))))
        bins, pdf = fit_distribution(timeseries, n_bins=Nb+1)
        barplot!(ax2, bins, pdf, color = pdf, colormap = [:white,CtpGray], strokecolor = :black, strokewidth = 1)

        ax3 = Axis(fig[2,1:2], title="Detrended signal")
        lines!(ax3, timestep, residual, color = :black, linewidth = 1.0)

        ax4 = Axis(fig[2,3], title="Empirical distribution")
        Nb = convert(Int64, ceil(abs(maximum(residual)-minimum(residual))/(3.49*std(residual)*(length(residual))^(-1.0/3.0))))
        bins, pdf = fit_distribution(residual, n_bins=Nb+1)
        barplot!(ax4, bins, pdf, color = pdf, colormap = [:white,CtpGray], strokecolor = :black, strokewidth = 1)

        savefig("cubic/residuals/residuals_1.png", fig)

        #---------------------------------------------------#
        # Residuals by removing IMFS[11+n]                  #
        #---------------------------------------------------#
        
        # List of titles
        titles = ["Original signal w/ trend (IMFS[11+10])", 
                  "Original signal w/ trend (IMFS[11+10+9])",
                  "Original signal w/ trend (IMFS[11+10+9+8])",
                  "Original signal w/ trend (IMFS[11+...+7])",
                  "Original signal w/ trend (IMFS[11+...+6])",
                  "Original signal w/ trend (IMFS[11+...+5])",
                  "Original signal w/ trend (IMFS[11+...+4])",
                  "Original signal w/ trend (IMFS[11+...+3])",
                  "Original signal w/ trend (IMFS[11+...+2])",
                  "Original signal w/ trend (IMFS[11+...+1])",
                 ]

        # Loop over the IMFs
        for n in 1:(size(imfs,1)-1)
                # Compute the trend and the residuals
                trend = sum(imfs[(end-n):end,:], dims=1)[:]
                residual = timeseries .- trend

                # Plot the figure
                fig = Figure()

                ax1 = Axis(fig[1,1:2], title=titles[n])
                lines!(ax1, timestep, timeseries, color = :black, linewidth = 1.0)
                lines!(ax1, timestep, trend, color = :red, linewidth = 3.0)

                ax2 = Axis(fig[1,3], title="Empirical distribution")
                Nb = convert(Int64, ceil(abs(maximum(timeseries)-minimum(timeseries))/(3.49*std(timeseries)*(length(timeseries))^(-1.0/3.0))))
                bins, pdf = fit_distribution(timeseries, n_bins=Nb+1)
                barplot!(ax2, bins, pdf, color = pdf, colormap = [:white,CtpGray], strokecolor = :black, strokewidth = 1)

                ax3 = Axis(fig[2,1:2], title="Detrended signal")
                lines!(ax3, timestep, residual, color = :black, linewidth = 1.0)

                ax4 = Axis(fig[2,3], title="Empirical distribution")
                Nb = convert(Int64, ceil(abs(maximum(residual)-minimum(residual))/(3.49*std(residual)*(length(residual))^(-1.0/3.0))))
                bins, pdf = fit_distribution(residual, n_bins=Nb+1)
                barplot!(ax4, bins, pdf, color = pdf, colormap = [:white,CtpGray], strokecolor = :black, strokewidth = 1)

                # Export the figure
                savefig("cubic/residuals/residuals_$(n+1).png", fig)
        end

        #---------------------------------------------------#
        # Small scale distributions                         #
        #---------------------------------------------------#
        
        fig = Figure()

        ax1 = Axis(fig[1,1:2], title="IMFS[1] (smallest scale)")
        lines!(ax1, timestep, imfs[1,:], color = :black, linewidth = 1.0)

        ax2 = Axis(fig[1,3], title="Empirical distribution")
        Nb = convert(Int64, ceil(abs(maximum(imfs[1,:])-minimum(imfs[1,:]))/(3.49*std(imfs[1,:])*(length(imfs[1,:]))^(-1.0/3.0))))
        bins, pdf = fit_distribution(imfs[1,:], n_bins=Nb+1)
        barplot!(ax2, bins, pdf, color = pdf, colormap = [:white,CtpGray], strokecolor = :black, strokewidth = 1)

        ax3 = Axis(fig[2,1:2], title="IMFS[2] (second smallest scale)")
        lines!(ax3, timestep, imfs[2,:], color = :black, linewidth = 1.0)

        ax4 = Axis(fig[2,3], title="Empirical distribution")
        Nb = convert(Int64, ceil(abs(maximum(imfs[2,:])-minimum(imfs[2,:]))/(3.49*std(imfs[2,:])*(length(imfs[2,:]))^(-1.0/3.0))))
        bins, pdf = fit_distribution(imfs[2,:], n_bins=Nb+1)
        barplot!(ax4, bins, pdf, color = pdf, colormap = [:white,CtpGray], strokecolor = :black, strokewidth = 1)

        ax5 = Axis(fig[3,1:2], title="IMFS[3] (third smallest scale)")
        lines!(ax5, timestep, imfs[3,:], color = :black, linewidth = 1.0)

        ax6 = Axis(fig[3,3], title="Empirical distribution")
        Nb = convert(Int64, ceil(abs(maximum(imfs[3,:])-minimum(imfs[3,:]))/(3.49*std(imfs[3,:])*(length(imfs[3,:]))^(-1.0/3.0))))
        bins, pdf = fit_distribution(imfs[3,:], n_bins=Nb+1)
        barplot!(ax6, bins, pdf, color = pdf, colormap = [:white,CtpGray], strokecolor = :black, strokewidth = 1)

        savefig("cubic/residuals/modes.png", fig)

        #---------------------------------------------------#
        # Detailed small scale oscillations                 #
        #---------------------------------------------------#
        
        # Final timestep
        Nt = 200
                
        fig = Figure()

        ax1 = Axis(fig[1,1], title="Original signal")
        lines!(ax1, timestep[1:Nt], timeseries[1:Nt], color = :black, linewidth = 1.0)

        ax2 = Axis(fig[2,1], title="IMFS[1] (smallest scale)")
        lines!(ax2, timestep[1:Nt], imfs[1,1:Nt], color = :black, linewidth = 1.0)

        ax3 = Axis(fig[3,1], title="IMFS[2] (second smallest scale)")
        lines!(ax3, timestep[1:Nt], imfs[2,1:Nt], color = :black, linewidth = 1.0)

        savefig("cubic/residuals/detailed.png", fig)
end
