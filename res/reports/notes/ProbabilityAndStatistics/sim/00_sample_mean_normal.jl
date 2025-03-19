include("../../../../inc/PlottingTools.jl")
include("../../../../inc/PotentialLearning.jl")

# Look for distributions to sample from
# https://juliastats.org/Distributions.jl/stable/univariate/#Continuous-Distributions 

# Define the number of samples of the random variable
N_samples = convert(Int64, 5e3)

# Define the number of bins for the sampling distribution of the sample mean statistic 
N_bins = convert(Int64, 5e1)

# Define the number of observations in each sample 
N_obs = convert(Int64,1e3)

# Define the number of bins for the distribution of observations 
N_bins_obs = convert(Int64, 2e1)

# Define the std of the sampling distribution
σ = 1.000

# Create arrays to store the observations and statistics of each sample 
z = Matrix{Float64}(undef, N_samples, N_obs)
T = Matrix{Float64}(undef, N_samples, 3)

# Create arrays to store the histograms sample 
bins = Matrix{Float64}(undef, N_samples, N_bins_obs)
pdf = Matrix{Float64}(undef, N_samples, N_bins_obs)

# Generate the samples 
using ProgressMeter
@showprogress for m in 1:N_samples
        # Generate the sample
        z[m,:] = σ*randn(N_obs)

        # Fit a distribution to the observations in sample
        bins[m,:], pdf[m,:] = fit_distribution(z[m,:], n_bins=N_bins_obs+1)
        
        # Compute the sample mean statistic 
        T[m,1] = mean(z[m,:])

        # Compute the std of the samples
        T[m,2] = std(z[m,:])

        # Compute the variance of the samples
        T[m,3] = var(z[m,:])
end

# Define different numbers of samples at which to plot the sampling distributions
samples = [100, 200, 500, N_samples]

# Create the figure and customise the axes for the observations distribution
fig, top_axs = rowfig(length(samples), 1,
                      size = [4000,1500],
                      bg_out = :transparent)
# Create the figure and customise the axes for the sampling distribution of the sample mean
nullfig, bottom_axs = rowfig(length(samples), 2,
                  fig = fig)

# Loop over the numbers of samples
@showprogress for s in 1:length(samples)
        # Remove the y-label from the top and bottom axes
        top_axs[s].ylabelvisible = false
        bottom_axs[s].ylabelvisible = false
       # Plot the distribution of the observations at different numbers of samples
        for n in 1:samples[s]
                lines!(top_axs[s], bins[n,:], pdf[n,:], color = (:black, 0.05))
        end
        # Fit a distribution to the sample mean statistic at different numbers of samples
        bins_μ, pdf_μ = fit_distribution(T[1:samples[s],1], n_bins=N_bins+1)
        # Set the ticks for the bottom axes
        set_ticks(bottom_axs[s], bins_μ, pdf_μ, n_ticks=1)
        # Plot the sampling distribution of the sample mean statistic
        lines!(bottom_axs[s], bins_μ, pdf_μ, color = (:red, 1.0), linewidth = 5)
end

# Export the figure
save("../fig/sampling_distribution.png", fig)

# Compute the mean of the sampling distribution of the sample mean statistic at increasing numbers of samples
μ = [mean(T[1:m,1]) for m in 1:N_samples]

# Create the figure and customise axes for the statistic mean series
fig, ax = mkfig(size = [2000,1000], bg_out = :transparent,
                lab = [L"\textbf{number of samples (\mathbf{m})}",L"\mathbf{\mu_{\overline{X}}}"],
                ticks_lab_trunc = [0,1]
               )
set_ticks(ax, LinRange(1, N_samples, N_samples), μ./1e-3, n_ticks=5)
# Plot the mean statistic series
lines!(ax, LinRange(1, N_samples, N_samples), μ./1e-3, color = :red, linewidth = 5)
# Plot the population mean
lines!(ax, [1, N_samples], [0,0], color = :black, linewidth = 3, linestyle = :dash)
text!(800, 1, text = L"\mathbf{\mu_{X}}", fontsize = 50)
# Add the scale label
Label(fig[1, 1, Top()], halign = :left, L"\mathbf{\times 10^{-3}}", fontsize = 50)

# Export the figure
save("../fig/fig2.3.2.2.png", fig)

# Compute the std of the sampling distribution of the sample mean statistic at increasing numbers of samples
s = [std(T[1:m,1]) for m in 1:N_samples]

# Create the figure and customise axes for the statistic std series
fig, ax = mkfig(size = [2000,1000], bg_out = :transparent,
                lab = [L"\textbf{number of samples (\mathbf{m})}",L"\mathbf{\sigma_{\overline{X}}}"],
                x_ticks = [1,834,1667,2500,3334,4167,5000],
                ticks_lab_trunc = [0,1]
               )
# Plot the mean statistic series
lines!(ax, LinRange(1, N_samples, N_samples), s./1e-2, color = :red, linewidth = 5)
# Plot the population mean
lines!(ax, [1, N_samples], [σ/sqrt(N_obs),σ/sqrt(N_obs)]./1e-2, color = :black, linewidth = 3, linestyle = :dash)
text!(800, 2.5, text = L"\mathbf{\frac{\sigma_{X}}{\sqrt{n}}}", fontsize = 50)
# Add the scale label
Label(fig[1, 1, Top()], halign = :left, L"\mathbf{\times 10^{-2}}", fontsize = 50)

# Export the figure
save("../fig/fig2.3.2.3.png", fig)


# Compute the standard error of the sample mean (https://en.wikipedia.org/wiki/Standard_error#Standard_error_of_the_mean)
#=
display(std(μ))
display(σ/sqrt(N_obs))
=#
