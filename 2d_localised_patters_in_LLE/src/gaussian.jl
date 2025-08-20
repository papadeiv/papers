# Define the Gaussian distribution

function gaussian(x, μ, σ)
        D = σ^2/2
        C = sqrt(1/(4*pi*D))
        return C*exp(-(x-μ)^2/(4*D))
end
