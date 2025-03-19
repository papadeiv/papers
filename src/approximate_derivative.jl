# Approximate the derivative using the forward finite difference
function approximate_derivative(data_x, data_f)
        df = Float64[]
        for n in 1:(length(data_f)-1)
                h = data_x[n+1]-data_x[n]
                push!(df, (data_f[n+1]-data_f[n])/h)
        end
        return df
end
