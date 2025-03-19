using Statistics

# Computes the temporal variance of a timeseries over a sliding window of specified width
function variance(time, data, width::Float64)
        sliding_window = trunc(Int64, width*(length(data)))
        time_ews = time[(sliding_window+1):end] 
        ews = [var(data[n:n+sliding_window]) for n in 1:(length(data) - sliding_window)]
        return time_ews, ews
end
