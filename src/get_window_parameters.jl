# Assemble the sliding window through the user's input (width given as a relative length of the total timeseries)
function get_window_parameters(Nt::Int64, width::Float64)
        # Compute the length of the window
        Nw = convert(Int64, floor(width*Nt))
        # Compute the number of strides
        Ns = (Nt - Nw) + 1::Int64
        return [Nw, Ns]
end
