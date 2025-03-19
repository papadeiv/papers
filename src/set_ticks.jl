using LaTeXStrings, CairoMakie, Makie.Colors
CairoMakie.activate!(; px_per_unit = 2)

# Customise the ticks for each axes in the plot
function set_ticks(axis, x_data, y_data; n_ticks = 0)
        # Extract the range in the data
        x_min = x_data[argmin(x_data)] 
        x_max = x_data[argmax(x_data)]
        y_min = y_data[argmin(y_data)]
        y_max = y_data[argmax(y_data)]

        # Check if there are additional ticks specified by the user
        if n_ticks == 0
                axis.xticks = [x_min, x_max]
                axis.yticks = [y_min, y_max] 
        else
                # Find the values of the intermediate ticks
                axis.xticks = LinRange(x_min, x_max, n_ticks+2)
                axis.yticks = LinRange(y_min, y_max, n_ticks+2)
        end

        # Return the axis structure
        return axis
end
