# Compute and plot the ensemble's average escape rate
function ensemble_escape_rate(time::Float64, distribution, Ne::Int64, treshold::Int64)
        # Extract the number of N-tippings from the ensemble
        N_exit = length(distribution) + 1

        # Initialise the ensemble escape rate
        rate = 0.0::Float64

        # Define a boolean variable to check the condition of satisftying the treshold 
        criterion = 0::Int64

        # Compute the escape rate based on the specified treshold (minimum number of N-tipping events) 
        if N_exit < treshold
                criterion = 1::Int64
                rate = (N_exit/Ne)/time
        else
                rate = (treshold/Ne)/distribution[treshold-1]
        end

        # Return both the rate and the boolean criterion
        return [rate, criterion]
end
