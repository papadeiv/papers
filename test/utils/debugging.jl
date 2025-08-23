"""
Short description. 

Author: Davide Papapicco
Affil: U. of Auckland
Date: 21-08-2025
"""

#------------------#
#                  # 
#   debugging.jl   #                     
#                  #
#------------------#

"""
    debug()

???
???
"""
function debug(raw_input)
        display(raw_input)
end

# Plot with formatting
function debug(label::String, variable)
        println(label, variable)
end

# Plot a sequence of elements with formatting
function debug(label::Vector{String}, variable)
        # Loop over the vector of entries until the penultimate element
        for n in 1:(length(label)-1)
                print(label[n], variable[n], ",  ")
        end

        # Print the last element
        println(label[end], variable[end])
end

# Plot a numbered sequence of elements with formatting
function debug(idx::Int, label::Vector{String}, variable)
        # Print the integer index
        print(idx, ")   ")

        # Loop over the vector of entries until the penultimate element
        for n in 1:(length(label)-1)
                print(label[n], variable[n], ",  ")
        end

        # Print the last element
        println(label[end], variable[end])
end
