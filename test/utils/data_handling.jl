"""
I/O utility functions to simplify the process of reading and writing data.

Author: Davide Papapicco
Affil: U. of Auckland
Date: 21-08-2025
"""

#----------------------#
#                      # 
#   data_handling.jl   #                     
#                      #
#----------------------#

"""
    readCSV(filename)

Import data from a .csv `filename` without header and as a Matrix.
"""
function readCSV(filename)
        # Read the data in dataframe form
        df = DataFrame(CSV.File(filename; delim=',', header=false))

        # Define the time of datatype to be returned based on the number of columns in the dataframe
        if ncol(df) == 1 
                # Convert the dataframe in vector form
                global data = Vector{Float64}(undef, nrow(df))
                data = df[:,1]
        else
                # Convert the dataframe in matrix form
                data = Matrix{Float64}(undef, nrow(df), ncol(df))
                for n in 1:ncol(df)
                        data[:,n] = df[:,n]
                end
        end

        # Export the array
        return data 
end

"""
    writeCSV(data, filename)

Export a Matrix `data` into a .csv `filename` without header.
"""
function writeCSV(data, filename)
        # Create the directory if it does not exists
        mkpath(dirname(filename))
        # Write the input data to a csv file
        CSV.write(filename, Tables.table(data), delim=',', writeheader=false)
        return 
end
