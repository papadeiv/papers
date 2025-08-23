# Standard output method

function writeout(data, filename)
        # Create the directory if it does not exists
        mkpath(dirname(filename))
        # Write the input data to a csv file
        CSV.write(filename, Tables.table(data), delim=',', writeheader=false)
        return 
end
