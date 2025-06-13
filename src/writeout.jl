function writeout(data, filename)
        CSV.write(filename, Tables.table(data), delim=',', writeheader=false)
        return 
end
