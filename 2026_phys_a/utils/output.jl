"""
Export data int txt, csv and mat filetype.

Author: Davide Papapicco
Affil: U. of Auckland
Date: 02-12-2025
"""

"""
$(TYPEDSIGNATURES)

Export 'data' to local path and name it 'filename'.

## Keyword arguments
* `path=../../res/data`: path of the local tree to export the data

## Example
"""

function writeout(data, filename; path="../../res/data/")
        # Create the export directory if it does not exists
        fullpath = path * filename 
        mkpath(dirname(fullpath))

        # Write the input data to a csv file
        CSV.write(fullpath, Tables.table(data), delim=',', writeheader=false)
end
