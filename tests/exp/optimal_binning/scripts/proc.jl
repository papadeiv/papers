"""
    Postprocessing script

In here we define the quantities related to the computation of EWSs from raw data.
"""

# Number of bins in the histogram
Nb = convert(Int64, floor(0.001*Nt))
