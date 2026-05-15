"""
    Main script

Run this file to execute the simulation, analyse and plot the results.
"""

# Import the necessary packages and local modules
include("inc.jl")

# Import the simulation's scripts
include("./scripts/proc.jl")
include("./scripts/plot.jl")

# Define the main algorithm
function main()
        # Import the AMOC timeseries
        data = NCDataset("../../res/data/amoc/AMOC_transport_depth_0-1000m_monthly.nc")
        time = data["time"][:]
        observable = data["Transport"][:]

        # Convert the analysis of the non-autonomous drift into an ensemble problem
        ensemble = preprocess_solution(time, observable, window_size)
        tipping = ensemble.tipping_point
        Ne = length(ensemble.trajectories)
        Nt = length(ensemble.trajectories[1])

        # Extract windowed subseries
        t = copy(ensemble.timesteps[1])
        u = copy(ensemble.trajectories[1]) 

        # Detrend it 
        PyEMD = pyimport("PyEMD")
        emd = PyEMD.EMD()
        imfs = Array(emd(u))
        display(size(imfs))

        # Plot the timeseries decomposition and their distribution
        plot_imfs(t, u, imfs)
end

# Execute the main
main()
