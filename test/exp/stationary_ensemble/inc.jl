using ProgressMeter 
using Revise

# Avoid re-loading SimpleIO
if !isdefined(Main, :SimpleIO)
        include("../../src/SimpleIO.jl")
        using .SimpleIO
end

# Avoid re-loading EscapeProblem 
if !isdefined(Main, :EscapeProblem)
        include("../../src/EscapeProblem.jl")
        using .EscapeProblem
end

# Avoid re-loading SystemAnalysis 
if !isdefined(Main, :SystemAnalysis)
        include("../../src/SystemAnalysis.jl")
        using .SystemAnalysis
end
