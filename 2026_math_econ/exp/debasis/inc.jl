using CairoMakie, Makie.Colors, LaTeXStrings
using ProgressMeter, Revise
using Printf

using NonlinearSolve 
using DynamicalSystems

# Avoid re-loading PlottingTools
if !isdefined(Main, :PlottingTools)
        include("../../src/PlottingTools.jl")
        using .PlottingTools
end
