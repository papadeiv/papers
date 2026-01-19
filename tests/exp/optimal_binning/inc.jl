using CairoMakie, Makie.Colors, LaTeXStrings
using Distributions, Statistics, StatsBase
using Polynomials, LinearAlgebra, LsqFit 
using ProgressMeter, Revise
using DataFrames, CSV
using Printf

# Avoid re-loading SystemAnalysis 
if !isdefined(Main, :SystemAnalysis)
        include("../../src/SystemAnalysis.jl")
        using .SystemAnalysis
end

# Avoid re-loading StatisticalMethods
if !isdefined(Main, :StatisticalMethods)
        include("../../src/StatisticalMethods.jl")
        using .StatisticalMethods
end

# Avoid re-loading PlottingTools
if !isdefined(Main, :PlottingTools)
        include("../../src/PlottingTools.jl")
        using .PlottingTools
end
