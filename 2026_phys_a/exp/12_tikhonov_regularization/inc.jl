using BifurcationKit, LinearAlgebra, Polynomials, Integrals, Statistics
using CairoMakie, Makie.Colors, LaTeXStrings
using ProgressMeter

if !isdefined(Main, :DataInterface)
        include("../../src/DataInterface.jl")
        using .DataInterface
end

if !isdefined(Main, :SystemAnalysis)
        include("../../src/SystemAnalysis.jl")
        using .SystemAnalysis
end

if !isdefined(Main, :StatisticalMethods)
        include("../../src/StatisticalMethods.jl")
        using .StatisticalMethods
end
