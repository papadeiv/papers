# Import the necessary packages and local modules
using LinearAlgebra, DifferentialEquations
using Roots, ForwardDiff
using DocStringExtensions

include("../../utils/equilibria.jl")
include("../../utils/preprocess.jl")
include("../../utils/evolution.jl")

# System parameters
μ0 = -(1.00::Float64)
μf = 0.00::Float64
ε = 1e-2
σ = 0.010::Float64
D = (σ^2)/2.0
dt = 1e-2
Ne = 10::Int
T = 10.0::Float64
Nt = 1e3

# Dynamical system  
f(x, y, μ) = - μ - x^2
g(x, y, μ) = - y 
Λ(t) = ε
η(x, y) = σ
