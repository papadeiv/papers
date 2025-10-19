using Test

#-------------------------------#
#    1-dimensional evolution    #
#-------------------------------#

# Import the utilities for this test 
include("../utils/evolve_1d.jl")

# Run test set
@testset "evolve.jl" begin
        @test estimate_theta([1.0, 2.0, 3.0]) ≈ 2.0
        @test_throws ArgumentError estimate_theta([])
        #=
        # Initial condition 
        equilibria = get_equilibria(f, μ0, domain=[-10,10])
        x0 = [equilibria.stable[2], μ0]

        # Test 1: evolve 1-d in linear parameter range with give stepstize 
        sample_path = evolve(f, η, Λ, x0, endparameter=μf, stepsize=dt, particles=Ne)
        μ = sample_path.parameter
        display(μ[1])
        display(μ[end])

        # Test 2: evolve 1-d in linear parameter range within given number of steps  
        sample_path = evolve(f, η, Λ, x0, endparameter=μf, steps=Nt, particles=Ne)
        μ = sample_path.parameter
        display(μ[1])
        display(μ[end])

        # Test 3: evolve 1-d in time range from -T to T
        sample_path = evolve(f, η, Λ, x0, timerange=[-T,T], stepsize=dt, particles=Ne)
        t = sample_path.time
        display(t[1])
        display(t[end])

        # Test 4: evolve 1-d in time range from 0 to T
        sample_path = evolve(f, η, Λ, x0, timerange=T, stepsize=dt, particles=Ne)
        t = sample_path.time
        display(t[1])
        display(t[end])

        # Test 5: identify the tipping point
        sample_path = evolve(f, η, Λ, x0, endparameter=2*μf, stepsize=dt, particles=Ne)
        t = sample_path.time
        u = (sample_path.state)[2]
        tipping = find_tipping(u)
        display(t[tipping.index])
        =#
end

#-------------------------------#
#    2-dimensional evolution    #
#-------------------------------#

# Import the utilities for this test 
include("../utils/evolve_2d.jl")

# Run test set
@testset "evolve.jl" begin
        #=
        # Initial condition 
        equilibria = get_equilibria(f, g, μ0, guesses=[[sqrt(-μ0),0.0]])
        x0 = [equilibria.stable[1]; μ0]

        # Test 1: evolve 2-d in linear parameter range with given stepsize 
        sample_path = evolve([f, g], [η, η], Λ, x0, endparameter=μf, stepsize=dt, particles=Ne)
        μ = sample_path.parameter
        display(μ[1])
        display(μ[end])

        # Test 2: evolve 2-d in linear parameter range within given number of steps 
        sample_path = evolve([f, g], [η, η], Λ, x0, endparameter=μf, steps=Nt, particles=Ne)
        μ = sample_path.parameter
        display(μ[1])
        display(μ[end])

        # Test 3: evolve 2-d in time range from -T to T
        sample_path = evolve([f, g], [η, η], Λ, x0, timerange=[-T,T], stepsize=dt, particles=Ne)
        t = sample_path.time
        display(t[1])
        display(t[end])

        # Test 4: evolve 2-d in time range from 0 to T
        sample_path = evolve([f, g], [η, η], Λ, x0, timerange=T, stepsize=dt, particles=Ne)
        t = sample_path.time
        display(t[1])
        display(t[end])

        # Test 5: identify the tipping point
        sample_path = evolve([f, g], [η, η], Λ, x0, endparameter=μf+1.0, stepsize=dt, particles=Ne)
        t = sample_path.time
        u = (sample_path.state)[2]
        for component in eachcol(u)
                tipping = find_tipping(component)
                display(t[tipping.index])
        end
        =#
end
