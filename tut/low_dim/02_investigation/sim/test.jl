include("../../../../inc/PlottingTools.jl")
include("../../../../inc/IO.jl")

Nμ = convert(Int64,1e2)
μ = LinRange(0 ,1 ,Nμ)

x_inf = -1
x_sup = -x_inf 
x = LinRange(x_inf, x_sup, 1000)

printstyled("Running the loop using $(Threads.nthreads()) threads")
@showprogress Threads.@threads for n in 1:Nμ
        fig, ax = mkfig(size = [1200, 800],
                        limits = ((x_inf, x_sup), (-1, 1))
                       )
        for m in 1:n
                lines!(ax, x, μ[m].*x, linewidth = 8)
        end
        save("../fig/test/$n.png", fig)
end

printstyled("Running the loop using a single threads")
@showprogress for n in 1:Nμ
        fig, ax = mkfig(size = [1200, 800],
                        limits = ((x_inf, x_sup), (-1, 1))
                       )
        for m in 1:n
                lines!(ax, x, μ[m].*x, linewidth = 8)
        end
        save("../fig/test/$n.png", fig)
end

