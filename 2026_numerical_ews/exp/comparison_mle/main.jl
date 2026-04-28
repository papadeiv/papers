"""
    Main script

Run this file to execute the simulation, analyse and plot the results.
"""

# Import the necessary packages and local modules
include("inc.jl")

# Import the simulation's scripts
include("./scripts/sim.jl")
include("./scripts/proc.jl")

# Main algorithm 
function main()
        # Loop over the fixed parameters
        for μ in μ_set[1]
                # Solve the ensemble problem 
                equilibria = get_equilibria(f, μ, domain=[-10,10])
                x0 = [equilibria.stable[2], μ]
                ensemble = evolve(f, η, Λ, x0, stepsize=dt, steps=Nt, particles=Ne)

                # Loop over the ensemble's sample paths
                printstyled("Solving the inverse problems over the ensemble for μ = $μ\n"; bold=true, underline=true, color=:light_blue)
                @showprogress for n in 1:length(ensemble.state)
                        # Extract the trajectory and center it
                        u = ensemble.state[n] .- x0[1] 

                        #---------------------------#
                        #  F-P nonlinear inversion  #          
                        #---------------------------#
         
                        solution = fit_potential(u, noise=σ, transformation=[0.0,1.0,8.0])
                        V_NLLS = shift_potential(x0[1], solution.fit, μ)
                        display(solution.fit) error[n,1] = get_error(V_NLLS, equilibria, μ)

                        #----------------#
                        #  EM Quasi MLE  #          
                        #----------------#
                        
                        Xn = u[1:end-1]
                        Y  = (u[2:end] .- u[1:end-1])./dt
                        Φ = hcat(ones(length(Xn)), Xn, Xn.^2)
                        β = Φ\Y
                        c = [-β[1], -β[2]/2, -β[3]/3]
                        V_EM = shift_potential(x0[1], c, μ)
                        error[n,2] = get_error(V_EM, equilibria, μ)

                        #------------------------#
                        #  F-P linear inversion  #          
                        #------------------------#
         
                        n_bins = convert(Int64, ceil(abs(maximum(u)-minimum(u))/(3.49*std(u)*(length(u))^(-1.0/3.0))))
                        bins, pdf = fit_distribution(u, n_bins=n_bins)
                        idx = findall(x -> x > 0.0, pdf)
                        y = [pdf[n] for n in idx]
                        x = [bins[n] for n in idx]
                        N=1/sqrt(2*pi*D)
                        U_discrete = -D .* log.(y ./ N)
                        Φ = hcat(ones(length(x)), x, x.^2, x.^3)
                        W = Diagonal(y)
                        c = ((Φ'*W*Φ)\(Φ'*W*U_discrete))[2:4]
                        V_LLS = shift_potential(x0[1], c, μ)
                        error[n,3] = get_error(V_LLS, equilibria, μ)

                        #-------------------#
                        #  Generalised MoM  #          
                        #-------------------#
                        
                        m = [mean(u.^k) for k in 0:5]
                        A = [
                             m[2]  2*m[3]  3*m[4];
                             m[3]  2*m[4]  3*m[5];
                             m[4]  2*m[5]  3*m[6]
                            ]
                        b = D.*[m[1]; 2*m[2]; 3*m[3]]
                        c = A\b
                        V_MOM = shift_potential(x0[1], c, μ)
                        error[n,4] = get_error(V_MOM, equilibria, μ)
                end

                # Plot distribution of errors
                fig = Figure()
                methods = ["F-P NLLS", "E-M QMLE", "F-P LLS", "G MoM"]
                ax = Axis(fig[1,1], limits = (-0.05, 0.5, nothing, nothing), xlabel = L"\mathbf{||V - U||_2/||V||_2}", yticks = ((1:4).*(4*μ_idx), methods), title = "μ=$μ")
                maximum_element = maximum(error[:,1])
                for n in size(error, 2):-1:1
                        # Filter out the outliers
                        error_method = error[:,n]
                        p90 = quantile(error_method, 0.9)
                        filtered_error = error_method[error_method .<= p90]
                        density!(ax, filtered_error, offset = (4*μ_idx)*n, color = :x, colormap = :hawaii, colorrange = (0,maximum_element), strokewidth = 1, strokecolor = :black)
                end
                savefig("../../res/fig/μ=$μ.png", fig)
                global μ_idx +=1
        end
end

# Execute the main 
main()
