"""
    Main script

Run this file to execute the simulation, analyse and plot the results.
"""

# Import the necessary packages and local modules
include("inc.jl")

# Import the simulation's scripts
include("./scripts/sim.jl")
include("./scripts/proc.jl")
include("./scripts/plot.jl")

# Samples generation for good and bad solutions
function generate_samples()
        # Define the initial condition 
        equilibria = get_equilibria(f, μ, domain=[-10,10])
        x0 = [equilibria.stable[2], μ]

        # Solve the ensemble problem 
        ensemble = evolve(f, η, Λ, x0, stepsize=dt, steps=Nt, particles=Ne)

        # Loop over the ensemble's sample paths
        printstyled("Looping over the ensemble\n"; bold=true, underline=true, color=:light_blue)
        for n in 1:length(ensemble.state)
                # Extract the trajectory until the tipping and center it 
                u = ensemble.state[n] .- x0[1] 
                t = ensemble.time

                # Solve the nonlinear least-squares problem to fit a cubic potential
                solution = fit_potential(u, noise=σ, transformation=[0.0,1.0,8.0])

                # Compute the error in the reconstruction
                Vs = shift_potential(x0[1], solution.fit)
                error = get_error(Vs, equilibria)

                # Export the "good" and "bad" solutions 
                if error < 0.05 && glb_idx_good <= 15
                        println("good n.$glb_idx_good) error = $error")
                        writeout(u, "good/$glb_idx_good.csv")
                        global glb_idx_good += 1
                elseif error > 0.2 && glb_idx_bad <= 15
                        println("bad n.$glb_idx_bad) error = $error")
                        writeout(u, "bad/$glb_idx_bad.csv")
                        global glb_idx_bad += 1
                end
        end

        # Return the stable (observed) equilibrium and number of samples
        return (
                equilibrium = x0[1],
                n_samples = min(glb_idx_good-1, glb_idx_bad-1)
               )
end

# Main algorithm
function main(equilibrium, N_exp)
         # Loop over the exported samples
        for n in 1:N_exp
                println("Sample n. $n")
                # Figures for potential and distribution 
                fig1 = Figure(size = (1000, 500))
                ax11 = Axis(fig1[1,1])
                ax12 = Axis(fig1[1,2], limits = (-2,2,-2,2))
                fig2 = Figure(size = (1000, 500))
                ax21 = Axis(fig2[1,1])
                ax22 = Axis(fig2[1,2], limits = (-2,2,-2,2))
                domain = LinRange(-2,2,1000)

                # Import the good and bad samples
                good = readin("/good/$n.csv")
                bad = readin("/bad/$n.csv")

                # Plot the histogram of the samples
                n_bins = convert(Int64, ceil(abs(maximum(good)-minimum(good))/(3.49*std(good)*(length(good))^(-1.0/3.0))))
                bins, pdf = fit_distribution(good, n_bins=n_bins)
                barplot!(ax11, bins, pdf, color = pdf, colormap = [(CtpYellow,0.5),(CtpMauve,0.5)], strokecolor = :black, strokewidth = 1.0)

                n_bins = convert(Int64, ceil(abs(maximum(bad)-minimum(bad))/(3.49*std(bad)*(length(bad))^(-1.0/3.0))))
                bins, pdf = fit_distribution(bad, n_bins=n_bins)
                barplot!(ax21, bins, pdf, color = pdf, colormap = [(CtpYellow,0.5),(CtpMauve,0.5)], strokecolor = :black, strokewidth = 1.0)

                # Plot the ground truth potential
                lines!(ax12, domain, [U(x, μ) for x in domain], color = :black, linewidth = 3.0)
                lines!(ax22, domain, [U(x, μ) for x in domain], color = :black, linewidth = 3.0)

                #----------------#
                #      NLLS      #          
                #----------------#
                
                # Good solution
                solution = fit_potential(good, noise=σ, transformation=[0.0,1.0,8.0])
                V_NLLS = shift_potential(equilibrium, solution.fit)
                lines!(ax12, domain, [V_NLLS(x) for x in domain], color = :red, linewidth = 2.0)

                I = (-Inf,Inf)
                if xs(solution.fit) > xu(solution.fit)
                        I = (xu(solution.fit), +Inf) 
                else
                        I = (-Inf, xu(solution.fit)) 
                end
                n_bins = convert(Int64, ceil(abs(maximum(good)-minimum(good))/(3.49*std(good)*(length(good))^(-1.0/3.0))))
                bins, pdf = fit_distribution(good, n_bins=n_bins)
                binning = LinRange(minimum(bins), maximum(bins), 1000)
                integral = IntegralProblem(ρ, I, solution.fit)
                quadrature = solve(integral, QuadGKJL(; order=20000); maxiters=10000)
                p(x, c) = ρ(x, c)/(quadrature.u)
                lines!(ax11, binning, [p(x, solution.fit) for x in binning], color = :red, linewidth = 4.0)

                # Bad solution
                solution = fit_potential(bad, noise=σ, transformation=[0.0,1.0,8.0])
                V_NLLS = shift_potential(equilibrium, solution.fit)
                lines!(ax22, domain, [V_NLLS(x) for x in domain], color = :red, linewidth = 2.0)

                I = (-Inf,Inf)
                if xs(solution.fit) > xu(solution.fit)
                        I = (xu(solution.fit), +Inf) 
                else
                        I = (-Inf, xu(solution.fit)) 
                end
                n_bins = convert(Int64, ceil(abs(maximum(bad)-minimum(bad))/(3.49*std(bad)*(length(bad))^(-1.0/3.0))))
                bins, pdf = fit_distribution(bad, n_bins=n_bins)
                binning = LinRange(minimum(bins), maximum(bins), 1000)
                integral = IntegralProblem(ρ, I, solution.fit)
                quadrature = solve(integral, QuadGKJL(; order=20000); maxiters=10000)
                p(x, c) = ρ(x, c)/(quadrature.u)
                lines!(ax21, binning, [p(x, solution.fit) for x in binning], color = :red, linewidth = 4.0)

                #----------------#
                #       LLS      #
                #----------------#
                
                # Good
                Xn = good[1:end-1]
                Y  = (good[2:end] .- good[1:end-1])./dt
                Φ = hcat(ones(length(Xn)), Xn, Xn.^2)
                β = Φ\Y
                c = [-β[1], -β[2]/2, -β[3]/3]
                V_LLS = shift_potential(equilibrium, c)
                lines!(ax12, domain, [V_LLS(x) for x in domain], color = :blue, linewidth = 2.0)

                I = (-Inf,Inf)
                if xs(c) > xu(c)
                        I = (xu(c), +Inf) 
                else
                        I = (-Inf, xu(c)) 
                end
                n_bins = convert(Int64, ceil(abs(maximum(good)-minimum(good))/(3.49*std(good)*(length(good))^(-1.0/3.0))))
                bins, pdf = fit_distribution(good, n_bins=n_bins)
                binning = LinRange(minimum(bins), maximum(bins), 1000)
                integral = IntegralProblem(ρ, I, solution.fit)
                quadrature = solve(integral, QuadGKJL(; order=20000); maxiters=10000)
                p(x, c) = ρ(x, c)/(quadrature.u)
                lines!(ax11, binning, [p(x, c) for x in binning], color = :blue, linewidth = 2.0)

                # Bad 
                Xn = bad[1:end-1]
                Y  = (bad[2:end] .- bad[1:end-1])./dt
                Φ = hcat(ones(length(Xn)), Xn, Xn.^2)
                β = Φ\Y
                c = [-β[1], -β[2]/2, -β[3]/3]
                V_LLS = shift_potential(equilibrium, c)
                lines!(ax22, domain, [V_LLS(x) for x in domain], color = :blue, linewidth = 2.0)

                I = (-Inf,Inf)
                if xs(c) > xu(c)
                        I = (xu(c), +Inf) 
                else
                        I = (-Inf, xu(c)) 
                end
                n_bins = convert(Int64, ceil(abs(maximum(bad)-minimum(bad))/(3.49*std(bad)*(length(bad))^(-1.0/3.0))))
                bins, pdf = fit_distribution(bad, n_bins=n_bins)
                binning = LinRange(minimum(bins), maximum(bins), 1000)
                integral = IntegralProblem(ρ, I, solution.fit)
                quadrature = solve(integral, QuadGKJL(; order=20000); maxiters=10000)
                p(x, c) = ρ(x, c)/(quadrature.u)
                lines!(ax21, binning, [p(x, c) for x in binning], color = :blue, linewidth = 2.0)

                #--------------------------------#
                # Invariant measure (histograms) #
                #--------------------------------#
                
                # Good
                n_bins = convert(Int64, ceil(abs(maximum(good)-minimum(good))/(3.49*std(good)*(length(good))^(-1.0/3.0))))
                bins, pdf = fit_distribution(good, n_bins=n_bins)
                idx = findall(x -> x > 0.0, pdf)
                y = [pdf[n] for n in idx]
                x = [bins[n] for n in idx]
                N=1/sqrt(2*pi*D)

                #=
                U_discrete = -D.*log.(y./N)
                c = Polynomials.fit(x, U_discrete, 3).coeffs[2:4]
                =#

                U_discrete = -D .* log.(y./N)
                Φ = hcat(ones(length(x)), x, x.^2, x.^3)
                W = Diagonal(y)
                c_full = (Φ'*W*Φ)\(Φ'*W*U_discrete)
                c = c_full[2:4]

                V_hist = shift_potential(equilibrium, c)
                lines!(ax12, domain, [V_hist(x) for x in domain], color = :darkviolet, linewidth = 2.0)

                I = (-Inf,Inf)
                if xs(c) > xu(c)
                        I = (xu(c), +Inf) 
                else
                        I = (-Inf, xu(c)) 
                end
                binning = LinRange(minimum(bins), maximum(bins), 1000)
                integral = IntegralProblem(ρ, I, c)
                quadrature = solve(integral, QuadGKJL(; order=20000); maxiters=10000)
                p(x, c) = ρ(x, c)/(quadrature.u)
                lines!(ax11, binning, [p(x, c) for x in binning], color = :darkviolet, linewidth = 2.0)

                # Bad 
                n_bins = convert(Int64, ceil(abs(maximum(bad)-minimum(bad))/(3.49*std(bad)*(length(bad))^(-1.0/3.0))))
                bins, pdf = fit_distribution(bad, n_bins=n_bins)
                idx = findall(x -> x > 0.0, pdf)
                y = [pdf[n] for n in idx]
                x = [bins[n] for n in idx]
                N=1/sqrt(2*pi*D)
                U_discrete = -D .* log.(y ./ N)
                Φ = hcat(ones(length(x)), x, x.^2, x.^3)
                W = Diagonal(y)
                c = ((Φ'*W*Φ)\(Φ'*W*U_discrete))[2:4]
                V_hist = shift_potential(equilibrium, c)
                lines!(ax22, domain, [V_hist(x) for x in domain], color = :darkviolet, linewidth = 2.0)

                I = (-Inf,Inf)
                if xs(c) > xu(c)
                        I = (xu(c), +Inf) 
                else
                        I = (-Inf, xu(c)) 
                end
                binning = LinRange(minimum(bins), maximum(bins), 1000)
                integral = IntegralProblem(ρ, I, c)
                quadrature = solve(integral, QuadGKJL(; order=20000); maxiters=10000)
                p(x, c) = ρ(x, c)/(quadrature.u)
                lines!(ax21, binning, [p(x, c) for x in binning], color = :darkviolet, linewidth = 2.0)

                #-------------------------#
                # Invariant measure (KDE) #
                #-------------------------#
                
                # Good
                kde_est = kde(good)
                bins = kde_est.x
                pdf = kde_est.density
                idx = findall(x -> x > 0.0, pdf)
                y = [pdf[n] for n in idx]
                x = [bins[n] for n in idx]
                N=1/sqrt(2*pi*D)
                U_discrete = -D .* log.(y ./ N)
                Φ = hcat(ones(length(x)), x, x.^2, x.^3)
                W = Diagonal(y)
                c = ((Φ'*W*Φ)\(Φ'*W*U_discrete))[2:4]
                V_kde = shift_potential(equilibrium, c)
                lines!(ax12, domain, [V_kde(x) for x in domain], color = :darkgreen, linewidth = 2.0)

                #------------------------#
                #         Moments        #
                #------------------------#
                
                # Good
                m = [mean(good.^k) for k in 0:5]
                A = [
                     m[2]  2*m[3]  3*m[4];
                     m[3]  2*m[4]  3*m[5];
                     m[4]  2*m[5]  3*m[6]
                    ]
                b = D.*[m[1]; 2*m[2]; 3*m[3]]
                c = A\b
                V_mom = shift_potential(equilibrium, c)
                lines!(ax12, domain, [V_mom(x) for x in domain], color = :gold, linewidth = 2.0)

                # Bad 
                m = [mean(bad.^k) for k in 0:5]
                A = [
                     m[2]  2*m[3]  3*m[4];
                     m[3]  2*m[4]  3*m[5];
                     m[4]  2*m[5]  3*m[6]
                    ]
                b = D.*[m[1]; 2*m[2]; 3*m[3]]
                c = A\b
                V_mom = shift_potential(equilibrium, c)
                lines!(ax22, domain, [V_mom(x) for x in domain], color = :gold, linewidth = 2.0)

                #--------------------#
                #         MLE        #
                #--------------------#
                
                # Export figures
                savefig("../../res/fig/mle/good/$n.png", fig1)
                savefig("../../res/fig/mle/bad/$n.png", fig2)
        end
end

# Execute the functions 
output = generate_samples()
main(output.equilibrium, output.n_samples)
