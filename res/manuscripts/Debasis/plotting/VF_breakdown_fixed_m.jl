include("../../../../inc/IO.jl")
include("../../../../inc/PlottingTools.jl")

###############################
#          PARAMETERS         # 
###############################

a = 0.0001::Float64
B = 1.0000::Float64
p = 0.0000::Float64
ρ = 0.5520::Float64
q = 0.0000::Float64
α = 0.6500::Float64
Θ = 0.8790::Float64
τ = 0.2000::Float64
ψ = 1.0000::Float64
g = 0.6000::Float64
L0 = 0.0673::Float64
Iy = 0.1000::Float64
γm = 0.8500::Float64
γw = 0.7050::Float64
hm0 = 0.2574::Float64
hw0 = 0.2455::Float64
zw0 = 1000.0000::Float64
zt0 = 1.0000::Float64

##############################
#          FUNCTIONS         # 
##############################

# Family potential lifetime income (Eq. (9) in the list)
yt(σI, mT) = ((zt0*zw0)^(1 - α))*((hw0^α) + (mT^(1 - α))*(hm0^α))

# Male consumption (Eq. (20) in the list)
cm(σI, mT) = yt(σI, mT)/(1.0::Float64 + (σI/mT) + ρ)

# Adult-age human capital (Eq. (11) in the list)
hm(m, w, σI, mT) = B*((m + Θ*τ)/(m + Θ*τ + g))
hw(m, w, σI, mT) = B*((w + Θ*τ)/(w + Θ*τ + g))

# System of equations ((24)~(25) in the list) to determine the parent's choice for child's education
function f(e, σI, mT)
        [(e[1] + Θ*τ)*(e[1] + Θ*τ + g) - (1.0::Float64 + q)*(γm/γw)*(e[2] + Θ*τ)*(e[2] + Θ*τ + g),
         (e[2] + Θ*τ)*(e[2] + Θ*τ + g) - ((2.0::Float64*γw*g)/(ρ*(1.0::Float64 + q)))*(0.5::Float64*(e[1] + (1.0::Float64 + q)*e[2]) + a/yt(σI, mT) + τ/(2.0::Float64 + Iy))]
end

# Quantity-quality tradeoff (Eq. (22) in the list)
ot(m, w, σI, mT) = (ρ/(1.0::Float64 + (σI/mT) + ρ))/(a/yt(σI, mT) + 0.5*(m + (1.0::Float64 + q)*w) + τ/(2.0::Float64 + Iy))

# Utility function (Eq. (87) in the ver. 22 of the paper)
V(m, w, σI, mT) = (1.0::Float64 + σI/mT)*log(cm(σI, mT)) + (σI/mT)*(log(σI/mT)) + ρ*log(ot(m, w, σI, mT)) + (γm*log(hm(m, w, σI, mT)) + γw*log(hw(m, w, σI, mT)))

#############################
#          DOMAINS          # 
#############################

# Number of points in the mT domain 
NσI = convert(Int64,3e2)

# mT domain
σI_inf = 0.00001::Float64 
σI_sup = 1.0000::Float64 
domain_σI = LinRange(σI_inf, σI_sup, NσI)

domain_mT = [1.0, 1.5, 2.0, 3.0]
NmT = length(domain_mT) 

###########################
#          PLOTS          # 
###########################

# Axis for VF 
fig, ax = mkfig(size = [1800,1200],
                bg_out = :white,
                limits = ((domain_σI[1], domain_σI[end]), nothing),
                lab = [L"\mathbf{\sigma_I}", L"\mathbf{V_F}"],
                lab_pad = [-40.0,-40.0],
                x_ticks = [0.0001, 0.2, 0.4, 0.6, 0.8, 1.0],
                y_ticks = [0.4,1.0],
)

# Axis for first term 
#=
fig, ax1 = mkfig(fig = fig,
                 box_position = [2,1],
                 limits = ((domain_σI[1], domain_σI[end]), nothing),
                 lab = [L"\mathbf{\sigma_I}", L"\textbf{consumption}"],
                 toggle_lab = [false, true],
                 lab_pad = [-40.0,-40.0],
                 #x_ticks = [0.1, 1.0, 2.0, 3.0],
                 #y_ticks = [-1.5, 0.0],
                 toggle_ticks_lab = [false, true]
)

# Axis for second term 
fig, ax2 = mkfig(fig = fig,
                 box_position = [3,1],
                 limits = ((domain_σI[1], domain_σI[end]), nothing),
                 lab = [L"\mathbf{\sigma_I}", L"\textbf{human capital}"],
                 toggle_lab = [false, true],
                 lab_pad = [-40.0,-60.0],
                 #x_ticks = [0.1, 1.0, 2.0, 3.0],
                 #y_ticks = [-0.76428,-0.76425],
                 toggle_ticks_lab = [false, true],
                 ticks_lab_trunc = [1,5]
)

# Axis for third term 
fig, ax3 = mkfig(fig = fig,
                 box_position = [4,1],
                 limits = ((domain_σI[1], domain_σI[end]), nothing),
                 lab = [L"\mathbf{\sigma_I}", L"\textbf{fertility}"],
                 toggle_lab = [true, true],
                 lab_pad = [-40.0,-40.0],
                 #x_ticks = [0.1, 1.0, 2.0, 3.0],
                 #y_ticks = [-1.5,-0.5],
                 toggle_ticks_lab = [true, true]
)
=#

# Loop over the σI domain
using NonlinearSolve 
@showprogress for j in 1:NmT
        # Initialise the vector for the utility function
        VF = Vector{Float64}(undef, NσI)
        # Initialise vectors for its terms
        consumption = Vector{Float64}(undef, NσI)
        human_capital = Vector{Float64}(undef, NσI)
        fertility = Vector{Float64}(undef, NσI)
        
        # Loop over the mT domain
        for k in 1:NσI
                # Define the function associated to the system
                F(e, p) = f(e, domain_σI[k], domain_mT[j])
                # Solve numerically the system (24)~(25)
                problem = NonlinearProblem(F, [0.75, 0.75], 0.0) 
                solution = solve(problem, abstol=1e-12, reltol=1e12)
                em = (solution.u)[1]
                ew = (solution.u)[2]

                # Compute the value of the utility function at the current mT
                VF[k] = V(em, ew, domain_σI[k], domain_mT[j])
                # Compute the values of its 3 terms at the current mT
                consumption[k] = (1.0::Float64 + domain_σI[k]/domain_mT[j])*log(cm(domain_σI[k], domain_mT[j])) + (domain_σI[k]/domain_mT[j])*(log(domain_σI[k]/domain_mT[j]))
                human_capital[k] = (γm*log(hm(em, ew, domain_σI[k], domain_mT[j])) + γw*log(hw(em, ew, domain_σI[k], domain_mT[j])))
                fertility[k] = ρ*log(ot(em, ew, domain_σI[k], domain_mT[j]))
        end

        if j==NmT
                display(VF[1])
                display(VF[end])
        end

        # Plot the utility function 
        lines!(ax, domain_σI, VF, color = j, colormap = :darktest, colorrange = (1,NmT), linewidth = 4.5)
        # Plot the total consumption
        #lines!(ax1, domain_σI, consumption, color = j, colormap = :Spectral_5, colorrange = (1,NmT), linewidth = 4.5)
        # Plot the human capital
        #lines!(ax2, domain_σI, human_capital, color = j, colormap = :Spectral_5, colorrange = (1,NmT), linewidth = 4.5)
        # Plot the fertility
        #lines!(ax3, domain_σI, fertility, color = j, colormap = :Spectral_5, colorrange = (1,NmT), linewidth = 4.5)
end

# Export the figure
save("../fig/fig4.png", fig)
