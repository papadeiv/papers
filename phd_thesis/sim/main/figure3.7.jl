using DynamicalSystems, DifferentialEquations
using LaTeXStrings, CairoMakie, Makie.Colors
using Statistics, SpecialFunctions, ProgressMeter

#                  Arnold and Boxler process 
# (stochastic transcritical bifurcation with multiplicative noise)

# Define the parameters of the process
σ = sqrt(0.8)

# Define initial state
x0 = [0.1]

# Define temporal evolution quantities
T = 10.00
δt = 1e-3

# Define the deterministic dynamics
function iip_det!(f, x, y, t)
        f[1] = y[1]*x[1] - (x[1])^2 + (0.5)*(σ^2)*x[1]
        return nothing
end

# Define the stochastic dynamics
function iip_stoc!(f, x, y, t)
        f[1] = +σ*x[1]
        return nothing
end

# Assemble the SDS
normal_form = SDEProblem(iip_det!, iip_stoc!, x0, (0.0, T), 0.0)

# Define parameter values
Ne = 1000
Ny = 100
counter = 1
y_values = LinRange(0.1,3.0,Ny)

# Array for ensemble paths' endpoints 
Xt = fill(Float64[], 0) 
# Arrays for storing the ensemble's entire simple sample paths for a specific value of y
St = fill(Float64[], 0)
t = Float64[]

println("Simulating the ensemble sample paths")
# Loop over the ensemble sample paths
@showprogress for y in y_values
        local xt = Float64[]
        for j in 1:Ne
                local normal_form = SDEProblem(iip_det!, iip_stoc!, x0, (0.0, T), y)
                local sol = solve(normal_form, EM(), dt=δt)
                # Store the ensemble path endpoint 
                push!(xt, sol[1,end])
                # Store the y-valued entire sample path
                if counter == 34 
                        push!(St, sol[1,:])
                        if j == 1
                                global t = sol.t
                        end
                end
        end
        push!(Xt, xt)
        global counter += 1
end

# Post processing to get the ensemble endpoint matrix
Et = Array{Float64, 2}(undef, Ne, Ny)
for j in 1:Ny
        ensemble_at_y = Xt[j]
        for k in 1:Ne
                Et[k,j] = ensemble_at_y[k] 
        end
end
 
# Plot the ensemble's endpoints parametric paths
CairoMakie.activate!(; px_per_unit = 3)
fig1 = Figure(; size = (1200, 600), backgroundcolor = :transparent)
ax1 = Axis(fig1[1,1],
    # Background
    backgroundcolor = :transparent,
    xgridvisible = false,
    ygridvisible = false,
    limits = ((0.1,3), (0,3)),
    # Title
    #title = L"r = r_c = %$R",
    titlevisible = false,
    titlesize = 22,
    titlealign = :center,
    titlegap = 4.0,
    # x-axis
    xlabel = L"\mu",
    xlabelvisible = true,
    xlabelsize = 24,
    xlabelcolor = :black,
    xlabelpadding = -20.0,
    xticks = [0,1,2,3],
    xticksvisible = true,
    xticksize = 10,
    xticklabelsvisible = true,
    xticklabelsize = 24,
    xtickformat = "{:.0f}",
    xscale = identity, #log10,
    xaxisposition = :bottom,
    # y-axis
    ylabel = L"moments",
    ylabelvisible = true,
    ylabelsize = 24,
    ylabelcolor = :black,
    ylabelpadding = 0,
    yticks = [0,1,2,3],
    yticksvisible = true,
    yticksize = 10,
    yticklabelsvisible = true,
    yticklabelsize = 24,
    ytickformat = "{:.0f}",
    yscale = identity,
    yaxisposition = :left,
)

# Plot the ensemble's average parametric path
At = Vector{Float64}(undef, Ny)
for j in 1:Ny
        At[j] = mean(Et[:,j])
end
lines!(ax1, y_values, At, color = (:blue, 0.5), linewidth=2, label = L"\mu_h(x)")

# Plot the ensemble's variance at specific values of y
Vt = Vector{Float64}(undef, Ny)
for j in 1:Ny
        Vt[j] = var(Et[:,j]; mean=At[j]) 
end
lines!(ax1, y_values, Vt, color = (:red, 0.5), linewidth=2, label = L"v_h(x)")

# Plot the theoretical mean from the closed formula given in (33)
fig2 = Figure(; size = (1200, 600), backgroundcolor = :transparent)
ax2 = Axis(fig2[1,1],
    # Background
    backgroundcolor = :transparent,
    xgridvisible = false,
    ygridvisible = false,
    limits = ((y_values[1],y_values[end]), (-0.25,8)),
    # Title
    #title = L"r = r_c = %$R",
    titlevisible = false,
    titlesize = 22,
    titlealign = :center,
    titlegap = 4.0,
    # x-axis
    xlabel = L"\mu",
    xlabelvisible = true,
    xlabelsize = 20,
    xlabelcolor = :black,
    xlabelpadding = -20.0,
    xticks = [0,1,2,3],
    xticksvisible = true,
    xticksize = 6,
    xticklabelsvisible = true,
    xticklabelsize = 14,
    xtickformat = "{:.0f}",
    xscale = identity, #log10,
    xaxisposition = :bottom,
    # y-axis
    ylabel = L"moments",
    ylabelvisible = true,
    ylabelsize = 20,
    ylabelcolor = :black,
    ylabelpadding = 10,
    yticks = [0,1,3,4],
    yticksvisible = true,
    yticksize = 6,
    yticklabelsvisible = true,
    yticklabelsize = 14,
    ytickformat = "{:.2f}",
    yscale = identity,
    yaxisposition = :left,
)
y_range = [0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.4,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.48,0.49,0.5,0.51,0.52,0.53,0.54,0.55,0.56,0.57,0.58,0.59,0.6,0.61,0.62,0.63,0.64,0.65,0.66,0.67,0.68,0.69,0.7,0.71,0.72,0.73,0.74,0.75,0.76,0.77,0.78,0.79,0.8,0.81,0.82,0.83,0.84,0.85,0.86,0.87,0.88,0.89,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,1.,1.01,1.02,1.03,1.04,1.05,1.06,1.07,1.08,1.09,1.1,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.2,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,1.3,1.31,1.32,1.33,1.34,1.35,1.36,1.37,1.38,1.39,1.4,1.41,1.42,1.43,1.44,1.45,1.46,1.47,1.48,1.49,1.5,1.51,1.52,1.53,1.54,1.55,1.56,1.57,1.58,1.59,1.6,1.61,1.62,1.63,1.64,1.65,1.66,1.67,1.68,1.69,1.7,1.71,1.72,1.73,1.74,1.75,1.76,1.77,1.78,1.79,1.8,1.81,1.82,1.83,1.84,1.85,1.86,1.87,1.88,1.89,1.9,1.91,1.92,1.93,1.94,1.95,1.96,1.97,1.98,1.99,2.,2.01,2.02,2.03,2.04,2.05,2.06,2.07,2.08,2.09,2.1,2.11,2.12,2.13,2.14,2.15,2.16,2.17,2.18,2.19,2.2,2.21,2.22,2.23,2.24,2.25,2.26,2.27,2.28,2.29,2.3,2.31,2.32,2.33,2.34,2.35,2.36,2.37,2.38,2.39,2.4,2.41,2.42,2.43,2.44,2.45,2.46,2.47,2.48,2.49,2.5,2.51,2.52,2.53,2.54,2.55,2.56,2.57,2.58,2.59,2.6,2.61,2.62,2.63,2.64,2.65,2.66,2.67,2.68,2.69,2.7,2.71,2.72,2.73,2.74,2.75,2.76,2.77,2.78,2.79,2.8,2.81,2.82,2.83,2.84,2.85,2.86,2.87,2.88,2.89,2.9,2.91,2.92,2.93,2.94,2.95,2.96,2.97,2.98,2.99]
theoretical_M = [0.288334,0.28031,0.272708,0.265501,0.258662,0.252169,0.246001,0.240137,0.23456,0.229253,0.2242,0.219386,0.214798,0.210424,0.206252,0.202271,0.198471,0.194843,0.191379,0.188068,0.184906,0.181883,0.178993,0.176231,0.17359,0.171064,0.168649,0.166339,0.164131,0.162019,0.16,0.15807,0.156225,0.154462,0.152777,0.151169,0.149633,0.148167,0.146769,0.145437,0.144167,0.142958,0.141808,0.140715,0.139677,0.138693,0.137761,0.136878,0.136045,0.135259,0.13452,0.133825,0.133175,0.132567,0.132001,0.131476,0.130991,0.130545,0.130137,0.129767,0.129434,0.129137,0.128875,0.128648,0.128456,0.128298,0.128173,0.128081,0.128022,0.127995,0.128,0.128037,0.128104,0.128203,0.128333,0.128494,0.128684,0.128906,0.129157,0.129439,0.12975,0.130092,0.130464,0.130865,0.131297,0.131758,0.13225,0.132772,0.133324,0.133907,0.13452,0.135164,0.135838,0.136544,0.137281,0.13805,0.138851,0.139683,0.140548,0.141446,0.142377,0.143342,0.14434,0.145373,0.14644,0.147543,0.148681,0.149855,0.151066,0.152314,0.1536,0.154924,0.156287,0.15769,0.159133,0.160617,0.162142,0.16371,0.165321,0.166976,0.168675,0.170421,0.172212,0.174051,0.175938,0.177874,0.17986,0.181898,0.183987,0.18613,0.188328,0.190581,0.19289,0.195258,0.197685,0.200173,0.202722,0.205334,0.208012,0.210755,0.213566,0.216446,0.219397,0.22242,0.225518,0.228691,0.231942,0.235272,0.238684,0.242179,0.24576,0.249428,0.253186,0.257035,0.260978,0.265018,0.269156,0.273396,0.27774,0.282189,0.286748,0.291419,0.296205,0.301108,0.306132,0.311279,0.316554,0.321959,0.327497,0.333173,0.33899,0.344951,0.351061,0.357322,0.36374,0.370319,0.377063,0.383975,0.391062,0.398327,0.405775,0.413412,0.421242,0.429271,0.437505,0.445948,0.454606,0.463487,0.472595,0.481937,0.49152,0.50135,0.511435,0.521781,0.532395,0.543287,0.554462,0.56593,0.577698,0.589776,0.602171,0.614894,0.627954,0.64136,0.655122,0.66925,0.683756,0.698651,0.713944,0.729649,0.745778,0.762342,0.779354,0.796829,0.814779,0.833218,0.852161,0.871624,0.891621,0.912169,0.933283,0.954982,0.977282,1.0002,1.02376,1.04798,1.07287,1.09846,1.12478,1.15183,1.17965,1.20825,1.23767,1.26793,1.29904,1.33105,1.36398,1.39785,1.43269,1.46854,1.50543,1.54338,1.58244,1.62264,1.66401,1.70659,1.75042,1.79553,1.84198,1.88979,1.93902,1.98971,2.04191,2.09566,2.15102,2.20803,2.26675,2.32724,2.38954,2.45373,2.51986,2.588,2.65821,2.73055,2.8051,2.88194,2.96112,3.04274,3.12688,3.21361,3.30301,3.39519,3.49024,3.58823,3.68929,3.7935,3.90097,4.01182,4.12615,4.24409,4.36574,4.49125,4.62074,4.75433,4.89219,5.03444,5.18123,5.33273,5.48909,5.65048]
lines!(ax2, y_range, theoretical_M, color = :indigo, linewidth=2, label = L"m_s^y(x)")

# Plot the theoretical mean as computed from the closed formila of the F-P distribution given in (32) and your own normalization constant
quadrature_M = [0.177468,0.183342,0.189372,0.195556,0.201894,0.208382,0.215018,0.221801,0.228726,0.235792,0.242995,0.250331,0.257798,0.265392,0.273108,0.280944,0.288896,0.296959,0.30513,0.313404,0.321778,0.330248,0.338809,0.347458,0.35619,0.365002,0.37389,0.38285,0.391878,0.400971,0.410126,0.419339,0.428606,0.437925,0.447293,0.456707,0.466163,0.475659,0.485193,0.494762,0.504364,0.513997,0.523658,0.533346,0.543059,0.552794,0.562551,0.572328,0.582123,0.591935,0.601763,0.611605,0.62146,0.631328,0.641207,0.651097,0.660996,0.670904,0.68082,0.690744,0.700674,0.710611,0.720553,0.7305,0.740453,0.750409,0.76037,0.770334,0.780302,0.790273,0.800246,0.810222,0.8202,0.83018,0.840162,0.850146,0.860132,0.870119,0.880107,0.890096,0.900086,0.910078,0.92007,0.930063,0.940056,0.95005,0.960045,0.970041,0.980036,0.990033,1.00003,1.01003,1.02002,1.03002,1.04002,1.05002,1.06002,1.07001,1.08001,1.09001,1.10001,1.11001,1.12001,1.13001,1.14001,1.15001,1.16,1.17,1.18,1.19,1.2,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,1.3,1.31,1.32,1.33,1.34,1.35,1.36,1.37,1.38,1.39,1.4,1.41,1.42,1.43,1.44,1.45,1.46,1.47,1.48,1.49,1.5,1.51,1.52,1.53,1.54,1.55,1.56,1.57,1.58,1.59,1.6,1.61,1.62,1.63,1.64,1.65,1.66,1.67,1.68,1.69,1.7,1.71,1.72,1.73,1.74,1.75,1.76,1.77,1.78,1.79,1.8,1.81,1.82,1.83,1.84,1.85,1.86,1.87,1.88,1.89,1.9,1.91,1.92,1.93,1.94,1.95,1.96,1.97,1.98,1.99,2.,2.01,2.02,2.03,2.04,2.05,2.06,2.07,2.08,2.09,2.1,2.11,2.12,2.13,2.14,2.15,2.16,2.17,2.18,2.19,2.2,2.21,2.22,2.23,2.24,2.25,2.26,2.27,2.28,2.29,2.3,2.31,2.32,2.33,2.34,2.35,2.36,2.37,2.38,2.39,2.4,2.41,2.42,2.43,2.44,2.45,2.46,2.47,2.48,2.49,2.5,2.51,2.52,2.53,2.54,2.55,2.56,2.57,2.58,2.59,2.6,2.61,2.62,2.63,2.64,2.65,2.66,2.67,2.68,2.69,2.7,2.71,2.72,2.73,2.74,2.75,2.76,2.77,2.78,2.79,2.8,2.81,2.82,2.83,2.84,2.85,2.86,2.87,2.88,2.89,2.9,2.91,2.92,2.93,2.94,2.95,2.96,2.97,2.98,2.99]
lines!(ax2, y_range, quadrature_M, color = :red, linewidth=2, label = L"\int_0^{\infty}x\,p_s^y(x)dx\,,\;N_y=\int_0^{\infty}x^{\frac{2y}{\sigma^2}-1}e^{-\frac{2x}{\sigma^2}}")
lines!(ax1, y_range, quadrature_M, color = :black, linewidth=3, linestyle = :dash)


# Plot the same as above but with the normalization constant provided by Arnold and Boxler
quadrature_M_AB = [0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.4,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.48,0.49,0.5,0.51,0.52,0.53,0.54,0.55,0.56,0.57,0.58,0.59,0.6,0.61,0.62,0.63,0.64,0.65,0.66,0.67,0.68,0.69,0.7,0.71,0.72,0.73,0.74,0.75,0.76,0.77,0.78,0.79,0.8,0.81,0.82,0.83,0.84,0.85,0.86,0.87,0.88,0.89,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,1.,1.01,1.02,1.03,1.04,1.05,1.06,1.07,1.08,1.09,1.1,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.2,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,1.3,1.31,1.32,1.33,1.34,1.35,1.36,1.37,1.38,1.39,1.4,1.41,1.42,1.43,1.44,1.45,1.46,1.47,1.48,1.49,1.5,1.51,1.52,1.53,1.54,1.55,1.56,1.57,1.58,1.59,1.6,1.61,1.62,1.63,1.64,1.65,1.66,1.67,1.68,1.69,1.7,1.71,1.72,1.73,1.74,1.75,1.76,1.77,1.78,1.79,1.8,1.81,1.82,1.83,1.84,1.85,1.86,1.87,1.88,1.89,1.9,1.91,1.92,1.93,1.94,1.95,1.96,1.97,1.98,1.99,2.,2.01,2.02,2.03,2.04,2.05,2.06,2.07,2.08,2.09,2.1,2.11,2.12,2.13,2.14,2.15,2.16,2.17,2.18,2.19,2.2,2.21,2.22,2.23,2.24,2.25,2.26,2.27,2.28,2.29,2.3,2.31,2.32,2.33,2.34,2.35,2.36,2.37,2.38,2.39,2.4,2.41,2.42,2.43,2.44,2.45,2.46,2.47,2.48,2.49,2.5,2.51,2.52,2.53,2.54,2.55,2.56,2.57,2.58,2.59,2.6,2.61,2.62,2.63,2.64,2.65,2.66,2.67,2.68,2.69,2.7,2.71,2.72,2.73,2.74,2.75,2.76,2.77,2.78,2.79,2.8,2.81,2.82,2.83,2.84,2.85,2.86,2.87,2.88,2.89,2.9,2.91,2.92,2.93,2.94,2.95,2.96,2.97,2.98,2.99]
lines!(ax2, y_range, quadrature_M_AB, color = :black, linewidth=2, label = L"\int_0^{\infty}xp_s^y(x)dx\,,\;N_y=\Gamma(\frac{2y}{\sigma^2})(\frac{\sigma^2}{2})^2")

# Plot the theoretical variance from the closed formula given in (33)
theoretical_V = [0.0754697,0.0730056,0.0713198,0.0703604,0.0700807,0.0704385,0.0713962,0.0729193,0.0749768,0.0775408,0.0805856,0.0840881,0.0880271,0.0923832,0.0971389,0.102278,0.107786,0.113649,0.119854,0.12639,0.133247,0.140414,0.147883,0.155645,0.163692,0.172018,0.180615,0.189478,0.198599,0.207975,0.2176,0.227469,0.237577,0.247921,0.258497,0.2693,0.280328,0.291576,0.303043,0.314724,0.326617,0.33872,0.351029,0.363543,0.376258,0.389173,0.402286,0.415594,0.429096,0.442789,0.456672,0.470742,0.484999,0.499439,0.514063,0.528867,0.54385,0.559012,0.574349,0.589861,0.605546,0.621402,0.637429,0.653624,0.669986,0.686513,0.703205,0.72006,0.737075,0.754251,0.771584,0.789074,0.80672,0.824519,0.84247,0.860572,0.878822,0.897221,0.915765,0.934453,0.953285,0.972256,0.991368,1.01062,1.03,1.04952,1.06917,1.08895,1.10886,1.1289,1.14906,1.16934,1.18974,1.21026,1.2309,1.25165,1.27252,1.29349,1.31457,1.33575,1.35704,1.37843,1.39991,1.42149,1.44316,1.46492,1.48677,1.5087,1.53071,1.55279,1.57495,1.59718,1.61948,1.64185,1.66427,1.68676,1.70929,1.73188,1.75451,1.77718,1.7999,1.82264,1.84542,1.86822,1.89104,1.91388,1.93673,1.95959,1.98245,2.0053,2.02815,2.05098,2.0738,2.09659,2.11935,2.14207,2.16475,2.18738,2.20995,2.23247,2.25491,2.27728,2.29957,2.32176,2.34386,2.36586,2.38774,2.4095,2.43113,2.45262,2.47397,2.49516,2.51618,2.53703,2.5577,2.57818,2.59845,2.6185,2.63833,2.65793,2.67728,2.69637,2.71519,2.73373,2.75198,2.76992,2.78754,2.80482,2.82176,2.83834,2.85455,2.87037,2.88578,2.90078,2.91534,2.92946,2.9431,2.95627,2.96894,2.98109,2.99271,3.00378,3.01427,3.02419,3.03349,3.04217,3.05021,3.05758,3.06427,3.07025,3.07551,3.08002,3.08377,3.08673,3.08887,3.09019,3.09064,3.09023,3.08891,3.08667,3.08349,3.07934,3.0742,3.06805,3.06086,3.05262,3.0433,3.03287,3.02132,3.00862,2.99476,2.97971,2.96346,2.94598,2.92726,2.90727,2.88601,2.86346,2.8396,2.81442,2.78791,2.76007,2.73089,2.70036,2.66849,2.63526,2.6007,2.5648,2.52759,2.48907,2.44926,2.40819,2.3659,2.32241,2.27778,2.23204,2.18527,2.13751,2.08886,2.03938,1.98917,1.93835,1.88701,1.8353,1.78336,1.73134,1.67943,1.6278,1.57668,1.52629,1.47689,1.42876,1.38219,1.33752,1.29511,1.25534,1.21865,1.18549,1.15636,1.13182,1.11245,1.09888,1.09182,1.092,1.10024,1.11741,1.14445,1.18239,1.23232,1.29544,1.37302,1.46645,1.57722,1.70692,1.85729,2.03019,2.22763,2.45175,2.70489,2.98955,3.3084,3.66435,4.0605,4.5002,4.98703,5.52487,6.11787,6.77049,7.48753,8.27414]
lines!(ax2, y_range, theoretical_V, color = :indigo, linewidth=2, linestyle = :dash, label = L"v_s^y(x)")

# Plot the theoretical variance as computed by quadrature from the closed formula of the F-P distribution given in (32) and your own normalization constant
quadrature_V = [0.0816297,0.0822988,0.0833417,0.0846868,0.0862817,0.0880877,0.0900753,0.0922213,0.0945075,0.0969191,0.0994438,0.102071,0.104793,0.107602,0.110491,0.113454,0.116486,0.119583,0.12274,0.125954,0.12922,0.132536,0.135898,0.139304,0.142751,0.146235,0.149756,0.153311,0.156897,0.160512,0.164156,0.167825,0.171518,0.175234,0.178971,0.182728,0.186457,0.190295,0.194104,0.197928,0.201765,0.205615,0.209477,0.21335,0.217233,0.221126,0.225027,0.228937,0.232854,0.236778,0.240708,0.244645,0.248586,0.252533,0.256484,0.26044,0.264399,0.268362,0.272329,0.276298,0.28027,0.284245,0.288221,0.2922,0.296181,0.300164,0.304148,0.308134,0.312121,0.316109,0.320098,0.324089,0.32808,0.332072,0.336065,0.340059,0.344053,0.348047,0.352043,0.356038,0.360035,0.364031,0.368028,0.372025,0.376022,0.38002,0.384018,0.388016,0.392015,0.396013,0.400012,0.40401,0.408009,0.412008,0.416008,0.420007,0.424006,0.428005,0.432005,0.436004,0.440004,0.444003,0.448003,0.452003,0.456002,0.460002,0.464002,0.468002,0.472002,0.476001,0.480001,0.484001,0.488001,0.492001,0.496001,0.500001,0.504001,0.508001,0.512,0.516,0.52,0.524,0.528,0.532,0.536,0.54,0.544,0.548,0.552,0.556,0.56,0.564,0.568,0.572,0.576,0.58,0.584,0.588,0.592,0.596,0.6,0.604,0.608,0.612,0.616,0.62,0.624,0.628,0.632,0.636,0.64,0.644,0.648,0.652,0.656,0.66,0.664,0.668,0.672,0.676,0.68,0.684,0.688,0.692,0.696,0.7,0.704,0.708,0.712,0.716,0.72,0.724,0.728,0.732,0.736,0.74,0.744,0.748,0.752,0.756,0.76,0.764,0.768,0.772,0.776,0.78,0.784,0.788,0.792,0.796,0.8,0.804,0.808,0.812,0.816,0.82,0.824,0.828,0.832,0.836,0.84,0.844,0.848,0.852,0.856,0.86,0.864,0.868,0.872,0.876,0.88,0.884,0.888,0.892,0.896,0.9,0.904,0.908,0.912,0.916,0.92,0.924,0.928,0.932,0.936,0.94,0.944,0.948,0.952,0.956,0.96,0.964,0.968,0.972,0.976,0.98,0.984,0.988,0.992,0.996,1.,1.004,1.008,1.012,1.016,1.02,1.024,1.028,1.032,1.036,1.04,1.044,1.048,1.052,1.056,1.06,1.064,1.068,1.072,1.076,1.08,1.084,1.088,1.092,1.096,1.1,1.104,1.108,1.112,1.116,1.12,1.124,1.128,1.132,1.136,1.14,1.144,1.148,1.152,1.156,1.16,1.164,1.168,1.172,1.176,1.18,1.184,1.188,1.192,1.196]
lines!(ax2, y_range, quadrature_V, color = :red, linewidth=2, linestyle = :dash, label = L"\int_0^{\infty}(x-\mu)^2 p_s^y(x)dx\,,\;N_y=\int_0^{\infty}x^{\frac{2y}{\sigma^2}-1}e^{-\frac{2x}{\sigma^2}}")
lines!(ax1, y_range, quadrature_V, color = :black, linewidth=2, linestyle = :dash)

# Plot the same as above but with the normalization constant provided by Arnold and Boxler
quadrature_V_AB = [0.0399986,0.0439993,0.0479997,0.0519998,0.0559999,0.06,0.064,0.068,0.072,0.076,0.08,0.084,0.088,0.092,0.096,0.1,0.104,0.108,0.112,0.116,0.12,0.124,0.128,0.132,0.136,0.14,0.144,0.148,0.152,0.156,0.16,0.164,0.168,0.172,0.176,0.18,0.183999,0.187999,0.192,0.196,0.2,0.204,0.208,0.212,0.216,0.22,0.224,0.228,0.232,0.236,0.24,0.244,0.248,0.252,0.256,0.26,0.264,0.268,0.272,0.276,0.28,0.284,0.288,0.292,0.296,0.3,0.304,0.308,0.312,0.316,0.32,0.324,0.328,0.332,0.336,0.34,0.344,0.348,0.352,0.356,0.36,0.364,0.368,0.372,0.376,0.38,0.384,0.388,0.392,0.396,0.4,0.404,0.408,0.412,0.416,0.42,0.424,0.428,0.432,0.436,0.44,0.444,0.448,0.452,0.456,0.46,0.464,0.468,0.472,0.476,0.48,0.484,0.488,0.492,0.496,0.5,0.504,0.508,0.512,0.516,0.52,0.524,0.528,0.532,0.536,0.54,0.544,0.548,0.552,0.556,0.56,0.564,0.568,0.572,0.576,0.58,0.584,0.588,0.592,0.596,0.6,0.604,0.608,0.612,0.616,0.62,0.624,0.628,0.632,0.636,0.64,0.644,0.648,0.652,0.656,0.66,0.664,0.668,0.672,0.676,0.68,0.684,0.688,0.692,0.696,0.7,0.704,0.708,0.712,0.716,0.72,0.724,0.728,0.732,0.736,0.74,0.744,0.748,0.752,0.756,0.76,0.764,0.768,0.772,0.776,0.78,0.784,0.788,0.792,0.796,0.8,0.804,0.808,0.812,0.816,0.82,0.824,0.828,0.832,0.836,0.84,0.844,0.848,0.852,0.856,0.86,0.864,0.868,0.872,0.876,0.88,0.884,0.888,0.892,0.896,0.9,0.904,0.908,0.912,0.916,0.92,0.924,0.928,0.932,0.936,0.94,0.944,0.948,0.952,0.956,0.96,0.964,0.968,0.972,0.976,0.98,0.984,0.988,0.992,0.996,1.,1.004,1.008,1.012,1.016,1.02,1.024,1.028,1.032,1.036,1.04,1.044,1.048,1.052,1.056,1.06,1.064,1.068,1.072,1.076,1.08,1.084,1.088,1.092,1.096,1.1,1.104,1.108,1.112,1.116,1.12,1.124,1.128,1.132,1.136,1.14,1.144,1.148,1.152,1.156,1.16,1.164,1.168,1.172,1.176,1.18,1.184,1.188,1.192,1.196]
lines!(ax2, y_range, quadrature_V_AB, color = :black, linewidth=2, linestyle = :dash, label = L"\int_0^{\infty}xp_s^y(x)dx\,,\;N_y=\Gamma(\frac{2y}{\sigma^2})(\frac{\sigma^2}{2})^2")
axislegend(ax2, position = :lt, labelsize = 20)

# Plot the critical manifold
y_pos = range(y_values[1], y_values[end], length=10)
stable = 0.0 .+ 1.0*y_pos
unstable = 0.0 .+ 0.0*y_pos
#lines!(ax1, y_pos, stable, color = :black, linewidth = 1.5)
#lines!(ax1, y_pos, unstable, color = :black, linewidth = 1.5, linestyle = :dash)

# Plot the ensemble' endpoints distribution at specific y-values and overlay the actual F-P distribution for the same values
fig3 = Figure(; size = (1200, 600), backgroundcolor = :transparent)
Y = y_values[11]
ax3 = Axis(fig3[1,1],
    # Background
    backgroundcolor = :transparent,
    xgridvisible = false,
    ygridvisible = false,
    limits = ((0,6), (0,1)),
    # Title
    title = (L"\mu=%$(round(Y; digits=3))"),
    titlevisible = true,
    titlesize = 22,
    titlealign = :center,
    titlegap = 4.0,
    # x-axis
    xlabel = L"x",
    xlabelvisible = false,
    xlabelsize = 20,
    xlabelcolor = :black,
    xlabelpadding = -20.0,
    xticks = [0,2,4,6],
    xticksvisible = true,
    xticksize = 10,
    xticklabelsvisible = false,
    xticklabelsize = 14,
    xtickformat = "{:.0f}",
    xscale = identity, #log10,
    xaxisposition = :bottom,
    # y-axis
    ylabel = L"density",
    ylabelvisible = true,
    ylabelsize = 24,
    ylabelcolor = :black,
    ylabelpadding = -10,
    yticks = [0,1],
    yticksvisible = true,
    yticksize = 10,
    yticklabelsvisible = true,
    yticklabelsize = 24,
    ytickformat = "{:.1f}",
    yscale = identity,
    yaxisposition = :left,
)
hist!(ax3, Et[:,11], bins = 100, normalization = :pdf, color = :red, strokecolor = :black, strokewidth = 1)
xrange = [0.01,0.11,0.21,0.31,0.41,0.51,0.61,0.71,0.81,0.91,1.01,1.11,1.21,1.31,1.41,1.51,1.61,1.71,1.81,1.91,2.01,2.11,2.21,2.31,2.41,2.51,2.61,2.71,2.81,2.91,3.01,3.11,3.21,3.31,3.41,3.51,3.61,3.71,3.81,3.91,4.01,4.11,4.21,4.31,4.41,4.51,4.61,4.71,4.81,4.91,5.01,5.11,5.21,5.31,5.41,5.51,5.61,5.71,5.81]
p = [2.64592,1.97502,1.52064,1.17614,0.911458,0.707107,0.548953,0.426378,0.33129,0.257478,0.200154,0.15562,0.121012,0.0941119,0.0731991,0.0569384,0.0442934,0.0344589,0.0268097,0.0208595,0.0162307,0.0126297,0.00982793,0.00764801,0.00595181,0.00463194,0.00360486,0.0028056,0.0021836,0.00169954,0.00132281,0.00102961,0.000801413,0.000623802,0.000485562,0.000377963,0.000294211,0.000229021,0.000178278,0.000138779,0.000108033,0.0000840996,0.000065469,0.0000509661,0.0000396763,0.0000308877,0.000024046,0.00001872,0.0000145737,0.0000113459,8.83302*10^-6,6.87676*10^-6,5.35379*10^-6,4.16813*10^-6,3.24507*10^-6,2.52645*10^-6,1.96697*10^-6,1.5314*10^-6,1.19229*10^-6]
lines!(ax3, xrange, p, color = :blue, linewidth = 2)


Y = y_values[21]
ax4 = Axis(fig3[1,2],
    # Background
    backgroundcolor = :transparent,
    xgridvisible = false,
    ygridvisible = false,
    limits = ((0,6), (0,1)),
    # Title
    title = (L"\mu=%$(round(Y; digits=3))"),
    titlevisible = true,
    titlesize = 22,
    titlealign = :center,
    titlegap = 4.0,
    # x-axis
    xlabel = L"x",
    xlabelvisible = false,
    xlabelsize = 20,
    xlabelcolor = :black,
    xlabelpadding = -20.0,
    xticks = [0,2,4,6],
    xticksvisible = true,
    xticksize = 10,
    xticklabelsvisible = false,
    xticklabelsize = 14,
    xtickformat = "{:.0f}",
    xscale = identity, #log10,
    xaxisposition = :bottom,
    # y-axis
    ylabel = L"density",
    ylabelvisible = false,
    ylabelsize = 20,
    ylabelcolor = :black,
    ylabelpadding = 10,
    yticks = [0,1],
    yticksvisible = true,
    yticksize = 10,
    yticklabelsvisible = false,
    yticklabelsize = 14,
    ytickformat = "{:.1f}",
    yscale = identity,
    yaxisposition = :left,
)
hist!(ax4, Et[:,21], bins = 100, normalization = :pdf, color = :red, strokecolor = :black, strokewidth = 1)
p = [0.191844,0.829048,1.02492,1.05437,1.00274,0.912751,0.807882,0.701275,0.600082,0.507885,0.426138,0.35504,0.294087,0.242407,0.198977,0.16274,0.132685,0.107883,0.0875017,0.0708162,0.0572002,0.0461203,0.037127,0.0298436,0.0239569,0.0192077,0.0153825,0.0123061,0.00983543,0.00785367,0.00626592,0.00499522,0.00397927,0.00316775,0.00252008,0.0020036,0.00159205,0.00126434,0.00100356,0.00079618,0.000631358,0.000500434,0.000396491,0.000314012,0.000248594,0.000196732,0.000155636,0.000123082,0.0000973066,0.000076905,0.0000607629,0.0000479953,0.0000379001,0.0000299204,0.0000236148,0.0000186335,0.0000146995,0.0000115935,9.14171*10^-6]
lines!(ax4, xrange, p, color = :blue, linewidth = 2)

Y = y_values[31]
ax5 = Axis(fig3[1,3],
    # Background
    backgroundcolor = :transparent,
    xgridvisible = false,
    ygridvisible = false,
    limits = ((0,6), (0,1)),
    # Title
    title = (L"\mu=%$(round(Y; digits=3))"),
    titlevisible = true,
    titlesize = 22,
    titlealign = :center,
    titlegap = 4.0,
    # x-axis
    xlabel = L"x",
    xlabelvisible = false,
    xlabelsize = 20,
    xlabelcolor = :black,
    xlabelpadding = -20.0,
    xticks = [0,2,4,6],
    xticksvisible = true,
    xticksize = 10,
    xticklabelsvisible = false,
    xticklabelsize = 14,
    xtickformat = "{:.0f}",
    xscale = identity, #log10,
    xaxisposition = :bottom,
    # y-axis
    ylabel = L"density",
    ylabelvisible = false,
    ylabelsize = 20,
    ylabelcolor = :black,
    ylabelpadding = 10,
    yticks = [0,1],
    yticksvisible = true,
    yticksize = 10,
    yticklabelsvisible = false,
    yticklabelsize = 14,
    ytickformat = "{:.1f}",
    yscale = identity,
    yaxisposition = :left,
)
hist!(ax5, Et[:,31], bins = 100, normalization = :pdf, color = :red, strokecolor = :black, strokewidth = 1)
p = [0.00914661,0.228839,0.454253,0.621534,0.725409,0.774748,0.78181,0.758443,0.714749,0.658768,0.596591,0.532636,0.469962,0.410571,0.355666,0.305863,0.261366,0.222097,0.187795,0.158089,0.132555,0.110748,0.0922272,0.0765764,0.0634093,0.0523754,0.0431622,0.0354942,0.0291309,0.0238646,0.019517,0.0159359,0.0129925,0.0105778,0.00860054,0.00698418,0.00566493,0.00458978,0.00371477,0.00300358,0.00242625,0.00195812,0.00157897,0.00127219,0.00102422,0.000823963,0.000662393,0.000532142,0.000427224,0.000342778,0.000274859,0.00022027,0.000176425,0.000141232,0.000113002,0.0000903692,0.0000722352,0.0000577136,0.0000460909]
lines!(ax5, xrange, p, color = :blue, linewidth = 2)

Y = y_values[41]
ax6 = Axis(fig3[1,4],
    # Background
    backgroundcolor = :transparent,
    xgridvisible = false,
    ygridvisible = false,
    limits = ((0,6), (0,1)),
    # Title
    title = (L"\mu=%$(round(Y; digits=3))"),
    titlevisible = true,
    titlesize = 22,
    titlealign = :center,
    titlegap = 4.0,
    # x-axis
    xlabel = L"x",
    xlabelvisible = false,
    xlabelsize = 20,
    xlabelcolor = :black,
    xlabelpadding = -20.0,
    xticks = [0,2,4,6],
    xticksvisible = true,
    xticksize = 6,
    xticklabelsvisible = false,
    xticklabelsize = 14,
    xtickformat = "{:.0f}",
    xscale = identity, #log10,
    xaxisposition = :bottom,
    # y-axis
    ylabel = L"density",
    ylabelvisible = false,
    ylabelsize = 20,
    ylabelcolor = :black,
    ylabelpadding = 10,
    yticks = [0,1],
    yticksvisible = true,
    yticksize = 10,
    yticklabelsvisible = false,
    yticklabelsize = 14,
    ytickformat = "{:.1f}",
    yscale = identity,
    yaxisposition = :left,
)
hist!(ax6, Et[:,41], bins = 100, normalization = :pdf, color = :red, strokecolor = :black, strokewidth = 1)
p = [0.000332981,0.048104,0.153213,0.278705,0.399071,0.499963,0.575093,0.623402,0.646911,0.649219,0.63452,0.606989,0.570436,0.52814,0.482796,0.436524,0.390926,0.347155,0.305994,0.267923,0.23319,0.201868,0.173899,0.149137,0.12738,0.10839,0.0919117,0.0776898,0.065474,0.0550268,0.0461277,0.0385749,0.0321862,0.0267988,0.0222689,0.01847,0.0152921,0.0126398,0.010431,0.00859533,0.00707259,0.00581173,0.0047695,0.00390936,0.00320059,0.00261738,0.00213816,0.00174488,0.00142255,0.00115867,0.000942882,0.000766618,0.000622786,0.000505534,0.000410039,0.000332336,0.000269165,0.000217851,0.000176201]
lines!(ax6, xrange, p, color = :blue, linewidth = 2)

Y = y_values[51]
ax7 = Axis(fig3[2,1],
    # Background
    backgroundcolor = :transparent,
    xgridvisible = false,
    ygridvisible = false,
    limits = ((0,6), (0,1)),
    # Title
    title = (L"\mu=%$(round(Y; digits=3))"),
    titlevisible = true,
    titlesize = 22,
    titlealign = :center,
    titlegap = 4.0,
    # x-axis
    xlabel = L"x",
    xlabelvisible = true,
    xlabelsize = 24,
    xlabelcolor = :black,
    xlabelpadding = -20.0,
    xticks = [0,2,4,6],
    xticksvisible = true,
    xticksize = 10,
    xticklabelsvisible = true,
    xticklabelsize = 24,
    xtickformat = "{:.0f}",
    xscale = identity, #log10,
    xaxisposition = :bottom,
    # y-axis
    ylabel = L"density",
    ylabelvisible = true,
    ylabelsize = 24,
    ylabelcolor = :black,
    ylabelpadding = -10,
    yticks = [0,1],
    yticksvisible = true,
    yticksize = 10,
    yticklabelsvisible = true,
    yticklabelsize = 24,
    ytickformat = "{:.1f}",
    yscale = identity,
    yaxisposition = :left,
)
hist!(ax7, Et[:,51], bins = 100, normalization = :pdf, color = :red, strokecolor = :black, strokewidth = 1)
p = [9.82002*10^-6,0.00823489,0.0421443,0.102009,0.179308,0.263636,0.345809,0.419006,0.478926,0.523473,0.552279,0.566193,0.566848,0.556291,0.53672,0.51029,0.478992,0.444582,0.408549,0.372109,0.336222,0.301611,0.268797,0.238129,0.20981,0.183934,0.160505,0.139463,0.120699,0.104075,0.089433,0.076604,0.0654179,0.0557075,0.0473126,0.0400822,0.0338766,0.0285678,0.0240401,0.0201894,0.0169232,0.0141598,0.0118271,0.00986258,0.00821146,0.00682653,0.00566705,0.00469807,0.00388966,0.00321631,0.00265631,0.00219127,0.00180563,0.00148626,0.00122212,0.00100392,0.000823891,0.000675519,0.000553372]
lines!(ax7, xrange, p, color = :blue, linewidth = 2)

Y = y_values[61]
ax8 = Axis(fig3[2,2],
    # Background
    backgroundcolor = :transparent,
    xgridvisible = false,
    ygridvisible = false,
    limits = ((0,6), (0,1)),
    # Title
    title = (L"\mu=%$(round(Y; digits=3))"),
    titlevisible = true,
    titlesize = 22,
    titlealign = :center,
    titlegap = 4.0,
    # x-axis
    xlabel = L"x",
    xlabelvisible = true,
    xlabelsize = 24,
    xlabelcolor = :black,
    xlabelpadding = -20.0,
    xticks = [0,2,4,6],
    xticksvisible = true,
    xticksize = 10,
    xticklabelsvisible = true,
    xticklabelsize = 24,
    xtickformat = "{:.0f}",
    xscale = identity, #log10,
    xaxisposition = :bottom,
    # y-axis
    ylabel = L"density",
    ylabelvisible = false,
    ylabelsize = 24,
    ylabelcolor = :black,
    ylabelpadding = 10,
    yticks = [0,1],
    yticksvisible = true,
    yticksize = 10,
    yticklabelsvisible = false,
    yticklabelsize = 24,
    ytickformat = "{:.1f}",
    yscale = identity,
    yaxisposition = :left,
)
hist!(ax8, Et[:,61], bins = 100, normalization = :pdf, color = :red, strokecolor = :black, strokewidth = 1)
p = [2.48733*10^-7,0.00120758,0.00992326,0.0319465,0.0689132,0.118884,0.177787,0.240749,0.303055,0.360721,0.410766,0.451261,0.481242,0.500558,0.509678,0.509516,0.501261,0.486242,0.465822,0.441316,0.413938,0.384768,0.354735,0.324612,0.295025,0.266455,0.239263,0.2137,0.189922,0.168012,0.147991,0.129832,0.113473,0.0988254,0.0857823,0.0742276,0.0640395,0.0550954,0.0472752,0.0404632,0.0345503,0.0294346,0.0250224,0.0212278,0.0179734,0.0151894,0.0128137,0.010791,0.00907268,0.00761597,0.0063835,0.00534271,0.00446538,0.00372711,0.00310688,0.00258664,0.00215093,0.00178654,0.00148223]
lines!(ax8, xrange, p, color = :blue, linewidth = 2)

Y = y_values[71]
ax9 = Axis(fig3[2,3],
    # Background
    backgroundcolor = :transparent,
    xgridvisible = false,
    ygridvisible = false,
    limits = ((0,6), (0,1)),
    # Title
    title = (L"\mu=%$(round(Y; digits=3))"),
    titlevisible = true,
    titlesize = 22,
    titlealign = :center,
    titlegap = 4.0,
    # x-axis
    xlabel = L"x",
    xlabelvisible = true,
    xlabelsize = 24,
    xlabelcolor = :black,
    xlabelpadding = -20.0,
    xticks = [0,2,4,6],
    xticksvisible = true,
    xticksize = 10,
    xticklabelsvisible = true,
    xticklabelsize = 24,
    xtickformat = "{:.0f}",
    xscale = identity, #log10,
    xaxisposition = :bottom,
    # y-axis
    ylabel = L"density",
    ylabelvisible = false,
    ylabelsize = 24,
    ylabelcolor = :black,
    ylabelpadding = 10,
    yticks = [0,1],
    yticksvisible = true,
    yticksize = 10,
    yticklabelsvisible = false,
    yticklabelsize = 24,
    ytickformat = "{:.1f}",
    yscale = identity,
    yaxisposition = :left,
)
hist!(ax9, Et[:,71], bins = 100, normalization = :pdf, color = :red, strokecolor = :black, strokewidth = 1)
p = [5.53514*10^-9,0.000155579,0.00205279,0.00878981,0.0232691,0.0470993,0.080304,0.12153,0.16848,0.218384,0.268413,0.315983,0.358951,0.395713,0.425225,0.446964,0.460864,0.467226,0.466627,0.459835,0.447732,0.431245,0.411297,0.38877,0.364472,0.339126,0.313356,0.287689,0.262554,0.23829,0.215153,0.193326,0.172928,0.154027,0.136645,0.120768,0.106358,0.0933528,0.0816778,0.0712478,0.0619719,0.0537571,0.0465105,0.0401416,0.0345632,0.0296931,0.0254544,0.0217759,0.0185923,0.0158441,0.0134776,0.0114446,0.00970201,0.00821149,0.00693918,0.00585524,0.0049335,0.00415108,0.00348807]
lines!(ax9, xrange, p, color = :blue, linewidth = 2)

Y = y_values[81]
ax10 = Axis(fig3[2,4],
    # Background
    backgroundcolor = :transparent,
    xgridvisible = false,
    ygridvisible = false,
    limits = ((0,6), (0,1)),
    # Title
    title = (L"\mu=%$(round(Y; digits=3))"),
    titlevisible = true,
    titlesize = 22,
    titlealign = :center,
    titlegap = 4.0,
    # x-axis
    xlabel = L"x",
    xlabelvisible = true,
    xlabelsize = 24,
    xlabelcolor = :black,
    xlabelpadding = -20.0,
    xticks = [0,2,4,6],
    xticksvisible = true,
    xticksize = 10,
    xticklabelsvisible = true,
    xticklabelsize = 24,
    xtickformat = "{:.0f}",
    xscale = identity, #log10,
    xaxisposition = :bottom,
    # y-axis
    ylabel = L"density",
    ylabelvisible = false,
    ylabelsize = 24,
    ylabelcolor = :black,
    ylabelpadding = 10,
    yticks = [0,1],
    yticksvisible = true,
    yticksize = 10,
    yticklabelsvisible = false,
    yticklabelsize = 24,
    ytickformat = "{:.1f}",
    yscale = identity,
    yaxisposition = :left,
)
hist!(ax10, Et[:,81], bins = 100, normalization = :pdf, color = :red, strokecolor = :black, strokewidth = 1)
p = [1.10344*10^-10,0.0000179559,0.000380414,0.0021665,0.0070385,0.0167159,0.0324935,0.0549571,0.0839068,0.118439,0.157122,0.198209,0.239844,0.280239,0.317807,0.351245,0.379581,0.402183,0.418738,0.429218,0.433835,0.432985,0.4272,0.417103,0.40336,0.386652,0.36764,0.346949,0.325151,0.302757,0.280209,0.25788,0.236081,0.215054,0.194989,0.17602,0.158239,0.141698,0.126415,0.112384,0.0995774,0.08795,0.0774456,0.0679996,0.0595416,0.0519989,0.0452977,0.0393654,0.0341313,0.0295278,0.0254912,0.0219616,0.0188837,0.0162067,0.013884,0.0118734,0.010137,0.00864039,0.00735323]
lines!(ax10, xrange, p, color = :blue, linewidth = 2)

# Export the figures
save("../../results/fast_slow/fig2.7.2.png", fig1)
save("../../results/fast_slow/Moments.png", fig2)
save("../../results/fast_slow/fig2.7.1.png", fig3)
