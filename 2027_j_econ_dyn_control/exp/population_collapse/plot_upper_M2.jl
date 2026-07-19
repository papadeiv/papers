#=
plot_M2.jl

Plots (all CairoMakie, all exported as PDF):
  1) The partition of the (G0, L0) plane into Omega_1 ... Omega_6
  2) The 3D region M2 = { (G0,L0,A0) : A0 < min(U_j, A0star_j) on Omega_j }, all
     regions plotted together in one call (no region-by-region loop)
  3) A fixed-G0 slice of M2 in the (L0, A0) plane

Structured so each Omega_j still has its own small set of formula functions
(region formula + classifier check) that plot_M2_3d loops over internally to
build the combined surface -- but there is no interactive stepping/pausing;
everything renders and saves directly.
=#

using Roots

# ---------------------------------------------------------------------------
# Parameters (Lagerlof's calibration)
# ---------------------------------------------------------------------------
const τ  = 0.15
const ρ  = 0.879
const α  = 0.6
const θ  = 1.0
const L̂  = 11.42

const γ  = 0.225
const c̃  = 1.0

# ---------------------------------------------------------------------------
# Region thresholds (from eq:partition_omega)
# ---------------------------------------------------------------------------
const G0_thresh  = ρ^2 * τ / (1 - ρ)          # e0 switch
const L0_thresh1 = ρ / ((1 - ρ) * θ)          # sub-split inside Omega1/Omega2 vs Omega4/Omega5
const L0_thresh2 = L̂ / θ                      # partition_L0 switch

curve(L0) = (ρ^4 * τ) / ((1 - ρ)^3 * θ^2) / L0^2   # only meaningful for L0 > 0

# ---------------------------------------------------------------------------
# Classifier: which Omega_j does (G0, L0) belong to?
# ---------------------------------------------------------------------------
function classify(G0::Float64, L0::Float64)::Int
    if G0 < G0_thresh
        if L0 < L0_thresh1
            return 1
        elseif L0 < L0_thresh2
            return 2
        else
            return 3
        end
    else
        if L0 < L0_thresh1
            return G0 < curve(L0) ? 4 : 5
        elseif L0 < L0_thresh2
            return 5
        else
            return 6
        end
    end
end

# ---------------------------------------------------------------------------
# Per-region (e0, e1, G1) triples, exactly as derived in eq:sub_partition_1..4
# (each takes (G0, L0) but some ignore G0 or L0, matching the derivation)
# ---------------------------------------------------------------------------
region1_e0e1G1(G0, L0) = (0.0, 0.0, τ * ρ * θ * L0)

region2_e0e1G1(G0, L0) = (0.0,
                          τ * (sqrt((1 - ρ) * ρ * θ * L0) - ρ),
                          τ * ρ * θ * L0)

region3_e0e1G1(G0, L0) = (0.0,
                          τ * (sqrt((1 - ρ) * ρ * L̂) - ρ),
                          τ * ρ * L̂)

function region4_e0e1G1(G0, L0)
    σ0 = sqrt(τ * (1 - ρ) * G0)
    e0 = σ0 - ρ * τ
    G1 = θ * L0 * σ0
    return (e0, 0.0, G1)
end

function region5_e0e1G1(G0, L0)
    σ0 = sqrt(τ * (1 - ρ) * G0)
    e0 = σ0 - ρ * τ
    G1 = θ * L0 * σ0
    e1 = sqrt(τ * (1 - ρ) * θ * L0 * σ0) - ρ * τ
    return (e0, e1, G1)
end

function region6_e0e1G1(G0, L0)
    σ0 = sqrt(τ * (1 - ρ) * G0)
    e0 = σ0 - ρ * τ
    G1 = L̂ * σ0
    e1 = sqrt(τ * (1 - ρ) * L̂ * σ0) - ρ * τ
    return (e0, e1, G1)
end

const REGION_FORMULA = Dict(
    1 => region1_e0e1G1,
    2 => region2_e0e1G1,
    3 => region3_e0e1G1,
    4 => region4_e0e1G1,
    5 => region5_e0e1G1,
    6 => region6_e0e1G1,
)

# ---------------------------------------------------------------------------
# Shared machinery: eq:upper (closed form) and eq:lower (implicit, root-find)
# ---------------------------------------------------------------------------

# U_j(G0,L0): closed-form eq:upper bound on A0
function upper_bound(e0, e1, G1, G0, L0)
    Φ  = c̃ * ((e1 + ρ*τ + G1) / (e1 + ρ*τ))^α          # = c̃ * h1^{-α}
    U0 = L0 * Φ^(1/(1-α)) / ((τ + e1) * (1 + G1))       # coefficient without gamma
    return γ * U0, U0
end

# V_j(A0; G0,L0): eq:lower RHS as a function of A0 (implicit fixed point)
function lower_rhs(A0, e0, e1, G1, G0, L0, U0)
    ζ0inv = ((e0 + ρ*τ + G0) / (e0 + ρ*τ))^α            # = c̃ / z0 * (A0/L0)^{1-α}... see below
    return U0 * (1 - c̃ * ζ0inv * (L0/A0)^(1-α))
end

# Find A0* solving A0 = V_j(A0). Falls back to Inf (unconstrained) or 0.0 (empty)
# if no sign change is found in the search bracket -- flag these cases for review.
function lower_star(e0, e1, G1, G0, L0, U0; lo=1e-6, hi=nothing)
    hi === nothing && (hi = 50 * max(U0, 1.0))
    g(A0) = A0 - lower_rhs(A0, e0, e1, G1, G0, L0, U0)
    glo, ghi = g(lo), g(hi)
    if sign(glo) == sign(ghi)
        # No sign change detected in bracket: decide by sign of g at lo
        return glo > 0 ? Inf : 0.0
    end
    return find_zero(g, (lo, hi), Bisection())
end

# A0_max(G0,L0) for a single region j (returns NaN if outside that region's domain)
function region_A0max(j::Int, G0::Float64, L0::Float64)
    if classify(G0, L0) != j
        return NaN
    end
    e0, e1, G1 = REGION_FORMULA[j](G0, L0)
    γU0, U0 = upper_bound(e0, e1, G1, G0, L0)
    A0star = lower_star(e0, e1, G1, G0, L0, U0)
    return min(γU0, A0star)
end

# A0_max(G0,L0) using whichever region actually applies at that point (no j needed)
function A0max_at(G0::Float64, L0::Float64)
    j = classify(G0, L0)
    e0, e1, G1 = REGION_FORMULA[j](G0, L0)
    γU0, U0 = upper_bound(e0, e1, G1, G0, L0)
    A0star = lower_star(e0, e1, G1, G0, L0, U0)
    return min(γU0, A0star), j
end

# ---------------------------------------------------------------------------
# PLOT 1: (G0, L0) partition into Omega_1 ... Omega_6  (CairoMakie)
# ---------------------------------------------------------------------------
using CairoMakie
function plot_omega_partition(; G0max=5.0, L0max=20.0, n=600, savepath="omega_partition.pdf")
    G0s = range(0, G0max, length=n)
    L0s = range(1e-6, L0max, length=n)   # avoid L0=0 in curve()
    Z = [classify(g, l) for g in G0s, l in L0s]   # Z[i,j] indexed by (G0,L0)

    # Plot the covering of Omega 
    fig = Figure(size=(800, 650))
    ax = Axis(fig[1,1], xlabel="L0", ylabel="G0", title="Partition of the (L0,G0)-plane")
    hm = heatmap!(ax, L0s, G0s, transpose(Z), colormap=:tab10, colorrange=(1,6))

    # Boundary lines (solid black)
    vlines!(ax, [L0_thresh1], color=:black, linewidth=2)
    vlines!(ax, [L0_thresh2], color=:black, linewidth=2)
    hlines!(ax, [G0_thresh], color=:black, linewidth=2)

    # Curved boundary G0 = curve(L0), only where it lies in view and L0 < L0_thresh1
    Lc = range(1e-6, L0_thresh1, length=300)
    Gc = [curve(l) for l in Lc]
    lines!(ax, Lc, Gc, color=:black, linewidth=2)

    xlims!(ax, 0, L0max)
    ylims!(ax, 0, G0max)

    save(savepath, fig)
end

# ---------------------------------------------------------------------------
# PLOT 2: 3D surface of A0_max(G0,L0), i.e. the upper boundary of M2, colored
# per region. Static CairoMakie render (no interactive rotation) -- all
# regions are computed and drawn together in a single call, saved as PDF.
# ---------------------------------------------------------------------------
using CairoMakie
function plot_M2_3d(; G0max=5.0, L0max=20.0, n=600, regions_to_run=1:6,
                     savepath="M2_surface.png")
    G0s = collect(range(1e-6, G0max, length=n))
    L0s = collect(range(1e-6, L0max, length=n))

    fig = Figure(size=(900,750))
    ax = Axis3(fig[1,1], xlabel="G0", ylabel="L0", zlabel="A0",
               title="M2 boundary surface (all regions)")

    colors = Dict(1=>:red, 2=>:orange, 3=>:gold, 4=>:green, 5=>:blue, 6=>:purple)

    for j in regions_to_run
        A0 = [region_A0max(j, g, l) for g in G0s, l in L0s]
        if all(isnan, A0)
            @warn "Region $j produced no valid points in the given (G0max,L0max) window"
            continue
        end
        surface!(ax, G0s, L0s, A0, color=fill(colors[j], size(A0)), nan_color=:transparent)
    end

    save(savepath, fig)
end

# ---------------------------------------------------------------------------
# PLOT 3: fixed-G0 slice of M2 in the (L0, A0) plane (CairoMakie)
#
# For a given G0, sweeps L0 and computes A0_max(L0) = boundary of M2 at that
# (G0,L0). The region under the curve (0 <= A0 < A0_max) is M2's slice.
# Colored by which Omega_j applies at each L0, so region transitions along
# the slice are visible; vertical dashed lines mark where classify() switches.
# ---------------------------------------------------------------------------
using CairoMakie
function plot_M2_slice(G0::Float64; L0max=20.0, n=2000, savepath="M2_slice_G0_$(G0).pdf")
    L0s = collect(range(1e-6, L0max, length=n))
    A0s = Vector{Float64}(undef, n)
    js  = Vector{Int}(undef, n)

    for (i, L0) in enumerate(L0s)
        A0max, j = A0max_at(G0, L0)
        A0s[i] = isfinite(A0max) ? A0max : NaN
        js[i]  = j
    end

    colors = Dict(1=>:red, 2=>:orange, 3=>:gold, 4=>:green, 5=>:blue, 6=>:purple)

    fig = Figure(size=(850,600))
    ax = Axis(fig[1,1], xlabel="L0", ylabel="A0",
              title="Slice of M2 at G0 = $(G0)")

    # Shade region under the boundary curve, segment by segment so each
    # segment can carry its own region color.
    switch_idx = findall(i -> i > 1 && js[i] != js[i-1], eachindex(js))
    seg_starts = vcat(1, switch_idx)
    seg_ends   = vcat(switch_idx .- 1, n)

    for (s, e) in zip(seg_starts, seg_ends)
        j = js[s]
        band!(ax, L0s[s:e], zeros(e-s+1), A0s[s:e], color=(colors[j], 0.35))
        lines!(ax, L0s[s:e], A0s[s:e], color=colors[j], linewidth=2)
    end

    # Mark region-transition boundaries with dashed vertical lines
    for i in switch_idx
        vlines!(ax, [L0s[i]], color=:black, linestyle=:dash, linewidth=1)
    end

    xlims!(ax, 0, L0max)
    ylims!(ax, 0, nothing)

    save(savepath, fig)
end

# ---------------------------------------------------------------------------
# Driver: renders and saves all three plots directly, no pausing/looping
# ---------------------------------------------------------------------------
function main(; G0max=5.0, L0max=20.0, regions_to_run=1:6, slice_G0=1.0)
    #println("Plotting (G0,L0) partition...")
    #plot_omega_partition(G0max=G0max, L0max=L0max)

    println("Plotting combined M2 surface (all regions)...")
    plot_M2_3d(G0max=G0max, L0max=L0max, regions_to_run=regions_to_run)

    #println("Plotting M2 slice at G0 = $(slice_G0)...")
    #plot_M2_slice(slice_G0; L0max=L0max)
end

# Execute the script
main()
