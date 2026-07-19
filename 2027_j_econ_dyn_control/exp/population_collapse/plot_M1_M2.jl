#=
plot_M2.jl

Plots (all CairoMakie, all exported as PDF):
  1) The partition of the (G0, L0) plane into Omega_1 ... Omega_6
  2) The pair of surfaces A0_lo(G0,L0) and A0_hi(G0,L0) that envelope M2, i.e.
     M2 = { (G0,L0,A0) : A0_lo(G0,L0) < A0 < A0_hi(G0,L0) }, all regions
     plotted together in one call (no region-by-region loop)
  3) A fixed-G0 slice of M2 in the (L0, A0) plane, shaded between the two
     bounding curves

IMPORTANT: eq:upper alone gives a simple ceiling on A0, but eq:lower is
implicit in A0 (A0 appears on both sides via z0) and does NOT reduce to a
simple ceiling -- solving it yields a bounded band (A0_lo, A0_hi), or empty,
never a floor-at-zero ceiling. See lower_band() below for the derivation.
This band is then intersected with eq:upper's ceiling to get the final
per-region envelope.

Structured so each Omega_j still has its own small set of formula functions
(region formula + classifier check) that plot_M2_3d loops over internally to
build the combined surfaces -- but there is no interactive stepping/pausing;
everything renders and saves directly.
=#

using Roots
using ProgressMeter

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
# M1 boundary: M1 = { (G0,A0,L0) in M : A0 < L0*(c_tilde*h(G0)^-alpha)^(1/(1-alpha)) },
# where h(G0) = (e0+rho*tau)/(e0+rho*tau+G0) by analogy with h1=h(G1)
# (eq:accessory_h1), and e0(G0) is the (region-independent) master formula
# from eq:accessory_e0.
# ---------------------------------------------------------------------------
e0_of(G0) = max(0.0, -ρ*τ + sqrt(τ * (1-ρ) * G0))
h_of(G0)  = (e0_of(G0) + ρ*τ) / (e0_of(G0) + ρ*τ + G0)

# At fixed G0, M1's boundary in the (L0,A0) plane is a straight line through
# the origin with this slope.
M1_slope(G0) = (c̃ * h_of(G0)^(-α))^(1/(1-α))

# ---------------------------------------------------------------------------
# Shared machinery: eq:upper (closed form) and eq:lower (implicit BAND)
#
# IMPORTANT: eq:lower is *not* a simple ceiling on A0. Writing
#   A0 < U0 * (1 - c_tilde * zeta0inv * (L0/A0)^(1-alpha))
# and rearranging gives  phi(A0) := A0 + B*A0^-(1-alpha) < U0 , where
#   B := U0 * c_tilde * zeta0inv * L0^(1-alpha).
# Since 1-alpha > 0, phi(A0) -> +Inf as A0 -> 0+ AND as A0 -> Inf, with a
# single interior minimum at A0star = [(1-alpha)*B]^(1/(2-alpha)). So the
# solution set of eq:lower alone is either EMPTY (if phi(A0star) >= U0) or a
# genuine bounded interval (A0_lo, A0_hi) straddling A0star. This band must
# then be intersected with eq:upper's ceiling U1 = gamma*U0.
# ---------------------------------------------------------------------------

# U_j(G0,L0): closed-form eq:upper bound on A0 (U1 = gamma*U0, U0 = bare coefficient)
function upper_bound(e0, e1, G1, G0, L0)
    Φ  = c̃ * ((e1 + ρ*τ + G1) / (e1 + ρ*τ))^α          # = c̃ * h1^{-α}
    U0 = L0 * Φ^(1/(1-α)) / ((τ + e1) * (1 + G1))       # coefficient without gamma
    return γ * U0, U0
end

# Solve eq:lower alone for the (A0_lo, A0_hi) band. Returns (NaN, NaN) if empty.
function lower_band(e0, e1, G1, G0, L0, U0)
    ζ0inv = ((e0 + ρ*τ + G0) / (e0 + ρ*τ))^α            # = c̃ / z0 * (A0/L0)^{1-α} coefficient
    B = U0 * c̃ * ζ0inv * L0^(1-α)

    if B <= 0
        return (0.0, Inf)   # degenerate: condition never binds
    end

    A0star = (B * (1-α))^(1/(2-α))
    φ(A0)  = A0 + B * A0^(-(1-α))
    φmin   = φ(A0star)

    if φmin >= U0
        return (NaN, NaN)   # eq:lower is never satisfiable here -> region empty
    end

    g(A0) = φ(A0) - U0

    # Lower root: bracket (eps, A0star). g(eps) > 0 (phi -> +Inf), g(A0star) < 0.
    A0_lo = find_zero(g, (1e-10, A0star), Bisection())

    # Upper root: bracket (A0star, hi), expanding hi until g flips positive.
    hi = max(10*A0star, 10*U0, 1.0)
    while g(hi) < 0
        hi *= 10
    end
    A0_hi = find_zero(g, (A0star, hi), Bisection())

    return (A0_lo, A0_hi)
end

# z0-validity window: eq:restricted_optimum's second branch only applies for
#   c_tilde <= z0 < c_tilde/(1-gamma).
# Since z0(A0) = zeta0(e0,G0) * (A0/L0)^(1-alpha) = (A0/L0)^(1-alpha) / zeta0inv,
# solving both endpoints for A0 gives:
#   A0_zmin = L0 * (c_tilde * zeta0inv)^(1/(1-alpha))
#   A0_zmax = L0 * (c_tilde * zeta0inv / (1-gamma))^(1/(1-alpha))
function z0_window(e0, G0, L0)
    ζ0inv = ((e0 + ρ*τ + G0) / (e0 + ρ*τ))^α
    A0_zmin = L0 * (c̃ * ζ0inv)^(1/(1-α))
    A0_zmax = L0 * (c̃ * ζ0inv / (1-γ))^(1/(1-α))
    return A0_zmin, A0_zmax
end

# Combine eq:upper and eq:lower into the final (A0_lo, A0_hi) band for region j.
# Returns (NaN, NaN) if outside region j's domain, or if the combined band is empty.
function region_A0_band(j::Int, G0::Float64, L0::Float64)
    if classify(G0, L0) != j
        return (NaN, NaN)
    end
    e0, e1, G1 = REGION_FORMULA[j](G0, L0)
    U1, U0 = upper_bound(e0, e1, G1, G0, L0)
    A0_lo, A0_hi = lower_band(e0, e1, G1, G0, L0, U0)
    isnan(A0_lo) && return (NaN, NaN)

    A0_hi_final = min(A0_hi, U1)
    A0_lo >= A0_hi_final && return (NaN, NaN)   # eq:upper cuts the band away entirely

    # Enforce the z0-validity window (branch 2 of eq:restricted_optimum only
    # applies while c_tilde <= z0 < c_tilde/(1-gamma)).
    A0_zmin, A0_zmax = z0_window(e0, G0, L0)
    A0_lo_final = max(A0_lo, A0_zmin)
    A0_hi_final = min(A0_hi_final, A0_zmax)
    A0_lo_final >= A0_hi_final && return (NaN, NaN)

    return (A0_lo_final, A0_hi_final)
end

# Band using whichever region actually applies at (G0,L0), plus the region index.
function A0band_at(G0::Float64, L0::Float64)
    j = classify(G0, L0)
    e0, e1, G1 = REGION_FORMULA[j](G0, L0)
    U1, U0 = upper_bound(e0, e1, G1, G0, L0)
    A0_lo, A0_hi = lower_band(e0, e1, G1, G0, L0, U0)
    if isnan(A0_lo)
        return (NaN, NaN, j)
    end
    A0_hi_final = min(A0_hi, U1)
    if A0_lo >= A0_hi_final
        return (NaN, NaN, j)
    end

    A0_zmin, A0_zmax = z0_window(e0, G0, L0)
    A0_lo_final = max(A0_lo, A0_zmin)
    A0_hi_final = min(A0_hi_final, A0_zmax)
    if A0_lo_final >= A0_hi_final
        return (NaN, NaN, j)
    end

    return (A0_lo_final, A0_hi_final, j)
end

# ---------------------------------------------------------------------------
# PLOT 1: (G0, L0) partition into Omega_1 ... Omega_6  (CairoMakie)
# ---------------------------------------------------------------------------
using CairoMakie
function plot_omega_partition(; G0max=5.0, L0max=20.0, n=600, savepath="./results/omega_partition.pdf")
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
# PLOT 2: 3D surfaces of A0_lo(G0,L0) and A0_hi(G0,L0), the pair that envelopes
# M2. Both surfaces drawn per region (semi-transparent upper, solid lower),
# all regions combined in one call, saved as PDF.
# ---------------------------------------------------------------------------
using CairoMakie
function plot_M2_3d(; G0max=5.0, L0max=20.0, n=600, regions_to_run=1:6,
                     savepath="./results/M2_surface.png", surface_alpha=0.05,
                     wire_stride=10, wire_alpha=0.00,
                     highlight_G0=nothing, highlight_n=1500, highlight_eps=nothing)
    G0s = collect(range(1e-6, G0max, length=n))
    L0s = collect(range(1e-6, L0max, length=n))

    fig = Figure(size=(900,750))
    ax = Axis3(fig[1,1], xlabel="G0", ylabel="L0", zlabel="A0",
               title="M2 enveloping surfaces (lower + upper, all regions)")

    colors = Dict(1=>:red, 2=>:orange, 3=>:gold, 4=>:green, 5=>:blue, 6=>:purple)

    # Coarser index subset for the wireframe overlay, so it reads as a grid
    # rather than a dense mesh of lines on top of the (already fine) surface.
    wire_idx_G0 = 1:wire_stride:length(G0s)
    wire_idx_L0 = 1:wire_stride:length(L0s)

    for j in regions_to_run
        A0_lo = Matrix{Float64}(undef, length(G0s), length(L0s))
        A0_hi = Matrix{Float64}(undef, length(G0s), length(L0s))
        for (ig, g) in enumerate(G0s), (il, l) in enumerate(L0s)
            lo, hi = region_A0_band(j, g, l)
            A0_lo[ig, il] = lo
            A0_hi[ig, il] = hi
        end

        if all(isnan, A0_lo)
            @warn "Region $j produced no valid (nonempty) points in the given window"
            continue
        end

        # Transparent fills for both the lower and upper surfaces.
        surface!(ax, G0s, L0s, A0_lo, color=fill((colors[j], surface_alpha), size(A0_lo)),
                 nan_color=:transparent)
        surface!(ax, G0s, L0s, A0_hi, color=fill((colors[j], surface_alpha), size(A0_hi)),
                 nan_color=:transparent)

        # Wireframe overlay on top of each, on a coarser grid, for shape legibility.
        wireframe!(ax, G0s[wire_idx_G0], L0s[wire_idx_L0], A0_lo[wire_idx_G0, wire_idx_L0],
                   color=(#=colors[j]=#:black, wire_alpha), linewidth=1)
        wireframe!(ax, G0s[wire_idx_G0], L0s[wire_idx_L0], A0_hi[wire_idx_G0, wire_idx_L0],
                   color=(#=colors[j]=#:black, wire_alpha), linewidth=1)
    end

    # ---------------------------------------------------------------------
    # Highlight: the slice of M2 at G0 = highlight_G0, drawn fully opaque.
    # Built as a thin vertical "ribbon" surface -- two rows at (essentially)
    # the same G0 value, with z = A0_lo(L0) on the first row and z = A0_hi(L0)
    # on the second, so the strip in between traces exactly the plane's
    # intersection with M2. Segmented by region (color + gaps at empty spots),
    # same logic as plot_M2_slice.
    # ---------------------------------------------------------------------
    if highlight_G0 !== nothing
        highlight_eps === nothing && (highlight_eps = 1e-4 * max(G0max, 1.0))

        L0s_h = collect(range(1e-6, L0max, length=highlight_n))
        A0lo_h = Vector{Float64}(undef, highlight_n)
        A0hi_h = Vector{Float64}(undef, highlight_n)
        js_h   = Vector{Int}(undef, highlight_n)
        for (i, L0) in enumerate(L0s_h)
            lo, hi, j = A0band_at(highlight_G0, L0)
            A0lo_h[i] = lo
            A0hi_h[i] = hi
            js_h[i]   = j
        end

        switch_idx = findall(i -> i > 1 && js_h[i] != js_h[i-1], eachindex(js_h))
        seg_starts = vcat(1, switch_idx)
        seg_ends   = vcat(switch_idx .- 1, highlight_n)

        xrow = [highlight_G0 - highlight_eps/2, highlight_G0 + highlight_eps/2]

        for (s, e) in zip(seg_starts, seg_ends)
            j = js_h[s]
            idx = s:e
            valid = .!isnan.(A0lo_h[idx])
            any(valid) || continue

            Z = Matrix{Float64}(undef, 2, length(idx))
            Z[1, :] = A0lo_h[idx]
            Z[2, :] = A0hi_h[idx]

            surface!(ax, xrow, L0s_h[idx], Z, color=fill(colors[j], size(Z)),
                     nan_color=:transparent)
        end
    end

    save(savepath, fig)
end

# ---------------------------------------------------------------------------
# PLOT 3: fixed-G0 slice of M2 in the (L0, A0) plane (CairoMakie)
#
# For a given G0, sweeps L0 and computes (A0_lo(L0), A0_hi(L0)), the band that
# envelopes M2 at that (G0,L0). Colored by which Omega_j applies at each L0,
# so region transitions along the slice are visible; vertical dashed lines
# mark where classify() switches. Points where the band is empty (NaN) are
# simply not shaded.
# ---------------------------------------------------------------------------
using CairoMakie
function plot_M2_slice(idx::Integer, G0::Float64; L0max=20.0, n=2000, savepath="./results/slices/$(idx).png")
    L0s   = collect(range(1e-6, L0max, length=n))
    A0lo  = Vector{Float64}(undef, n)
    A0hi  = Vector{Float64}(undef, n)
    js    = Vector{Int}(undef, n)

    for (i, L0) in enumerate(L0s)
        lo, hi, j = A0band_at(G0, L0)
        A0lo[i] = lo
        A0hi[i] = hi
        js[i]   = j
    end

    colors = Dict(1=>:red, 2=>:orange, 3=>:gold, 4=>:green, 5=>:blue, 6=>:purple)

    fig = Figure(size=(850, 600))
    ax = Axis(fig[1,1], xlabel="L0", ylabel="A0", title="Slices of M1 and M2 at G0 = $(G0)")

    # Shade the band between A0_lo and A0_hi, segment by segment so each
    # segment can carry its own region color; skip segments where the band
    # is empty (NaN).
    switch_idx = findall(i -> i > 1 && js[i] != js[i-1], eachindex(js))
    seg_starts = vcat(1, switch_idx)
    seg_ends   = vcat(switch_idx .- 1, n)

    for (s, e) in zip(seg_starts, seg_ends)
        j = js[s]
        idx = s:e
        valid = .!isnan.(A0lo[idx])
        any(valid) || continue
        band!(ax, L0s[idx], A0lo[idx], A0hi[idx], color=(colors[j], 0.35))
        lines!(ax, L0s[idx], A0lo[idx], color=colors[j], linewidth=3)
        lines!(ax, L0s[idx], A0hi[idx], color=colors[j], linewidth=3)
    end

    # Overlay M1's boundary at this G0: a straight line through the origin,
    # A0 = L0 * M1_slope(G0). Shade the area between A0=0 and this line (i.e.
    # M1 itself) in blue.
    slope = M1_slope(G0)
    A0_M1 = slope .* L0s
    band!(ax, L0s, zeros(n), A0_M1, color=(:blue, 0.35))
    lines!(ax, L0s, A0_M1, color=:blue, linewidth=3)

    # Mark region-transition boundaries with dashed vertical lines
    for i in switch_idx
        vlines!(ax, [L0s[i]], color=:black, linestyle=:dash, linewidth=2)
    end

    xlims!(ax, 0, L0max)
    ylims!(ax, 0, 225)

    save(savepath, fig)
end

# ---------------------------------------------------------------------------
# Driver: renders and saves all three plots directly, no pausing/looping
# ---------------------------------------------------------------------------
function main(; G0max=0.5, L0max=20.0, regions_to_run=1:6, slice_G0=0.5)
        #=
        println("Plotting (G0,L0) partition...")
        plot_omega_partition(G0max=10*G0max, L0max=L0max)

        println("Plotting combined M2 surface (all regions)...")
        plot_M2_3d(G0max=G0max, L0max=L0max, regions_to_run=regions_to_run, highlight_G0=slice_G0)
        =#

        println("Plotting M2 slice at G0 = $(slice_G0)...")
        @showprogress for (idx, G0) in enumerate(LinRange(0.00, slice_G0, 100)) 
                plot_M2_slice(idx, G0; L0max=L0max)
        end
end

# Execute the script
main()
