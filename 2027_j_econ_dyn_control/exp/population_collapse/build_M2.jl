# ---------------------------------------------------------------------------
# IMPORTANT: inequality (30) gives an upper bound A0. However inequality (31) 
# is implicit in A0; solving it yields either an interval (A0_low, A0_high), 
# or an empty set. The interval is then intersected with the solution of (30) 
# to get the final surface.
# ---------------------------------------------------------------------------

using Roots
using CairoMakie
using MAT

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
# Partitions of Ω from equations (42)
# ---------------------------------------------------------------------------
const G1 = ρ^2 * τ / (1 - ρ)
const L1 = ρ / ((1 - ρ) * θ)
const L2 = L̂ / θ
curve(L0) = (ρ^4 * τ) / ((1 - ρ)^3 * θ^2) / L0^2   # only meaningful for L0 > 0

# ---------------------------------------------------------------------------
# Classifier for Ωj
# ---------------------------------------------------------------------------
function classify(G0::Float64, L0::Float64)::Int
    if G0 < G1
        if L0 < L1
            return 1
        elseif L0 < L2
            return 2
        else
            return 3
        end
    else
        if L0 < L1
            return G0 < curve(L0) ? 4 : 5
        elseif L0 < L2
            return 5
        else
            return 6
        end
    end
end

# ---------------------------------------------------------------------------
# Triple of values (e0, e1, G1) in each Ωj as derived in equations (37)~(41) 
# ---------------------------------------------------------------------------
region1(G0, L0) = (0.0, 
                   0.0, 
                   τ * ρ * θ * L0)

region2(G0, L0) = (0.0,
                   τ * (sqrt((1 - ρ) * ρ * θ * L0) - ρ),
                   τ * ρ * θ * L0)

region3(G0, L0) = (0.0,
                   τ * (sqrt((1 - ρ) * ρ * L̂) - ρ),
                   τ * ρ * L̂)

function region4(G0, L0)
    σ0 = sqrt(τ * (1 - ρ) * G0)
    e0 = σ0 - ρ * τ
    G1 = θ * L0 * σ0
    return (e0, 0.0, G1)
end

function region5(G0, L0)
    σ0 = sqrt(τ * (1 - ρ) * G0)
    e0 = σ0 - ρ * τ
    G1 = θ * L0 * σ0
    e1 = sqrt(τ * (1 - ρ) * θ * L0 * σ0) - ρ * τ
    return (e0, e1, G1)
end

function region6(G0, L0)
    σ0 = sqrt(τ * (1 - ρ) * G0)
    e0 = σ0 - ρ * τ
    G1 = L̂ * σ0
    e1 = sqrt(τ * (1 - ρ) * L̂ * σ0) - ρ * τ
    return (e0, e1, G1)
end

const REGION_FORMULA = Dict(
    1 => region1,
    2 => region2,
    3 => region3,
    4 => region4,
    5 => region5,
    6 => region6,
)

# ---------------------------------------------------------------------------
# Construction of the surfaces (30) and (31) given triples (e0, e1, G1)
# ---------------------------------------------------------------------------

function upper_bound(e0, e1, G1, G0, L0)
    h1_inv = ((e1 + ρ*τ + G1) / (e1 + ρ*τ))^α
    A = L0 * (c̃*h1_inv)^(1/(1-α)) / ((τ + e1) * (1 + G1))
    return γ * A, A
end

function lower_band(e0, e1, G1, G0, L0, A)
    h0_inv = ((e0 + ρ*τ + G0) / (e0 + ρ*τ))^α
    B = A * c̃ * h0_inv * L0^(1-α)

    # Handle the degenerate case in which the c̃, h0_inv or L0 are negative
    if B <= 0
        return (0.0, Inf)
    end

    # Handle the degenerate case in which (31) is not satisfied (empty region)
    A0star = (B * (1-α))^(1/(2-α))
    φ(A0)  = A0 + B * A0^(-(1-α))
    φmin   = φ(A0star)
    if φmin >= A 
        return (NaN, NaN)
    end

    # Find the root for the lower branch
    g(A0) = φ(A0) - A 
    A0_low = find_zero(g, (1e-10, A0star), Bisection())

    # Find the root for the upper branch
    high = max(10*A0star, 10*A, 1.0)
    while g(high) < 0
        high *= 10
    end
    A0_high = find_zero(g, (A0star, high), Bisection())

    return (A0_low, A0_high)
end

# Impose subsistence condition to derive A0 
function z0_window(e0, G0, L0)
    h0_inv = ((e0 + ρ*τ + G0) / (e0 + ρ*τ))^α
    A0_zmin = L0 * (c̃ * h0_inv)^(1/(1-α))
    A0_zmax = L0 * (c̃ * h0_inv / (1-γ))^(1/(1-α))
    return A0_zmin, A0_zmax
end

# Interset the result from (30) and (31) to compute M_2^{(j)}
function region_A0_band(j::Int, G0::Float64, L0::Float64)
    if classify(G0, L0) != j
        return (NaN, NaN)
    end
    e0, e1, G1 = REGION_FORMULA[j](G0, L0)
    U1, U0 = upper_bound(e0, e1, G1, G0, L0)
    A0_low, A0_high = lower_band(e0, e1, G1, G0, L0, U0)
    isnan(A0_low) && return (NaN, NaN)

    A0_high_final = min(A0_high, U1)
    A0_low >= A0_high_final && return (NaN, NaN)

    # Enforce the subsistence condition 
    A0_zmin, A0_zmax = z0_window(e0, G0, L0)
    A0_low_final = max(A0_low, A0_zmin)
    A0_high_final = min(A0_high_final, A0_zmax)
    A0_low_final >= A0_high_final && return (NaN, NaN)

    return (A0_low_final, A0_high_final)
end

# Band using whichever region actually applies at (G0,L0), plus the region index.
function A0band_at(G0::Float64, L0::Float64)
    j = classify(G0, L0)
    e0, e1, G1 = REGION_FORMULA[j](G0, L0)
    U1, U0 = upper_bound(e0, e1, G1, G0, L0)
    A0_low, A0_high = lower_band(e0, e1, G1, G0, L0, U0)
    if isnan(A0_low)
        return (NaN, NaN, j)
    end
    A0_high_final = min(A0_high, U1)
    if A0_low >= A0_high_final
        return (NaN, NaN, j)
    end

    A0_zmin, A0_zmax = z0_window(e0, G0, L0)
    A0_low_final = max(A0_low, A0_zmin)
    A0_high_final = min(A0_high_final, A0_zmax)
    if A0_low_final >= A0_high_final
        return (NaN, NaN, j)
    end

    return (A0_low_final, A0_high_final, j)
end

# ---------------------------------------------------------------------------
# Plot the partition of ℝ2≥0 into the covering of Ωj
# ---------------------------------------------------------------------------
function plot_omega_partition(; G0max=5.0, L0max=20.0, n=600, savepath="./results/omega_partition.pdf")
    G0s = range(0, G0max, length=n)
    L0s = range(1e-6, L0max, length=n)   # avoid L0=0 in curve()
    Z = [classify(g, l) for g in G0s, l in L0s]   # Z[i,j] indexed by (G0,L0)

    # Plot the covering of Omega 
    fig = Figure(size=(800, 650))
    ax = Axis(fig[1,1], xlabel="L0", ylabel="G0", title="Partition of the (L0,G0)-plane")
    hm = heatmap!(ax, L0s, G0s, transpose(Z), colormap=:tab10, colorrange=(1,6))

    # Boundary lines (solid black)
    vlines!(ax, [L1], color=:black, linewidth=2)
    vlines!(ax, [L2], color=:black, linewidth=2)
    hlines!(ax, [G1], color=:black, linewidth=2)

    # Curved boundary G0 = curve(L0), only where it lies in view and L0 < L0_thresh1
    Lc = range(1e-6, L1, length=300)
    Gc = [curve(l) for l in Lc]
    lines!(ax, Lc, Gc, color=:black, linewidth=2)

    xlims!(ax, 0, L0max)
    ylims!(ax, 0, G0max)

    save(savepath, fig)
end

# ---------------------------------------------------------------------------
# Plot the boundary of M2
# ---------------------------------------------------------------------------
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

        surface!(ax, G0s, L0s, A0_lo, color=fill((colors[j], surface_alpha), size(A0_lo)),
                 nan_color=:transparent)
        surface!(ax, G0s, L0s, A0_hi, color=fill((colors[j], surface_alpha), size(A0_hi)),
                 nan_color=:transparent)

        wireframe!(ax, G0s[wire_idx_G0], L0s[wire_idx_L0], A0_lo[wire_idx_G0, wire_idx_L0],
                   color=(:black, wire_alpha), linewidth=1)
        wireframe!(ax, G0s[wire_idx_G0], L0s[wire_idx_L0], A0_hi[wire_idx_G0, wire_idx_L0],
                   color=(:black, wire_alpha), linewidth=1)
    end

    # SLice at Lagerlof's G0=0.048
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
# Plot the the projection onto (A0,L0)-plane for fixed G0
# ---------------------------------------------------------------------------
function plot_M2_slice(G0::Float64; L0max=20.0, n=2000, savepath="./results/M2_slice_G0_$(G0).pdf")
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

    fig = Figure(size=(850,600))
    ax = Axis(fig[1,1], xlabel="L0", ylabel="A0",
              title="Slice of M2 at G0 = $(G0)")

    switch_idx = findall(i -> i > 1 && js[i] != js[i-1], eachindex(js))
    seg_starts = vcat(1, switch_idx)
    seg_ends   = vcat(switch_idx .- 1, n)

    for (s, e) in zip(seg_starts, seg_ends)
        j = js[s]
        idx = s:e
        valid = .!isnan.(A0lo[idx])
        any(valid) || continue
        band!(ax, L0s[idx], A0lo[idx], A0hi[idx], color=(colors[j], 0.35))
        lines!(ax, L0s[idx], A0lo[idx], color=colors[j], linewidth=4)
        lines!(ax, L0s[idx], A0hi[idx], color=colors[j], linewidth=4)
    end

    for i in switch_idx
        vlines!(ax, [L0s[i]], color=:black, linestyle=:dash, linewidth=2)
    end

    xlims!(ax, 0, L0max)
    ylims!(ax, 0, nothing)

    save(savepath, fig)
end

# ---------------------------------------------------------------------------
# Export the M2 boundaries (A0_lo and A0_hi) at a fixed G0 to .mat, for import
# in MATLAB. Sweeps L0 exactly as plot_M2_slice does, and writes:
#   G0      -- scalar, the fixed G0 value
#   L0      -- vector, the L0 grid
#   A0_lo   -- vector, lower branch of the M2 boundary (NaN where empty)
#   A0_hi   -- vector, upper branch of the M2 boundary (NaN where empty)
#   region  -- vector of Int, which Omega_j applies at each L0 (for reference)
# ---------------------------------------------------------------------------
function export_M2_slice_mat(G0::Float64; L0max=20.0, n=2000,
                              savepath="./results/M2_slice_G0_$(G0).mat")
    L0s  = collect(range(1e-6, L0max, length=n))
    A0lo = Vector{Float64}(undef, n)
    A0hi = Vector{Float64}(undef, n)
    js   = Vector{Int}(undef, n)

    for (i, L0) in enumerate(L0s)
        lo, hi, j = A0band_at(G0, L0)
        A0lo[i] = lo
        A0hi[i] = hi
        js[i]   = j
    end

    matwrite(savepath, Dict(
        "G0"      => G0,
        "L0"      => L0s,
        "A0_low"  => A0lo,
        "A0_high" => A0hi,
        "region"  => js,
    ))
end

# ---------------------------------------------------------------------------
# Driver: main script
# ---------------------------------------------------------------------------
function main(; G0max=0.5, L0max=20.0, regions_to_run=1:6, slice_G0=0.048)
    println("Plotting ∪_j Ωj")
    plot_omega_partition(G0max=10*G0max, L0max=L0max)

    println("Plotting ∂M2")
    plot_M2_3d(G0max=G0max, L0max=L0max, regions_to_run=regions_to_run,
               highlight_G0=slice_G0)

    println("Plotting and exporting Π_G0(M2) (G0 = $(slice_G0))")
    plot_M2_slice(slice_G0; L0max=L0max)
    export_M2_slice_mat(slice_G0; L0max=L0max)
end

# Execute the script
main()
