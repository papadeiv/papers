# ---------------------------------------------------------------------------
# IMPORTANT: M_2^{(j)}(G0,L0) is a UNION of two branch-consistent intervals in
# A0, not a single interval satisfying both (30) and (31) at once. (30)/(31)
# come from two mutually exclusive branches of n(G0,A0,L0) (eq:lagerlof_optimum),
# selected by comparing z0(A0) against c̃ and c̃/(1-γ) -- only one branch is
# ever active for a given A0, so the two conditions must be OR'd (restricted
# each to its own branch's z0-window), not AND'd:
#
#   branch 1 (z0 >= c̃/(1-γ)): only (30) applies -> [A0_zmax, γ*U0)
#   branch 2 (c̃ <= z0 < c̃/(1-γ)): only (31) applies -> [A0_zmin,A0_zmax) ∩ (A0_low,A0_high)
#
#   M_2^{(j)}(G0,L0) = branch1 ∪ branch2
#
# The two pieces may or may not be contiguous; when A0_high < A0_zmax there is
# a genuine gap between them. Every consumer below carries both pieces
# through as (b1_lo,b1_hi,b2_lo,b2_hi), each NaN when empty.
#
# Also included: the boundary of M1 (collapse after 1 iteration), the
# dynamical system (eq:lagerlof_system) simulated with a plain iteration
# loop, and a plot classifying initial conditions by collapse time against
# ∂M1, ∂M2.
# ---------------------------------------------------------------------------

using Roots
using CairoMakie
using MAT
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
    tries = 0
    while g(high) < 0 && tries < 60
        high *= 10
        tries += 1
    end
    if g(high) < 0 || !isfinite(high)
        return (NaN, NaN)
    end
    A0_high = find_zero(g, (A0star, high), Bisection())
    if !isfinite(A0_low) || !isfinite(A0_high)
        return (NaN, NaN)
    end

    return (A0_low, A0_high)
end

# The window of A0 for which z0(A0) lies in [c̃, c̃/(1-γ)) (branch 2's window).
# Branch 1's window is simply A0 >= A0_zmax (z0 above the window).
function z0_window(e0, G0, L0)
    h0_inv = ((e0 + ρ*τ + G0) / (e0 + ρ*τ))^α
    A0_zmin = L0 * (c̃ * h0_inv)^(1/(1-α))
    A0_zmax = L0 * (c̃ * h0_inv / (1-γ))^(1/(1-α))
    return A0_zmin, A0_zmax
end

# Build the two branch-consistent intervals of M_2^{(j)} at (G0,L0):
#   b1 = [A0_zmax, γ*U0)                             (branch 1: only (30) applies)
#   b2 = [A0_zmin,A0_zmax) ∩ (A0_low,A0_high)         (branch 2: only (31) applies)
# M_2^{(j)}(G0,L0) = b1 ∪ b2. Each returned as (lo,hi), NaN if empty.
function region_A0_bands(j::Int, G0::Float64, L0::Float64)
    if classify(G0, L0) != j
        return (NaN, NaN, NaN, NaN)
    end
    e0, e1, G1 = REGION_FORMULA[j](G0, L0)
    U1, U0 = upper_bound(e0, e1, G1, G0, L0)
    A0_zmin, A0_zmax = z0_window(e0, G0, L0)

    # branch 1: (30) restricted to z0 >= c̃/(1-γ)
    b1_lo, b1_hi = A0_zmax, U1
    if !(b1_lo < b1_hi)
        b1_lo, b1_hi = NaN, NaN
    end

    # branch 2: (31) restricted to c̃ <= z0 < c̃/(1-γ)
    A0_low, A0_high = lower_band(e0, e1, G1, G0, L0, U0)
    if isnan(A0_low)
        b2_lo, b2_hi = NaN, NaN
    else
        b2_lo = max(A0_low, A0_zmin)
        b2_hi = min(A0_high, A0_zmax)
        if !(b2_lo < b2_hi)
            b2_lo, b2_hi = NaN, NaN
        end
    end

    return (b1_lo, b1_hi, b2_lo, b2_hi)
end

# Bands using whichever region actually applies at (G0,L0), plus the region index.
function A0_bands_at(G0::Float64, L0::Float64)
    j = classify(G0, L0)
    b1_lo, b1_hi, b2_lo, b2_hi = region_A0_bands(j, G0, L0)
    return (b1_lo, b1_hi, b2_lo, b2_hi, j)
end

# ---------------------------------------------------------------------------
# Observables e(G), h(G), z(G,A,L) from equations (lagerlof_education),
# (lagerlof_capital), (lagerlof_income)
# ---------------------------------------------------------------------------
e(G) = max(0.0, -ρ*τ + sqrt(τ * (1-ρ) * G))
h(G) = (e(G) + ρ*τ) / (e(G) + ρ*τ + G)
z(G, A, L) = h(G)^α * (A/L)^(1-α)

# ---------------------------------------------------------------------------
# Boundary of M1 (population collapse after 1 iteration)
# ---------------------------------------------------------------------------
M1_slope(G0) = (c̃ * h(G0)^(-α))^(1/(1-α))

# ---------------------------------------------------------------------------
# Dynamical system (Lagerlöf map, eq:lagerlof_system) and simulator
# ---------------------------------------------------------------------------

# n(G,A,L) from eq:lagerlof_optimum; L<=0 is absorbing (0 population stays 0),
# handled before z(G,A,L) is evaluated to avoid a division by L=0.
function n(G, A, L)
    L <= 0.0 && return 0.0
    G1 = (e(G) + ρ*τ) * min(θ*L, L̂)   # G_{t+1}, needed for e_{t+1}
    e1 = e(G1)
    z0 = z(G, A, L)
    if z0 >= c̃/(1-γ)
        return γ / (τ + e1)
    elseif z0 >= c̃
        return (1 - c̃/z0) / (τ + e1)
    else
        return 0.0
    end
end

# f(G_t,A_t,L_t) from eq:lagerlof_system
function lagerlof_map(G, A, L)
    G1 = (e(G) + ρ*τ) * min(θ*L, L̂)
    A1 = (1 + G1) * A
    L1 = n(G, A, L) * L
    return (G1, A1, L1)
end

# Simulate T iterations of f starting from (G0,A0,L0) with a plain loop;
# returns the full trajectory as a Vector of (G,A,L) tuples, T+1 entries
# (t=0 included as entry 1).
function simulate(G0::Float64, A0::Float64, L0::Float64; T::Int=60)
    traj = Vector{NTuple{3,Float64}}(undef, T+1)
    traj[1] = (G0, A0, L0)
    for t in 1:T
        traj[t+1] = lagerlof_map(traj[t]...)
    end
    return traj
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
# Plot the boundary of M2 (both branches, all regions)
# ---------------------------------------------------------------------------
function plot_M2_3d(; G0max=5.0, L0max=20.0, n=600, regions_to_run=1:6,
                     savepath="./results/M2_surface.png", surface_alpha=0.25,
                     wire_stride=10, wire_alpha=0.00,
                     highlight_G0=nothing, highlight_n=1500, highlight_eps=nothing)
    G0s = collect(range(1e-6, G0max, length=n))
    L0s = collect(range(1e-6, L0max, length=n))

    fig = Figure(size=(900,750))
    ax = Axis3(fig[1,1], xlabel="G0", ylabel="L0", zlabel="A0",
               title="M2 enveloping surfaces (branch 1 + branch 2, all regions)")

    colors = Dict(1=>:red, 2=>:orange, 3=>:gold, 4=>:green, 5=>:blue, 6=>:purple)

    wire_idx_G0 = 1:wire_stride:length(G0s)
    wire_idx_L0 = 1:wire_stride:length(L0s)

    for j in regions_to_run
        B1_lo = Matrix{Float64}(undef, length(G0s), length(L0s))
        B1_hi = Matrix{Float64}(undef, length(G0s), length(L0s))
        B2_lo = Matrix{Float64}(undef, length(G0s), length(L0s))
        B2_hi = Matrix{Float64}(undef, length(G0s), length(L0s))
        for (ig, g) in enumerate(G0s), (il, l) in enumerate(L0s)
            b1lo, b1hi, b2lo, b2hi = region_A0_bands(j, g, l)
            B1_lo[ig, il] = b1lo
            B1_hi[ig, il] = b1hi
            B2_lo[ig, il] = b2lo
            B2_hi[ig, il] = b2hi
        end

        if all(isnan, B1_lo) && all(isnan, B2_lo)
            @warn "Region $j produced no valid (nonempty) points in the given window"
            continue
        end

        for M in (B1_lo, B1_hi, B2_lo, B2_hi)
            all(isnan, M) && continue
            surface!(ax, G0s, L0s, M, color=fill((colors[j], surface_alpha), size(M)),
                     nan_color=:transparent)
            wireframe!(ax, G0s[wire_idx_G0], L0s[wire_idx_L0], M[wire_idx_G0, wire_idx_L0],
                       color=(:black, wire_alpha), linewidth=1)
        end
    end

    # SLice at Lagerlof's G0=0.048
    if highlight_G0 !== nothing
        highlight_eps === nothing && (highlight_eps = 1e-4 * max(G0max, 1.0))

        L0s_h = collect(range(1e-6, L0max, length=highlight_n))
        B1lo_h = Vector{Float64}(undef, highlight_n)
        B1hi_h = Vector{Float64}(undef, highlight_n)
        B2lo_h = Vector{Float64}(undef, highlight_n)
        B2hi_h = Vector{Float64}(undef, highlight_n)
        js_h   = Vector{Int}(undef, highlight_n)
        for (i, L0) in enumerate(L0s_h)
            b1lo, b1hi, b2lo, b2hi, j = A0_bands_at(highlight_G0, L0)
            B1lo_h[i] = b1lo; B1hi_h[i] = b1hi
            B2lo_h[i] = b2lo; B2hi_h[i] = b2hi
            js_h[i]   = j
        end

        switch_idx = findall(i -> i > 1 && js_h[i] != js_h[i-1], eachindex(js_h))
        seg_starts = vcat(1, switch_idx)
        seg_ends   = vcat(switch_idx .- 1, highlight_n)

        xrow = [highlight_G0 - highlight_eps/2, highlight_G0 + highlight_eps/2]

        for (s, e) in zip(seg_starts, seg_ends)
            j = js_h[s]
            idx = s:e
            for (Blo, Bhi) in ((B1lo_h, B1hi_h), (B2lo_h, B2hi_h))
                valid = .!isnan.(Blo[idx])
                any(valid) || continue
                Z = Matrix{Float64}(undef, 2, length(idx))
                Z[1, :] = Blo[idx]
                Z[2, :] = Bhi[idx]
                surface!(ax, xrow, L0s_h[idx], Z, color=fill(colors[j], size(Z)),
                         nan_color=:transparent)
            end
        end
    end

    save(savepath, fig)
end

# ---------------------------------------------------------------------------
# Plot the projection of M1, M2 onto (A0,L0)-plane for fixed G0
# ---------------------------------------------------------------------------
function plot_M2_slice(G0::Float64; L0max=20.0, n=2000, savepath="./results/M2_slice_G0_$(G0).pdf")
    L0s   = collect(range(1e-6, L0max, length=n))
    B1lo  = Vector{Float64}(undef, n)
    B1hi  = Vector{Float64}(undef, n)
    B2lo  = Vector{Float64}(undef, n)
    B2hi  = Vector{Float64}(undef, n)
    js    = Vector{Int}(undef, n)

    for (i, L0) in enumerate(L0s)
        b1lo, b1hi, b2lo, b2hi, j = A0_bands_at(G0, L0)
        B1lo[i] = b1lo; B1hi[i] = b1hi
        B2lo[i] = b2lo; B2hi[i] = b2hi
        js[i]   = j
    end

    colors = Dict(1=>:red, 2=>:orange, 3=>:gold, 4=>:green, 5=>:blue, 6=>:purple)

    fig = Figure(size=(850,600))
    ax = Axis(fig[1,1], xlabel="L0", ylabel="A0",
              title="Slice of M1, M2 at G0 = $(G0)")

    switch_idx = findall(i -> i > 1 && js[i] != js[i-1], eachindex(js))
    seg_starts = vcat(1, switch_idx)
    seg_ends   = vcat(switch_idx .- 1, n)

    for (s, e) in zip(seg_starts, seg_ends)
        j = js[s]
        idx = s:e
        for (Blo, Bhi) in ((B1lo, B1hi), (B2lo, B2hi))
            valid = .!isnan.(Blo[idx])
            any(valid) || continue
            band!(ax, L0s[idx], Blo[idx], Bhi[idx], color=(colors[j], 0.35))
            lines!(ax, L0s[idx], Blo[idx], color=colors[j], linewidth=4)
            lines!(ax, L0s[idx], Bhi[idx], color=colors[j], linewidth=4)
        end
    end

    # M1's boundary: straight line through the origin, A0 = L0*M1_slope(G0).
    # Shade the area between A0=0 and the line (i.e. M1 itself) in blue.
    slope = M1_slope(G0)
    A0_M1 = slope .* L0s
    band!(ax, L0s, zeros(n), A0_M1, color=(:blue, 0.35))
    lines!(ax, L0s, A0_M1, color=:blue, linewidth=4)

    for i in switch_idx
        vlines!(ax, [L0s[i]], color=:black, linestyle=:dash, linewidth=2)
    end

    xlims!(ax, 0, L0max)
    ylims!(ax, 0, nothing)

    save(savepath, fig)
end

# ---------------------------------------------------------------------------
# Export the M2 boundaries (both branches) at a fixed G0 to .mat, for import
# in MATLAB. Sweeps L0 exactly as plot_M2_slice does, and writes:
#   G0                  -- scalar, the fixed G0 value
#   L0                  -- vector, the L0 grid
#   branch1_low/high    -- branch 1 piece of the M2 boundary (NaN where empty)
#   branch2_low/high    -- branch 2 piece of the M2 boundary (NaN where empty)
#   region              -- vector of Int, which Omega_j applies at each L0
# ---------------------------------------------------------------------------
function export_M2_slice_mat(G0::Float64; L0max=20.0, n=2000,
                              savepath="./results/slices/$(G0).mat")
    L0s  = collect(range(1e-6, L0max, length=n))
    B1lo = Vector{Float64}(undef, n)
    B1hi = Vector{Float64}(undef, n)
    B2lo = Vector{Float64}(undef, n)
    B2hi = Vector{Float64}(undef, n)
    js   = Vector{Int}(undef, n)

    for (i, L0) in enumerate(L0s)
        b1lo, b1hi, b2lo, b2hi, j = A0_bands_at(G0, L0)
        B1lo[i] = b1lo; B1hi[i] = b1hi
        B2lo[i] = b2lo; B2hi[i] = b2hi
        js[i]   = j
    end

    matwrite(savepath, Dict(
        "G0"           => G0,
        "L0"           => L0s,
        "branch1_low"  => B1lo,
        "branch1_high" => B1hi,
        "branch2_low"  => B2lo,
        "branch2_high" => B2hi,
        "region"       => js,
    ))
end

# ---------------------------------------------------------------------------
# Plot collapse-time classification of a grid of initial conditions at fixed
# G0, with ∂M1 (dark blue) and ∂M2 (purple, shaded by Ωj, both branches)
# overlaid.
#   blue      : simulate(...) gives L=0 already after 1 iteration
#   deeppink  : L=0 for the first time after 2 iterations
# Branch 1 of ∂M2 is drawn solid, branch 2 dashed, same purple shade per Ωj.
# ---------------------------------------------------------------------------
function plot_collapse(idx::Integer, G0::Float64; L0max=20.0, A0max=225.0, n=150, T=60,
                savepath="./results/scan/$(idx).png")
    L0s = range(1e-6, L0max, length=n)
    A0s = range(1e-6, A0max, length=n)

    # Class[iL,iA] = 1.0 (collapse after 1 iter), 2.0 (after 2 iter), NaN (neither)
    Class = fill(NaN, n, n)

    for (iL, L0) in enumerate(L0s), (iA, A0) in enumerate(A0s)
        traj = simulate(G0, A0, L0; T=T)
        L1 = traj[2][3]   # L at t=1
        L2 = traj[3][3]   # L at t=2
        if L1 == 0.0
            Class[iL, iA] = 1.0
        elseif L2 == 0.0
            Class[iL, iA] = 2.0
        end
    end

    fig = Figure(size=(850,600))
    ax = Axis(fig[1,1], xlabel="L0", ylabel="A0", title = "G0 = $(G0)")

    # Squares centered at each grid point, tiled continuously (no gaps between
    # points, unlike scatter markers); transparent where neither class applies.
    heatmap!(ax, L0s, A0s, Class, colormap=[:royalblue1, :pink1],
             colorrange=(1,2), nan_color=:transparent)

    # ∂M1: thick dark-blue line (darker shade of the blue points above)
    L0line = collect(range(1e-6, L0max, length=1000))
    A0_M1  = M1_slope(G0) .* L0line
    lines!(ax, L0line, A0_M1, color=:navy, linewidth=4)

    # ∂M2: both branches, shaded purple by Ωj (branch 1 solid, branch 2 dashed)
    B1lo = Vector{Float64}(undef, length(L0line))
    B1hi = Vector{Float64}(undef, length(L0line))
    B2lo = Vector{Float64}(undef, length(L0line))
    B2hi = Vector{Float64}(undef, length(L0line))
    js   = Vector{Int}(undef, length(L0line))
    for (i, L0) in enumerate(L0line)
        b1lo, b1hi, b2lo, b2hi, j = A0_bands_at(G0, L0)
        B1lo[i] = b1lo; B1hi[i] = b1hi
        B2lo[i] = b2lo; B2hi[i] = b2hi
        js[i]   = j
    end

    switch_idx = findall(i -> i > 1 && js[i] != js[i-1], eachindex(js))
    seg_starts = vcat(1, switch_idx)
    seg_ends   = vcat(switch_idx .- 1, length(js))

    for (s, e) in zip(seg_starts, seg_ends)
        j = js[s]
        idx = s:e
        if any(.!isnan.(B1lo[idx]))
            #lines!(ax, L0line[idx], B1lo[idx], color=:purple, linewidth=4, linestyle=:solid)
            lines!(ax, L0line[idx], B1hi[idx], color=:purple, linewidth=4)
        end
        if any(.!isnan.(B2lo[idx]))
            lines!(ax, L0line[idx], B2lo[idx], color=:purple, linewidth=4)
            #lines!(ax, L0line[idx], B2hi[idx], color=:navy, linewidth=4, linestyle=:dash)
        end
    end

    xlims!(ax, 0, L0max)
    ylims!(ax, 0, A0max)

    save(savepath, fig)
end

# ---------------------------------------------------------------------------
# Driver: main script
# ---------------------------------------------------------------------------
function main(; G0max=0.5, L0max=20.0, regions_to_run=1:6, slice_G0=0.5)
    #=
    println("Plotting ∪_j Ωj")
    plot_omega_partition(G0max=10*G0max, L0max=L0max)

    println("Plotting ∂M2")
    plot_M2_3d(G0max=G0max, L0max=L0max, regions_to_run=regions_to_run,
               highlight_G0=slice_G0)
    =#

    println("Plotting and exporting Π_G0(M2) (G0 = $(slice_G0))")
    @showprogress for (idx, G0) in enumerate(LinRange(0.00, slice_G0, 100)) 
            export_M2_slice_mat(G0; L0max=L0max)
            plot_collapse(idx, G0; L0max=L0max)
    end
end

# Execute the script
main()
