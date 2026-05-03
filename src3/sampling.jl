# src3/sampling.jl — Sample (entry point, direction) for one γ from an
# EffectiveSource at its ICV inner exit surface.
#
# u = cos θ is drawn from the source's numerical dN/du(u) via inverse CDF;
# the azimuth around the local inward normal is uniform in [0, 2π).
# Position is uniform on the exit surface:
#   :barrel        → ICV cylindrical inner wall, z uniform on [z_RFR_bottom, z_gate]
#   :endcap_top    → ICV top inner head (2:1 oblate ellipsoid), uniform-by-area
#   :endcap_bottom → ICV bottom inner head (3:1 oblate ellipsoid), uniform-by-area
#
# Head geometry: the 2:1 / 3:1 aspect ratios match the cryostat
# (lz_cryo_geometry.csv); equator-z values are derived from the
# LXeDetector's z_ICV_top and z_LXe_bottom fields.

# Aspect ratios of the ICV inner heads (matches cryostat geometry CSV)
const ICV_TOP_ASPECT = 2.0
const ICV_BOT_ASPECT = 3.0

# ---------------------------------------------------------------------------
# Inverse-CDF sampling on a numerical dN/du(u) array
# ---------------------------------------------------------------------------

"""
    build_cdf(u_bins, dNdu) -> Vector{Float64}

Cumulative distribution of `dN/du`, normalised to 1 at the last bin.
`length(cdf) == length(u_bins)`. Trapezoid rule between adjacent bins.
"""
function build_cdf(u_bins::Vector{Float64}, dNdu::Vector{Float64})::Vector{Float64}
    n = length(u_bins)
    @assert length(dNdu) == n
    cdf = zeros(Float64, n)
    @inbounds for i in 2:n
        cdf[i] = cdf[i-1] + 0.5 * (dNdu[i-1] + dNdu[i]) * (u_bins[i] - u_bins[i-1])
    end
    total = cdf[end]
    total > 0 || error("build_cdf: dNdu integrates to zero")
    cdf ./= total
    cdf
end

"""
    sample_u(rng, u_bins, cdf) -> Float64

Draw u from the empirical distribution defined by `cdf` (built once via
`build_cdf`). Linearly interpolates between the bracketing bins.
"""
function sample_u(rng::AbstractRNG, u_bins::Vector{Float64},
                  cdf::Vector{Float64})::Float64
    r = rand(rng)
    idx = searchsortedfirst(cdf, r)
    if idx <= 1
        return u_bins[1]
    elseif idx > length(cdf)
        return u_bins[end]
    end
    # Linear interpolation between (cdf[idx-1], u_bins[idx-1]) and
    # (cdf[idx], u_bins[idx])
    c1 = cdf[idx-1]; c2 = cdf[idx]
    u1 = u_bins[idx-1]; u2 = u_bins[idx]
    if c2 == c1
        return u1
    end
    u1 + (r - c1) * (u2 - u1) / (c2 - c1)
end

# ---------------------------------------------------------------------------
# Per-region position + direction sampling
# ---------------------------------------------------------------------------

"""
    sample_barrel_entry(rng, det, cdf, u_bins) -> (x, y, z, dx, dy, dz, u)

Sample one γ entering the LXe through the ICV inner cylindrical wall.
Position uniform on the cylinder R = R_ICV_inner over z ∈ [z_RFR_bottom,
z_gate]; direction at polar angle θ from the local inward normal
(u = cos θ from `cdf`), azimuth uniform. The 7th return value is the
sampled `u` itself, exposed so callers can histogram the entry-angle
distribution as a Cut-1 diagnostic.
"""
function sample_barrel_entry(rng::AbstractRNG, det::LXeDetector,
                              cdf::Vector{Float64}, u_bins::Vector{Float64})
    R   = det.R_ICV_inner
    φ   = 2π * rand(rng)
    z   = det.z_RFR_bottom + (det.z_gate - det.z_RFR_bottom) * rand(rng)
    x   = R * cos(φ)
    y   = R * sin(φ)
    nx, ny, nz = -cos(φ), -sin(φ), 0.0
    u   = sample_u(rng, u_bins, cdf)
    ψ   = 2π * rand(rng)
    dx, dy, dz = rotate_direction(nx, ny, nz, u, ψ)
    (x, y, z, dx, dy, dz, u)
end

"""
    icv_top_inner_disk(det) -> GDisk

Build the GDisk representing the inner surface of the ICV top head.
Equator at z_ICV_top − R_ICV_inner / ICV_TOP_ASPECT.
"""
icv_top_inner_disk(det::LXeDetector) = GDisk(
    det.R_ICV_inner, 0.0,
    det.z_ICV_top - det.R_ICV_inner / ICV_TOP_ASPECT,
    ICV_TOP_ASPECT, :up
)

"""
    icv_bot_inner_disk(det) -> GDisk

Inner surface of the ICV bottom head. Equator at z_LXe_bottom +
R_ICV_inner / ICV_BOT_ASPECT.
"""
icv_bot_inner_disk(det::LXeDetector) = GDisk(
    det.R_ICV_inner, 0.0,
    det.z_LXe_bottom + det.R_ICV_inner / ICV_BOT_ASPECT,
    ICV_BOT_ASPECT, :down
)

"""
    sample_endcap_entry(rng, det, region, cdf, u_bins) -> (x, y, z, dx, dy, dz, u)

Sample one γ entering through the ICV inner top (`:endcap_top`) or
bottom (`:endcap_bottom`) head. Position uniform-by-area on the
ellipsoidal surface (delegated to GDisk.sample_inner_surface);
direction at polar angle θ from the local inward normal. The 7th return
value is the sampled `u` (Cut-1 diagnostic).

Note: the top head sits *above* the LXe surface, so γ from CTH start
in the gas region and must propagate to the liquid surface before any
LXe interaction. The MC handles that propagation; this function just
emits the (point, direction) pair.
"""
function sample_endcap_entry(rng::AbstractRNG, det::LXeDetector,
                              region::Symbol,
                              cdf::Vector{Float64},
                              u_bins::Vector{Float64})
    disk = if region === :endcap_top
        icv_top_inner_disk(det)
    elseif region === :endcap_bottom
        icv_bot_inner_disk(det)
    else
        error("sample_endcap_entry: unknown region $region")
    end
    x, y, z   = sample_inner_surface(rng, disk)
    nx, ny, nz = inward_normal(disk, x, y, z)
    u   = sample_u(rng, u_bins, cdf)
    ψ   = 2π * rand(rng)
    dx, dy, dz = rotate_direction(nx, ny, nz, u, ψ)
    (x, y, z, dx, dy, dz, u)
end

"""
    sample_entry(rng, det, eff::EffectiveSource) -> (x, y, z, dx, dy, dz, u)
    sample_entry(rng, det, eff::EffectiveSource, cdf::Vector{Float64})
        -> (x, y, z, dx, dy, dz, u)

Sample one (entry point, direction) for the given EffectiveSource by
dispatching on its region. Returns the sampled cos θ_inward as the 7th
element so callers can histogram it (Cut-1 diagnostic).

The 3-arg form rebuilds the CDF on every call — convenient for tests
and one-off scripts but **not** for production MC (allocates ~2N times
per thread, which crashes the GC at high N).

The 4-arg form takes a precomputed CDF (from `build_cdf(eff.u_bins,
eff.dNdu)`); production callers should build the CDF once per source
outside the event loop and pass it in.
"""
function sample_entry(rng::AbstractRNG, det::LXeDetector, eff::EffectiveSource)
    cdf = build_cdf(eff.u_bins, eff.dNdu)
    sample_entry(rng, det, eff, cdf)
end

function sample_entry(rng::AbstractRNG, det::LXeDetector,
                       eff::EffectiveSource, cdf::Vector{Float64})
    if eff.region === :barrel
        return sample_barrel_entry(rng, det, cdf, eff.u_bins)
    elseif eff.region === :endcap_top || eff.region === :endcap_bottom
        return sample_endcap_entry(rng, det, eff.region, cdf, eff.u_bins)
    elseif eff.region === :fc_barrel
        # FC barrel: read geometry from the (single) contribution's PCyl.
        p = eff.contributions[1].source.producer::PCyl
        return sample_fc_barrel_entry(rng, p.geom, cdf, eff.u_bins)
    elseif eff.region === :fc_annular
        # FC annular grid holder: read geometry from PAnnularDisk.
        p = eff.contributions[1].source.producer::PAnnularDisk
        return sample_fc_annular_entry(rng, p, cdf, eff.u_bins)
    else
        error("sample_entry: unknown region $(eff.region)")
    end
end

# ---------------------------------------------------------------------------
# Field-cage sampling primitives
# ---------------------------------------------------------------------------

"""
    sample_fc_barrel_entry(rng, g::GCyl, cdf, u_bins) -> (x, y, z, dx, dy, dz, u)

Sample one γ entering the LXe from a field-cage barrel source. Position
uniform on the inner cylindrical surface at radius `g.R_inner`, with z
in `[g.z_min, g.z_max]`. Direction at polar angle θ from the inward
normal (= −r̂), azimuth uniform.

Used for FCRN, FCRS, FCSE (R = 74.3 cm) and FCPT (R = 72.8 cm); all four
share `z ∈ [-13.75, 145.6] cm` (full FC barrel envelope).
"""
function sample_fc_barrel_entry(rng::AbstractRNG, g::GCyl,
                                 cdf::Vector{Float64}, u_bins::Vector{Float64})
    R   = g.R_inner
    φ   = 2π * rand(rng)
    z   = g.z_min + height(g) * rand(rng)
    x   = R * cos(φ)
    y   = R * sin(φ)
    nx, ny, nz = -cos(φ), -sin(φ), 0.0   # inward = −r̂
    u   = sample_u(rng, u_bins, cdf)
    ψ   = 2π * rand(rng)
    dx, dy, dz = rotate_direction(nx, ny, nz, u, ψ)
    (x, y, z, dx, dy, dz, u)
end

"""
    sample_fc_annular_entry(rng, p::PAnnularDisk, cdf, u_bins)
        -> (x, y, z, dx, dy, dz, u)

Sample one γ exiting the LXe-facing face of a grid-holder slab. Position
uniform on the annular face at z = `p.z_face`, with `r² ~ Uniform(R_in²,
R_out²)` and `φ ~ Uniform(0, 2π)`. Inward normal is `+ẑ` (FCBG, slab
opens downward) or `−ẑ` (FCTG, slab opens upward); set by
`p.normal_sign`. Direction at polar angle θ from the inward normal,
azimuth uniform.
"""
function sample_fc_annular_entry(rng::AbstractRNG, p::PAnnularDisk,
                                  cdf::Vector{Float64}, u_bins::Vector{Float64})
    φ   = 2π * rand(rng)
    r²  = p.R_in^2 + (p.R_out^2 - p.R_in^2) * rand(rng)
    r   = sqrt(r²)
    x   = r * cos(φ)
    y   = r * sin(φ)
    z   = p.z_face
    nx, ny, nz = 0.0, 0.0, Float64(p.normal_sign)   # ±ẑ inward
    u   = sample_u(rng, u_bins, cdf)
    ψ   = 2π * rand(rng)
    dx, dy, dz = rotate_direction(nx, ny, nz, u, ψ)
    (x, y, z, dx, dy, dz, u)
end
