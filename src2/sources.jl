# src2/sources.jl — Build individual GammaSource objects from a Cryostat.
#
# One `GammaSource` per (physical Ti volume, isotope). Each carries the
# self-shielded inward γ flux dN/du(u) at the source's own inner surface.
# No downstream attenuation is applied here — that lives in the (later)
# effective-source builder.

"""
    GammaSource

Self-shielded gamma source at a single physical Ti volume.

Fields:
  * `name`            — human label (e.g. "OCV_barrel_Bi214")
  * `producer`        — the underlying `PCyl` or `PDisk`
  * `isotope`         — `:Bi214` or `:Tl208`
  * `E_MeV`           — 2.448 (Bi-214) or 2.615 (Tl-208)
  * `produced_per_yr` — γ/yr produced inside the Ti volume
  * `exit_inward_per_yr` — γ/yr crossing the inner surface, inward-going
  * `u_bins`          — discretization of cos θ ∈ (0, 1]
  * `dNdu`            — γ/yr per unit u at the inner surface (sums to
                         `exit_inward_per_yr` under trapezoidal integration)
"""
struct GammaSource
    name::String
    producer::Union{PCyl, PDisk}
    isotope::Symbol
    E_MeV::Float64
    produced_per_yr::Float64
    exit_inward_per_yr::Float64
    u_bins::Vector{Float64}
    dNdu::Vector{Float64}
end

"""
    _trapz(y, x) -> Float64

Trapezoidal integral of y(x). Assumes uniform-or-arbitrary spacing of x.
"""
function _trapz(y::Vector{Float64}, x::Vector{Float64})::Float64
    s = 0.0
    @inbounds for i in 1:length(x)-1
        s += 0.5 * (y[i] + y[i+1]) * (x[i+1] - x[i])
    end
    s
end

"""
    make_gamma_source(p::PObject, isotope::Symbol, u_bins) -> GammaSource

Build the self-shielded `GammaSource` for one (physical object, isotope)
pair.
"""
function make_gamma_source(p::PObject, isotope::Symbol,
                           u_bins::Vector{Float64})::GammaSource
    if isotope === :Bi214
        E       = E_BI214_MEV
        rate_pr = gamma_rate_Bi214(p)
    elseif isotope === :Tl208
        E       = E_TL208_MEV
        rate_pr = gamma_rate_Tl208(p)
    else
        error("make_gamma_source: unknown isotope $isotope")
    end
    dNdu     = self_shielded_spectrum(p, rate_pr, E, u_bins)
    exit_inw = _trapz(dNdu, u_bins)
    name     = string(p.name, "_", isotope)
    GammaSource(name, p, isotope, E, rate_pr, exit_inw, u_bins, dNdu)
end

# ---------------------------------------------------------------------------
# Cryostat → list of physical objects → list of GammaSources
# ---------------------------------------------------------------------------

# Default u-binning: 100 uniform bins in (0, 1]. Avoid u=0 (divergent path).
const DEFAULT_U_BINS = collect(range(0.005, 0.995, length=100))

# bb0nu Table I late-chain Ti specific activities (mBq/kg, upper limits).
const TI_BB0NU_U238_LATE_MBQKG = 0.08
const TI_BB0NU_TH232_LATE_MBQKG = 0.22

"""
    pobjects_from_cryostat(c::Cryostat, mat_Ti; mc_only=true) -> Vector{PObject}

Wrap each `Cryostat` element (barrels, heads, MC-active extras) as a
`PCyl` or `PDisk` carrying the bb0nu Ti specific activities. With
`mc_only=true` (default) only extras flagged `include_in_mc = true`
are returned.
"""
function pobjects_from_cryostat(c::Cryostat, mat_Ti::Material;
                                mc_only::Bool=true)::Vector{PObject}
    ps = PObject[]
    aU  = TI_BB0NU_U238_LATE_MBQKG
    aTh = TI_BB0NU_TH232_LATE_MBQKG

    # Barrels (always MC-active)
    barrel_names = ("OCV_barrel", "ICV_barrel")
    for (i, g) in enumerate(c.barrels)
        push!(ps, PCyl(g, mat_Ti, aU, aTh; count=1, name=barrel_names[i]))
    end

    # Heads (always MC-active). c.heads is ordered (OCV top, OCV bot,
    # ICV top, ICV bot) by build_cryostat.
    head_names = ("OCV_top_head", "OCV_bottom_head",
                  "ICV_top_head", "ICV_bottom_head")
    for (i, g) in enumerate(c.heads)
        push!(ps, PDisk(g, mat_Ti, aU, aTh; count=1, name=head_names[i]))
    end

    # Extras
    for e in c.extras
        if mc_only && !e.include_in_mc
            continue
        end
        push!(ps, PCyl(e.shell, mat_Ti, aU, aTh;
                       count=e.count, name=e.name))
    end

    ps
end

"""
    build_individual_sources(c::Cryostat, mat_Ti, u_bins=DEFAULT_U_BINS;
                             mc_only=true) -> Vector{GammaSource}

For each MC-active physical Ti volume in the cryostat, build a
self-shielded `GammaSource` for both Bi-214 and Tl-208. With the
default `mc_only=true`, returns 2 × N_active = 34 sources for the
current LZ model.
"""
function build_individual_sources(c::Cryostat, mat_Ti::Material,
                                  u_bins::Vector{Float64}=DEFAULT_U_BINS;
                                  mc_only::Bool=true)::Vector{GammaSource}
    ps = pobjects_from_cryostat(c, mat_Ti; mc_only=mc_only)
    out = GammaSource[]
    for p in ps
        push!(out, make_gamma_source(p, :Bi214, u_bins))
        push!(out, make_gamma_source(p, :Tl208, u_bins))
    end
    out
end
