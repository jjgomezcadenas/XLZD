# src2/effective_sources.jl — Effective sources at the ICV inner surfaces.
#
# An `EffectiveSource` aggregates inward γ flux from a set of upstream
# `GammaSource` objects, each filtered through its own downstream Ti
# slab list, evaluated at a common ICV-inner exit surface.
#
# For LZ we define three effective sources per isotope:
#   CB  = Cryostat Barrel        (exit at ICV cylindrical inner wall)
#   CTH = Cryostat Top Head      (exit at ICV inner top head)
#   CBH = Cryostat Bottom Head   (exit at ICV inner bottom head)
#
# 3 regions × 2 isotopes (Bi-214, Tl-208) = 6 EffectiveSource objects.

"""
    SourceContribution(source, downstream_slabs)

One `GammaSource` plus the Ti slabs the gamma must traverse to reach
an effective-source exit surface. `downstream_slabs` is empty when the
individual source already sits on the effective-source exit (e.g. ICV
barrel feeding into the CB).
"""
struct SourceContribution
    source::GammaSource
    downstream_slabs::Vector{Slab}
end

"""
    EffectiveSource

Aggregated inward γ flux at a chosen exit surface. `dNdu(u)` integrates
to `total_per_yr` over the `u_bins`.
"""
struct EffectiveSource
    name::String
    region::Symbol
    isotope::Symbol
    E_MeV::Float64
    contributions::Vector{SourceContribution}
    u_bins::Vector{Float64}
    dNdu::Vector{Float64}
    total_per_yr::Float64
end

"""
    aggregate_dNdu(contributions, u_bins) -> Vector{Float64}

Sum each contribution's `source.dNdu` weighted by the downstream
transmission factor at the same u-grid.
"""
function aggregate_dNdu(contributions::Vector{SourceContribution},
                         u_bins::Vector{Float64})::Vector{Float64}
    out = zeros(Float64, length(u_bins))
    for c in contributions
        T = transmission_factor(c.downstream_slabs, u_bins)
        @. out += c.source.dNdu * T
    end
    out
end

"""
    build_effective_source(name, region, isotope, contributions, u_bins) -> EffectiveSource
"""
function build_effective_source(name::AbstractString, region::Symbol,
                                isotope::Symbol,
                                contributions::Vector{SourceContribution},
                                u_bins::Vector{Float64})::EffectiveSource
    dNdu  = aggregate_dNdu(contributions, u_bins)
    total = _trapz(dNdu, u_bins)
    E_MeV = isotope === :Bi214 ? E_BI214_MEV : E_TL208_MEV
    EffectiveSource(String(name), region, isotope, E_MeV,
                    contributions, u_bins, dNdu, total)
end

# ---------------------------------------------------------------------------
# LZ-specific registry: three effective sources per isotope
# ---------------------------------------------------------------------------

"""
    _icv_slab_thicknesses(c::Cryostat) -> (t_body, t_top, t_bot)

Read the three ICV wall thicknesses (cm) from the Cryostat. ICV is the
second barrel and the third/fourth head per `build_cryostat`'s push order.
"""
function _icv_slab_thicknesses(c::Cryostat)
    @assert length(c.barrels) >= 2 "Cryostat: missing ICV barrel"
    @assert length(c.heads)   >= 4 "Cryostat: missing ICV heads"
    t_body = c.barrels[2].wall_thickness
    t_top  = c.heads[3].wall_thickness   # ICV top head
    t_bot  = c.heads[4].wall_thickness   # ICV bottom head
    (t_body, t_top, t_bot)
end

"""
    build_effective_sources(individual_sources, c::Cryostat, mat_Ti;
                            u_bins=DEFAULT_U_BINS) -> Vector{EffectiveSource}

Build the 6 LZ effective sources (CB, CTH, CBH × Bi214, Tl208). The
individual GammaSources must come from `build_individual_sources` on
the same cryostat.
"""
function build_effective_sources(individual_sources::Vector{GammaSource},
                                  c::Cryostat, mat_Ti::Material;
                                  u_bins::Vector{Float64}=DEFAULT_U_BINS
                                  )::Vector{EffectiveSource}

    by_name = Dict{String, GammaSource}()
    for s in individual_sources
        by_name[s.name] = s
    end

    # Helper: lookup with informative error
    function get_src(name::String)
        haskey(by_name, name) ||
            error("build_effective_sources: missing GammaSource '$name'. " *
                  "Available: $(sort(collect(keys(by_name))))")
        by_name[name]
    end

    t_ICV_body, t_ICV_top, t_ICV_bot = _icv_slab_thicknesses(c)

    out = EffectiveSource[]

    for iso in (:Bi214, :Tl208)
        E   = iso === :Bi214 ? E_BI214_MEV : E_TL208_MEV
        μ   = mat_Ti.μ_lin(E)
        slab_body = Slab(μ, t_ICV_body, "ICV_body")
        slab_top  = Slab(μ, t_ICV_top,  "ICV_top_head")
        slab_bot  = Slab(μ, t_ICV_bot,  "ICV_bottom_head")

        # ── CB: Cryostat Barrel (exit at ICV barrel inner wall) ─────────
        cb_contribs = [
            # already at ICV inner surface — no downstream
            SourceContribution(get_src("ICV_barrel_$iso"),                  Slab[]),
            SourceContribution(get_src("ICV_main_flange_pair_$iso"),        Slab[]),
            # OCV-attached → must traverse ICV body wall
            SourceContribution(get_src("OCV_barrel_$iso"),                  [slab_body]),
            SourceContribution(get_src("OCV_upper_body_flange_pair_$iso"),  [slab_body]),
            SourceContribution(get_src("OCV_lower_body_flange_pair_$iso"),  [slab_body]),
            SourceContribution(get_src("OCV_HV_port_$iso"),                 [slab_body]),
            SourceContribution(get_src("OCV_cathode_side_port_$iso"),       [slab_body]),
        ]
        push!(out, build_effective_source("CB", :barrel, iso, cb_contribs, u_bins))

        # ── CTH: Cryostat Top Head (exit at ICV inner top head) ─────────
        cth_contribs = [
            SourceContribution(get_src("ICV_top_head_$iso"),                Slab[]),
            SourceContribution(get_src("OCV_top_head_$iso"),                [slab_top]),
            SourceContribution(get_src("OCV_top_reinforcing_ring_$iso"),    [slab_top]),
            SourceContribution(get_src("OCV_top_YBe_recess_$iso"),          [slab_top]),
            SourceContribution(get_src("OCV_top_conduit_ports_$iso"),       [slab_top]),
            SourceContribution(get_src("OCV_top_tie_bar_ports_$iso"),       [slab_top]),
        ]
        push!(out, build_effective_source("CTH", :endcap_top, iso, cth_contribs, u_bins))

        # ── CBH: Cryostat Bottom Head (exit at ICV inner bottom head) ──
        cbh_contribs = [
            SourceContribution(get_src("ICV_bottom_head_$iso"),             Slab[]),
            SourceContribution(get_src("OCV_bottom_head_$iso"),             [slab_bot]),
            SourceContribution(get_src("OCV_bottom_support_flange_$iso"),   [slab_bot]),
            SourceContribution(get_src("OCV_bottom_umbilical_port_$iso"),   [slab_bot]),
        ]
        push!(out, build_effective_source("CBH", :endcap_bottom, iso, cbh_contribs, u_bins))
    end

    out
end
