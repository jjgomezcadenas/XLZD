# src3/histograms.jl — Diagnostic histograms for the per-photon MC.
#
# Three histogram families:
#
#   * `StackHistogramSet`     — three diagnostic plots derived from the
#                                full PhotonStack: interaction-type
#                                frequency, per-photon LXe path length,
#                                region × interaction matrix.
#   * `ClusterHistogramSet`   — per-cluster Ec spectrum (per-cluster
#                                fill, every fast-pass event).
#   * `CutHistograms`         — the cut-flow set: one accumulator per
#                                stage of the analysis funnel. See its
#                                docstring for the 8 sub-fields.
#
# Bin indexing helper `_bin_idx` is shared. Constructors take kwargs so
# detector-specific binning can be tuned at the call site.

# ---------------------------------------------------------------------------
# Bin index helper
# ---------------------------------------------------------------------------

"Bin index in 1..n for `x ∈ [lo, hi)`, else 0."
@inline function _bin_idx(x::Real, lo::Real, hi::Real, n::Int)::Int
    x < lo  && return 0
    x >= hi && return 0
    return clamp(floor(Int, (x - lo) / (hi - lo) * n) + 1, 1, n)
end

# ===========================================================================
# StackHistogramSet — per-event quantities from the PhotonStack.
# ===========================================================================
"""
    StackHistogramSet

Per-event diagnostic histograms accumulated from the rows of
`PhotonStack`. Stripped to three plots that the cut-flow refactor keeps
on the `diagnostics.png` panel:

  - `interaction_type_freq` 4 bins (PHOTO/COMPTON/PAIR/BELOW_THRESH) —
                            one entry per stack row, so normalised the
                            bar heights match the cross-section ratios
                            at the source energy.
  - `path_length_LXe`       1D, [0, 500] cm — per-photon total distance
                            through any LXe region (TPC / Skin / Inert).
                            Sourced from `stack.path_length_LXe`,
                            accumulated by the tracker.
  - `region_interaction`    3×4 matrix — per-region × per-interaction
                            counts, one entry per stack row.

`region_interaction_counts[r, t]`:
  rows: 1=:TPC, 2=:Skin, 3=:Inert
  cols: 1=PHOTO, 2=COMPTON, 3=PAIR, 4=BELOW_THRESH
"""
struct StackHistogramSet
    interaction_type_freq::Vector{Int}          # 4 bins
    path_length_n_bins::Int
    path_length_max_cm::Float64
    path_length_LXe_counts::Vector{Int}
    region_interaction_counts::Matrix{Int}      # 3×4
end

const _STACK_REGIONS      = (:TPC, :Skin, :Inert)
const _STACK_INTERACTIONS = (INT_PHOTO, INT_COMPTON, INT_PAIR, INT_BELOW_THRESH)

function StackHistogramSet(;
        path_length_n_bins::Int = 100,
        path_length_max_cm::Real = 500.0)
    StackHistogramSet(
        zeros(Int, length(_STACK_INTERACTIONS)),
        path_length_n_bins, Float64(path_length_max_cm),
        zeros(Int, path_length_n_bins),
        zeros(Int, length(_STACK_REGIONS), length(_STACK_INTERACTIONS)),
    )
end

@inline function _bin_int!(v::Vector{Int}, k::Int)
    i = clamp(k, 0, length(v) - 1) + 1
    v[i] += 1
    nothing
end

"Return the column index (1..4) for an interaction symbol, or 0 if unknown."
@inline function _interaction_col(s::Symbol)
    s === INT_PHOTO        && return 1
    s === INT_COMPTON      && return 2
    s === INT_PAIR         && return 3
    s === INT_BELOW_THRESH && return 4
    return 0
end

"Return the row index (1..3) for a region symbol, or 0 if unknown."
@inline function _stack_region_row(s::Symbol)
    s === :TPC   && return 1
    s === :Skin  && return 2
    s === :Inert && return 3
    return 0
end

"""
    update_stack_histograms!(sh::StackHistogramSet, stack::PhotonStack,
                              params::MCParams)

Fold one event's stack into `sh`. Increments:
  - one entry per stack row in `interaction_type_freq` and
    `region_interaction_counts`;
  - one entry in `path_length_LXe_counts` from the photon's accumulated
    `stack.path_length_LXe`. Filled even when the stack is empty (the
    photon may have traversed gas / FC / outside without depositing).
"""
function update_stack_histograms!(sh::StackHistogramSet,
                                   stack::PhotonStack,
                                   params::MCParams)
    @inbounds for r in stack.rows
        ci = _interaction_col(r.interaction)
        if ci > 0
            sh.interaction_type_freq[ci] += 1
            ri = _stack_region_row(r.region)
            ri > 0 && (sh.region_interaction_counts[ri, ci] += 1)
        end
    end
    ip = _bin_idx(stack.path_length_LXe, 0.0,
                   sh.path_length_max_cm, sh.path_length_n_bins)
    ip > 0 && (sh.path_length_LXe_counts[ip] += 1)
    nothing
end

function merge_stack_histograms!(into::StackHistogramSet,
                                  src::StackHistogramSet)
    @assert into.path_length_n_bins == src.path_length_n_bins
    @. into.interaction_type_freq      += src.interaction_type_freq
    @. into.path_length_LXe_counts     += src.path_length_LXe_counts
    @. into.region_interaction_counts  += src.region_interaction_counts
    into
end

# ===========================================================================
# ClusterHistogramSet — per-event quantities from the cluster vector.
# ===========================================================================
"""
    ClusterHistogramSet

Per-cluster diagnostic histogram. Stripped to a single field after the
cut-flow refactor: per-cluster `Ec` (one entry per cluster across all
fast-pass events). The cut-flow side of the analysis lives in
`CutHistograms`.
"""
struct ClusterHistogramSet
    E_n_bins::Int
    E_max_MeV::Float64
    Ec_counts::Vector{Int}
end

function ClusterHistogramSet(; E_n_bins::Int = 270, E_max_MeV::Real = 2.7,
                              kwargs...)         # ignore legacy kwargs
    ClusterHistogramSet(E_n_bins, Float64(E_max_MeV), zeros(Int, E_n_bins))
end

"""
    update_cluster_histograms!(ch::ClusterHistogramSet,
                                clusters::Vector{Cluster}, params::MCParams)

Fold one event's clusters into `ch`. One entry per cluster is added to
`Ec_counts`. Events with no clusters contribute nothing.
"""
function update_cluster_histograms!(ch::ClusterHistogramSet,
                                     clusters::Vector{Cluster},
                                     params::MCParams)
    @inbounds for c in clusters
        ie = _bin_idx(c.ec, 0.0, ch.E_max_MeV, ch.E_n_bins)
        ie > 0 && (ch.Ec_counts[ie] += 1)
    end
    nothing
end

function merge_cluster_histograms!(into::ClusterHistogramSet,
                                    src::ClusterHistogramSet)
    @assert into.E_n_bins == src.E_n_bins
    @. into.Ec_counts += src.Ec_counts
    into
end

# ===========================================================================
# CutHistograms — one accumulator per cut stage of the analysis funnel.
# ===========================================================================
"""
    CutHistograms

Cut-flow diagnostic histograms, one set of bins per cut stage:

  Cut 1 (into detector):
    - `h_u_sampled`           1D, [0, 1], `u = cos θ_inward` of every
                              sampled photon (before any propagation).

  Cut 2 (pass skin + FV):
    - `first_interaction_r_z` 2D r × z heatmap of the *first*
                              interaction position for every event with
                              an interaction (filled in `fast_veto`).

  Cut 3 (SS / MS classification — filled for cut-2-pass events):
    - `dz_inclusive`          1D, [0, dz_max] cm, every consecutive
                              visible-cluster |Δz| (one entry per
                              consecutive pair, not per event).
    - `n_visible`             1D integer, # visible clusters per event.
    - `E_total`               1D, [0, E_max] MeV, Σ ec across the
                              visible clusters of the event.

  Cut 4 (SS in ROI — filled for SS_in_ROI ∪ SS_outside_ROI):
    - `ss_ec_pre_roi`         1D, [0, E_max] MeV, true cluster energy.
    - `ss_es_pre_roi`         1D, [0, E_max] MeV, smeared cluster
                              energy (what the ROI cut acts on).
    - `ss_r_z`                2D r × z, SS cluster centroid.

All bin edges and conditioning live with the histogram itself so the
plotter can read them off without consulting MCParams.
"""
struct CutHistograms
    # Cut 1
    u_n_bins::Int
    h_u_sampled::Vector{Int}

    # Cut 2 / Cut 4 share the (r, z) binning
    r_n_bins::Int
    r_max_cm::Float64
    z_n_bins::Int
    z_min_cm::Float64
    z_max_cm::Float64
    first_interaction_r_z::Matrix{Int}

    # Cut 3
    dz_n_bins::Int
    dz_max_cm::Float64
    dz_inclusive::Vector{Int}
    n_visible_n_bins::Int
    n_visible::Vector{Int}
    E_n_bins::Int
    E_max_MeV::Float64
    E_total::Vector{Int}

    # Cut 4
    ss_ec_pre_roi::Vector{Int}
    ss_es_pre_roi::Vector{Int}
    ss_r_z::Matrix{Int}
end

function CutHistograms(;
        u_n_bins::Int = 100,
        r_n_bins::Int = 100, r_max_cm::Real = 82.1,
        z_n_bins::Int = 100, z_min_cm::Real = -69.0, z_max_cm::Real = 145.6,
        dz_n_bins::Int = 200, dz_max_cm::Real = 5.0,        # zoom near 3mm threshold
        n_visible_n_bins::Int = 21,                          # 0..20
        E_n_bins::Int = 270, E_max_MeV::Real = 2.7)
    CutHistograms(
        u_n_bins,
        zeros(Int, u_n_bins),
        r_n_bins,  Float64(r_max_cm),
        z_n_bins,  Float64(z_min_cm), Float64(z_max_cm),
        zeros(Int, r_n_bins, z_n_bins),
        dz_n_bins, Float64(dz_max_cm),
        zeros(Int, dz_n_bins),
        n_visible_n_bins,
        zeros(Int, n_visible_n_bins),
        E_n_bins,  Float64(E_max_MeV),
        zeros(Int, E_n_bins),
        zeros(Int, E_n_bins),
        zeros(Int, E_n_bins),
        zeros(Int, r_n_bins, z_n_bins),
    )
end

# --- Cut 1 ---------------------------------------------------------------

"Bin one sampled `u = cos θ_inward` into `h_u_sampled`."
@inline function fill_cut1_u!(ch::CutHistograms, u::Float64)
    iu = _bin_idx(u, 0.0, 1.0, ch.u_n_bins)
    iu > 0 && (ch.h_u_sampled[iu] += 1)
    nothing
end

# --- Cut 2 ---------------------------------------------------------------

"Bin one first-interaction (x, y, z) into `first_interaction_r_z`."
@inline function fill_cut2_first_interaction!(ch::CutHistograms,
                                                x::Float64, y::Float64, z::Float64)
    r = sqrt(x*x + y*y)
    ir = _bin_idx(r, 0.0, ch.r_max_cm, ch.r_n_bins)
    iz = _bin_idx(z, ch.z_min_cm, ch.z_max_cm, ch.z_n_bins)
    if ir > 0 && iz > 0
        ch.first_interaction_r_z[ir, iz] += 1
    end
    nothing
end

# --- Cut 3 ---------------------------------------------------------------

"""
    fill_cut3!(ch, clusters, params)

Fill cut-3 histograms for one event whose outcome ∈ {`MS_rejected`,
`SS_outside_ROI`, `SS_in_ROI`, `companion_vetoed`} (i.e., events that
passed cut 2). Iterates only over visible clusters
(`c.ec > params.E_visible_keV / 1000`):

  - `dz_inclusive`: |Δz| between consecutive visible clusters in z-sort
    order, one entry per consecutive pair.
  - `n_visible`: count of visible clusters in the event.
  - `E_total`: Σ ec over visible clusters.
"""
function fill_cut3!(ch::CutHistograms,
                     clusters::Vector{Cluster},
                     params::MCParams)
    thr_MeV = params.E_visible_keV / 1000.0

    # Collect visible-cluster z and energies in z-sort order.
    vis_z = Float64[]
    E_tot = 0.0
    @inbounds for c in clusters
        if c.ec > thr_MeV
            push!(vis_z, c.zc)
            E_tot += c.ec
        end
    end
    n_vis = length(vis_z)

    _bin_int!(ch.n_visible, n_vis)
    if n_vis > 0
        iE = _bin_idx(E_tot, 0.0, ch.E_max_MeV, ch.E_n_bins)
        iE > 0 && (ch.E_total[iE] += 1)
    end
    if n_vis >= 2
        sort!(vis_z)
        @inbounds for i in 1:(n_vis - 1)
            adz = abs(vis_z[i+1] - vis_z[i])
            iz  = _bin_idx(adz, 0.0, ch.dz_max_cm, ch.dz_n_bins)
            iz > 0 && (ch.dz_inclusive[iz] += 1)
        end
    end
    nothing
end

# --- Cut 4 ---------------------------------------------------------------

"""
    fill_cut4!(ch, clusters, params)

Fill cut-4 histograms for one event whose outcome ∈ {`SS_in_ROI`,
`SS_outside_ROI`} — i.e., a single visible cluster, in FV, about to
have the ROI cut applied. Picks the (single) visible cluster.
"""
function fill_cut4!(ch::CutHistograms,
                     clusters::Vector{Cluster},
                     params::MCParams)
    thr_MeV = params.E_visible_keV / 1000.0
    @inbounds for c in clusters
        if c.ec > thr_MeV
            iec = _bin_idx(c.ec, 0.0, ch.E_max_MeV, ch.E_n_bins)
            ies = _bin_idx(c.es, 0.0, ch.E_max_MeV, ch.E_n_bins)
            iec > 0 && (ch.ss_ec_pre_roi[iec] += 1)
            ies > 0 && (ch.ss_es_pre_roi[ies] += 1)

            r  = sqrt(c.xc * c.xc + c.yc * c.yc)
            ir = _bin_idx(r,    0.0,           ch.r_max_cm, ch.r_n_bins)
            iz = _bin_idx(c.zc, ch.z_min_cm,   ch.z_max_cm, ch.z_n_bins)
            if ir > 0 && iz > 0
                ch.ss_r_z[ir, iz] += 1
            end
            return  # only one visible cluster by construction
        end
    end
    nothing
end

# --- Merge ---------------------------------------------------------------

"Sum `src` into `into` element-wise across every histogram field."
function merge_cut_histograms!(into::CutHistograms, src::CutHistograms)
    @assert into.u_n_bins  == src.u_n_bins
    @assert into.r_n_bins  == src.r_n_bins
    @assert into.z_n_bins  == src.z_n_bins
    @assert into.dz_n_bins == src.dz_n_bins
    @assert into.E_n_bins  == src.E_n_bins
    @. into.h_u_sampled              += src.h_u_sampled
    @. into.first_interaction_r_z    += src.first_interaction_r_z
    @. into.dz_inclusive             += src.dz_inclusive
    @. into.n_visible                += src.n_visible
    @. into.E_total                  += src.E_total
    @. into.ss_ec_pre_roi            += src.ss_ec_pre_roi
    @. into.ss_es_pre_roi            += src.ss_es_pre_roi
    @. into.ss_r_z                   += src.ss_r_z
    into
end
