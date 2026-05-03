# src3/histograms.jl — Diagnostic histograms for the per-photon MC.
#
# Three histogram families:
#
#   * `StackHistogramSet`     — per-event quantities derived from the
#                                full PhotonStack (chain depth, first
#                                interaction type, per-region deposit
#                                counts, inclusive energy, etc.)
#   * `ClusterHistogramSet`   — per-event quantities derived from the
#                                cluster vector (cluster energies,
#                                pair distances, r²-vs-z and D-vs-z
#                                heatmaps).
#   * `RejectionHistograms`   — diagnostic for fast_veto rejections
#                                (skin overflow / outside-FV first
#                                deposit). 2D r²×z + 1D E for each.
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

Per-event histograms accumulated from the rows of `PhotonStack`.

Bins (defaults; configurable via kwargs):
  - `ng_max`              21 bins (0..20)        chain depth
  - `first_interaction`    4 bins                PHOTO/COMPTON/PAIR/BELOW_THRESH
  - `n_photo`             21 bins (0..20)        per-event counts
  - `n_compton`           21 bins (0..20)
  - `n_pair`              11 bins (0..10)
  - `n_below_thresh`      11 bins (0..10)
  - `inclusive_edep`     270 bins on [0, 2.7] MeV
  - `E_first`            270 bins on [0, 2.7] MeV  (first :TPC deposit)
  - `Δz`                 100 bins on [0, 50] cm    (Δz of subsequent
                                                    :TPC deposits from
                                                    first :TPC deposit)
  - `region_interaction` 3×4 matrix (Symbol order documented below)

`region_interaction_counts[r, t]`:
  rows: 1=:TPC, 2=:Skin, 3=:Inert
  cols: 1=PHOTO, 2=COMPTON, 3=PAIR, 4=BELOW_THRESH
"""
struct StackHistogramSet
    # Integer-bucket counts
    ng_max_counts::Vector{Int}
    first_interaction_counts::Vector{Int}      # length 4
    n_photo_counts::Vector{Int}
    n_compton_counts::Vector{Int}
    n_pair_counts::Vector{Int}
    n_below_thresh_counts::Vector{Int}

    # Energy histograms
    E_n_bins::Int
    E_max_MeV::Float64
    inclusive_edep_counts::Vector{Int}
    E_first_counts::Vector{Int}

    # Spatial: Δz from first :TPC deposit to subsequent :TPC deposits
    Δz_n_bins::Int
    Δz_max_cm::Float64
    Δz_counts::Vector{Int}

    # Per-region × per-interaction counts (3×4, see docstring)
    region_interaction_counts::Matrix{Int}
end

const _STACK_REGIONS      = (:TPC, :Skin, :Inert)
const _STACK_INTERACTIONS = (INT_PHOTO, INT_COMPTON, INT_PAIR, INT_BELOW_THRESH)

function StackHistogramSet(;
        ng_max_n::Int=21,
        n_photo_n::Int=21, n_compton_n::Int=21,
        n_pair_n::Int=11,  n_below_thresh_n::Int=11,
        E_n_bins::Int=270, E_max_MeV::Real=2.7,
        Δz_n_bins::Int=100, Δz_max_cm::Real=50.0)
    StackHistogramSet(
        zeros(Int, ng_max_n),
        zeros(Int, 4),
        zeros(Int, n_photo_n),
        zeros(Int, n_compton_n),
        zeros(Int, n_pair_n),
        zeros(Int, n_below_thresh_n),
        E_n_bins, Float64(E_max_MeV),
        zeros(Int, E_n_bins),
        zeros(Int, E_n_bins),
        Δz_n_bins, Float64(Δz_max_cm),
        zeros(Int, Δz_n_bins),
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

Fold one event's stack into `sh`. Does nothing for an empty stack
(those events contribute to every per-event count as 0, which lands
in bin 0 of the integer-bucket histograms — i.e. the n_*=0 events).
"""
function update_stack_histograms!(sh::StackHistogramSet,
                                   stack::PhotonStack,
                                   params::MCParams)
    n = length(stack)

    # Chain depth — fill even when zero (records "no interaction" events).
    _bin_int!(sh.ng_max_counts, n)

    # Per-interaction-type counts (also filled when zero).
    n_photo = 0
    n_compton = 0
    n_pair = 0
    n_below = 0
    inclusive_E = 0.0
    @inbounds for r in stack.rows
        if r.interaction === INT_PHOTO
            n_photo += 1
        elseif r.interaction === INT_COMPTON
            n_compton += 1
        elseif r.interaction === INT_PAIR
            n_pair += 1
        elseif r.interaction === INT_BELOW_THRESH
            n_below += 1
        end
        inclusive_E += r.edep
        ri = _stack_region_row(r.region)
        ci = _interaction_col(r.interaction)
        if ri > 0 && ci > 0
            sh.region_interaction_counts[ri, ci] += 1
        end
    end
    _bin_int!(sh.n_photo_counts,         n_photo)
    _bin_int!(sh.n_compton_counts,       n_compton)
    _bin_int!(sh.n_pair_counts,          n_pair)
    _bin_int!(sh.n_below_thresh_counts,  n_below)

    # First-interaction-type bin (only when any rows).
    if n > 0
        ci = _interaction_col(stack.rows[1].interaction)
        ci > 0 && (sh.first_interaction_counts[ci] += 1)
    end

    # Inclusive edep — skip when zero to avoid clogging bin 1.
    if inclusive_E > 0.0
        ie = _bin_idx(inclusive_E, 0.0, sh.E_max_MeV, sh.E_n_bins)
        ie > 0 && (sh.inclusive_edep_counts[ie] += 1)
    end

    # First :TPC deposit: E_first and Δz to subsequent :TPC rows.
    first_tpc_idx = 0
    @inbounds for i in 1:n
        if stack.rows[i].region === :TPC
            first_tpc_idx = i
            break
        end
    end
    if first_tpc_idx > 0
        z_first = stack.rows[first_tpc_idx].z
        e_first = stack.rows[first_tpc_idx].edep
        ie = _bin_idx(e_first, 0.0, sh.E_max_MeV, sh.E_n_bins)
        ie > 0 && (sh.E_first_counts[ie] += 1)
        @inbounds for i in (first_tpc_idx + 1):n
            if stack.rows[i].region === :TPC
                Δz = abs(stack.rows[i].z - z_first)
                iz = _bin_idx(Δz, 0.0, sh.Δz_max_cm, sh.Δz_n_bins)
                iz > 0 && (sh.Δz_counts[iz] += 1)
            end
        end
    end
    nothing
end

function merge_stack_histograms!(into::StackHistogramSet,
                                  src::StackHistogramSet)
    @assert into.E_n_bins == src.E_n_bins
    @assert into.Δz_n_bins == src.Δz_n_bins
    @. into.ng_max_counts            += src.ng_max_counts
    @. into.first_interaction_counts += src.first_interaction_counts
    @. into.n_photo_counts           += src.n_photo_counts
    @. into.n_compton_counts         += src.n_compton_counts
    @. into.n_pair_counts            += src.n_pair_counts
    @. into.n_below_thresh_counts    += src.n_below_thresh_counts
    @. into.inclusive_edep_counts    += src.inclusive_edep_counts
    @. into.E_first_counts           += src.E_first_counts
    @. into.Δz_counts                += src.Δz_counts
    @. into.region_interaction_counts += src.region_interaction_counts
    into
end

# ===========================================================================
# ClusterHistogramSet — per-event quantities from the cluster vector.
# ===========================================================================
"""
    ClusterHistogramSet

Per-event histograms accumulated from the `Vector{Cluster}` produced by
`build_clusters`.

Bins (defaults; configurable via kwargs):
  - `Ec`              270 bins on [0, 2.7] MeV   per-cluster energy
  - `N_clusters`       21 bins (0..20)            cluster multiplicity
  - `Emax`, `Emin`,
    `Einc`            270 bins on [0, 2.7] MeV   per-event extremes/sum
  - `closest_D3`,
    `furthest_D3`     100 bins on [0, 200] cm    3D nearest/farthest pair
  - `closest_dz`,
    `furthest_dz`     100 bins on [0, 50]  cm    z-only pair distances
  - `r2_vs_z`        100×100 on [0, R_ICV²] × [z_min, z_max]
  - `D_vs_z`         100×100 on [0,  100  ] × [z_min, z_max]
"""
struct ClusterHistogramSet
    # Energies
    E_n_bins::Int
    E_max_MeV::Float64
    Ec_counts::Vector{Int}
    N_clusters_counts::Vector{Int}
    Emax_counts::Vector{Int}
    Emin_counts::Vector{Int}
    Einc_counts::Vector{Int}

    # 1D pair distances
    D3_n_bins::Int
    D3_max_cm::Float64
    closest_D3_counts::Vector{Int}
    furthest_D3_counts::Vector{Int}
    dz_n_bins::Int
    dz_max_cm::Float64
    closest_dz_counts::Vector{Int}
    furthest_dz_counts::Vector{Int}

    # 2D heatmaps
    r2_n_bins::Int
    r2_max_cm2::Float64
    z_n_bins::Int
    z_min_cm::Float64
    z_max_cm::Float64
    r2_vs_z_2d_counts::Matrix{Int}

    D_n_bins::Int
    D_max_cm::Float64
    D_vs_z_2d_counts::Matrix{Int}

    # Per-cluster (ec, Δz_to_nearest_other_cluster) — 2D scatter
    # for diagnosing MS-rejection regime.
    ec2d_n_bins::Int
    ec2d_max_MeV::Float64
    Ec_vs_dz_2d_counts::Matrix{Int}    # shape: (dz_n_bins, ec2d_n_bins)

    # SS pre-ROI spectra: cluster energy of the (single) visible cluster
    # for events that have passed all upstream cuts (skin, FV, SC) and are
    # about to be ROI-tested. `ec` is the true (deterministic) cluster
    # energy, `es` is the smeared energy that the ROI cut acts on. Filled
    # by `fill_ss_pre_roi!` from `run.jl`. Same binning as `Ec_counts`.
    ss_ec_pre_roi_counts::Vector{Int}
    ss_es_pre_roi_counts::Vector{Int}
end

function ClusterHistogramSet(;
        E_n_bins::Int=270, E_max_MeV::Real=2.7,
        N_clusters_n::Int=21,
        D3_n_bins::Int=100, D3_max_cm::Real=200.0,
        dz_n_bins::Int=100, dz_max_cm::Real=50.0,
        r2_n_bins::Int=100, r2_max_cm2::Real=82.1^2,
        z_n_bins::Int=100,  z_min_cm::Real=-69.0, z_max_cm::Real=145.6,
        D_n_bins::Int=100,  D_max_cm::Real=100.0,
        ec2d_n_bins::Int=100, ec2d_max_MeV::Real=2.7)
    ClusterHistogramSet(
        E_n_bins, Float64(E_max_MeV),
        zeros(Int, E_n_bins),
        zeros(Int, N_clusters_n),
        zeros(Int, E_n_bins),
        zeros(Int, E_n_bins),
        zeros(Int, E_n_bins),
        D3_n_bins, Float64(D3_max_cm),
        zeros(Int, D3_n_bins),
        zeros(Int, D3_n_bins),
        dz_n_bins, Float64(dz_max_cm),
        zeros(Int, dz_n_bins),
        zeros(Int, dz_n_bins),
        r2_n_bins, Float64(r2_max_cm2),
        z_n_bins,  Float64(z_min_cm), Float64(z_max_cm),
        zeros(Int, r2_n_bins, z_n_bins),
        D_n_bins,  Float64(D_max_cm),
        zeros(Int, D_n_bins, z_n_bins),
        ec2d_n_bins, Float64(ec2d_max_MeV),
        zeros(Int, dz_n_bins, ec2d_n_bins),
        zeros(Int, E_n_bins),
        zeros(Int, E_n_bins),
    )
end

"""
    update_cluster_histograms!(ch::ClusterHistogramSet,
                                clusters::Vector{Cluster}, params::MCParams)

Fold one event's clusters into `ch`. Empty cluster vectors only update
`N_clusters_counts[1]` (= 0 clusters); per-event Emax/Emin/Einc and
pair-distance histograms are skipped (no clusters → nothing to fill).
"""
function update_cluster_histograms!(ch::ClusterHistogramSet,
                                     clusters::Vector{Cluster},
                                     params::MCParams)
    nc = length(clusters)
    _bin_int!(ch.N_clusters_counts, nc)
    nc == 0 && return

    # Per-cluster fills + per-event extremes/sum.
    Emax = -Inf
    Emin =  Inf
    Einc = 0.0
    @inbounds for c in clusters
        ie = _bin_idx(c.ec, 0.0, ch.E_max_MeV, ch.E_n_bins)
        ie > 0 && (ch.Ec_counts[ie] += 1)
        Einc += c.ec
        Emax = max(Emax, c.ec)
        Emin = min(Emin, c.ec)

        # 2D heatmaps: r² vs z and D = √(x²+y²) vs z.
        r2 = c.xc * c.xc + c.yc * c.yc
        D  = sqrt(r2)
        ir2 = _bin_idx(r2, 0.0, ch.r2_max_cm2, ch.r2_n_bins)
        iD  = _bin_idx(D,  0.0, ch.D_max_cm,   ch.D_n_bins)
        iz  = _bin_idx(c.zc, ch.z_min_cm, ch.z_max_cm, ch.z_n_bins)
        if iz > 0
            ir2 > 0 && (ch.r2_vs_z_2d_counts[ir2, iz] += 1)
            iD  > 0 && (ch.D_vs_z_2d_counts[iD,  iz] += 1)
        end
    end
    iemax = _bin_idx(Emax, 0.0, ch.E_max_MeV, ch.E_n_bins)
    iemin = _bin_idx(Emin, 0.0, ch.E_max_MeV, ch.E_n_bins)
    ieinc = _bin_idx(Einc, 0.0, ch.E_max_MeV, ch.E_n_bins)
    iemax > 0 && (ch.Emax_counts[iemax] += 1)
    iemin > 0 && (ch.Emin_counts[iemin] += 1)
    ieinc > 0 && (ch.Einc_counts[ieinc] += 1)

    # Pair distances: closest and furthest among all C(n,2) pairs.
    if nc >= 2
        d3_min = Inf;  d3_max = -Inf
        dz_min = Inf;  dz_max = -Inf
        @inbounds for i in 1:(nc - 1)
            ci = clusters[i]
            for j in (i + 1):nc
                cj = clusters[j]
                dx = ci.xc - cj.xc
                dy = ci.yc - cj.yc
                dz = ci.zc - cj.zc
                d3 = sqrt(dx*dx + dy*dy + dz*dz)
                adz = abs(dz)
                d3 < d3_min && (d3_min = d3)
                d3 > d3_max && (d3_max = d3)
                adz < dz_min && (dz_min = adz)
                adz > dz_max && (dz_max = adz)
            end
        end
        i = _bin_idx(d3_min, 0.0, ch.D3_max_cm, ch.D3_n_bins); i > 0 && (ch.closest_D3_counts[i]  += 1)
        i = _bin_idx(d3_max, 0.0, ch.D3_max_cm, ch.D3_n_bins); i > 0 && (ch.furthest_D3_counts[i] += 1)
        i = _bin_idx(dz_min, 0.0, ch.dz_max_cm, ch.dz_n_bins); i > 0 && (ch.closest_dz_counts[i]  += 1)
        i = _bin_idx(dz_max, 0.0, ch.dz_max_cm, ch.dz_n_bins); i > 0 && (ch.furthest_dz_counts[i] += 1)

        # Per-cluster (ec, dz_to_nearest_other) — diagnostic for MS regime.
        @inbounds for i in 1:nc
            ci = clusters[i]
            dz_near = Inf
            for j in 1:nc
                j == i && continue
                adz = abs(ci.zc - clusters[j].zc)
                adz < dz_near && (dz_near = adz)
            end
            idx_dz = _bin_idx(dz_near, 0.0, ch.dz_max_cm,  ch.dz_n_bins)
            idx_ec = _bin_idx(ci.ec,  0.0, ch.ec2d_max_MeV, ch.ec2d_n_bins)
            if idx_dz > 0 && idx_ec > 0
                ch.Ec_vs_dz_2d_counts[idx_dz, idx_ec] += 1
            end
        end
    end
    nothing
end

function merge_cluster_histograms!(into::ClusterHistogramSet,
                                    src::ClusterHistogramSet)
    @assert into.E_n_bins  == src.E_n_bins
    @assert into.D3_n_bins == src.D3_n_bins
    @assert into.dz_n_bins == src.dz_n_bins
    @assert into.r2_n_bins == src.r2_n_bins
    @assert into.z_n_bins  == src.z_n_bins
    @assert into.D_n_bins  == src.D_n_bins
    @. into.Ec_counts             += src.Ec_counts
    @. into.N_clusters_counts     += src.N_clusters_counts
    @. into.Emax_counts           += src.Emax_counts
    @. into.Emin_counts           += src.Emin_counts
    @. into.Einc_counts           += src.Einc_counts
    @. into.closest_D3_counts     += src.closest_D3_counts
    @. into.furthest_D3_counts    += src.furthest_D3_counts
    @. into.closest_dz_counts     += src.closest_dz_counts
    @. into.furthest_dz_counts    += src.furthest_dz_counts
    @. into.r2_vs_z_2d_counts     += src.r2_vs_z_2d_counts
    @. into.D_vs_z_2d_counts      += src.D_vs_z_2d_counts
    @. into.Ec_vs_dz_2d_counts    += src.Ec_vs_dz_2d_counts
    @. into.ss_ec_pre_roi_counts  += src.ss_ec_pre_roi_counts
    @. into.ss_es_pre_roi_counts  += src.ss_es_pre_roi_counts
    into
end

"""
    fill_ss_pre_roi!(ch::ClusterHistogramSet, clusters::Vector{Cluster},
                      params::MCParams)

Bin the (ec, es) of the single visible cluster into `ss_ec_pre_roi_counts`
and `ss_es_pre_roi_counts`. Caller must guarantee the event reached the
pre-ROI SS state — i.e., `classify_event` returned `:SS_in_ROI` or
`:SS_outside_ROI`. Picks the visible cluster (ec > E_visible_keV/1000)
rather than `clusters[1]` to handle the case where a sub-visible cluster
is sorted ahead of the visible one.
"""
function fill_ss_pre_roi!(ch::ClusterHistogramSet,
                           clusters::Vector{Cluster},
                           params::MCParams)
    thr_MeV = params.E_visible_keV / 1000.0
    @inbounds for c in clusters
        if c.ec > thr_MeV
            iec = _bin_idx(c.ec, 0.0, ch.E_max_MeV, ch.E_n_bins)
            ies = _bin_idx(c.es, 0.0, ch.E_max_MeV, ch.E_n_bins)
            iec > 0 && (ch.ss_ec_pre_roi_counts[iec] += 1)
            ies > 0 && (ch.ss_es_pre_roi_counts[ies] += 1)
            return
        end
    end
end

# ===========================================================================
# RejectionHistograms — fast_veto rejection diagnostics (unchanged).
# ===========================================================================
"""
    RejectionHistograms

Two pairs of histograms: one pair (2-D r²×z + 1-D E) for skin-rejection
events and one pair for FV-rejection events. Each histogram counts the
*triggering* deposit position / energy.
"""
struct RejectionHistograms
    r2_n_bins::Int
    r2_max_cm2::Float64
    z_n_bins::Int
    z_min_cm::Float64
    z_max_cm::Float64
    E_n_bins::Int
    E_max_MeV::Float64
    skin_r2z_counts::Matrix{Int}
    skin_E_counts::Vector{Int}
    fv_r2z_counts::Matrix{Int}
    fv_E_counts::Vector{Int}
end

function RejectionHistograms(; r2_n_bins::Int=60, r2_max_cm2::Real=82.1^2,
                              z_n_bins::Int=100, z_min_cm::Real=-69.0,
                              z_max_cm::Real=145.6,
                              E_n_bins::Int=270,
                              E_max_MeV::Real=2.7)
    RejectionHistograms(
        r2_n_bins, Float64(r2_max_cm2),
        z_n_bins,  Float64(z_min_cm), Float64(z_max_cm),
        E_n_bins,  Float64(E_max_MeV),
        zeros(Int, r2_n_bins, z_n_bins),
        zeros(Int, E_n_bins),
        zeros(Int, r2_n_bins, z_n_bins),
        zeros(Int, E_n_bins),
    )
end

@inline function _fill_r2z!(M::Matrix{Int}, x::Float64, y::Float64, z::Float64,
                             rh::RejectionHistograms)
    r2 = x*x + y*y
    i = _bin_idx(r2, 0.0, rh.r2_max_cm2, rh.r2_n_bins)
    j = _bin_idx(z,  rh.z_min_cm, rh.z_max_cm, rh.z_n_bins)
    if i > 0 && j > 0
        M[i, j] += 1
    end
    nothing
end

@inline function fill_rejected_skin!(rh::RejectionHistograms,
                                      x::Float64, y::Float64, z::Float64,
                                      E_dep::Float64)
    _fill_r2z!(rh.skin_r2z_counts, x, y, z, rh)
    iE = _bin_idx(E_dep, 0.0, rh.E_max_MeV, rh.E_n_bins)
    iE > 0 && (rh.skin_E_counts[iE] += 1)
    nothing
end

@inline function fill_rejected_fv!(rh::RejectionHistograms,
                                    x::Float64, y::Float64, z::Float64,
                                    E_dep::Float64)
    _fill_r2z!(rh.fv_r2z_counts, x, y, z, rh)
    iE = _bin_idx(E_dep, 0.0, rh.E_max_MeV, rh.E_n_bins)
    iE > 0 && (rh.fv_E_counts[iE] += 1)
    nothing
end

function merge_rejection_histograms!(into::RejectionHistograms,
                                      src::RejectionHistograms)
    @. into.skin_r2z_counts += src.skin_r2z_counts
    @. into.skin_E_counts   += src.skin_E_counts
    @. into.fv_r2z_counts   += src.fv_r2z_counts
    @. into.fv_E_counts     += src.fv_E_counts
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
