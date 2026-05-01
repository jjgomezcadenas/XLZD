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
end

function ClusterHistogramSet(;
        E_n_bins::Int=270, E_max_MeV::Real=2.7,
        N_clusters_n::Int=21,
        D3_n_bins::Int=100, D3_max_cm::Real=200.0,
        dz_n_bins::Int=100, dz_max_cm::Real=50.0,
        r2_n_bins::Int=100, r2_max_cm2::Real=82.1^2,
        z_n_bins::Int=100,  z_min_cm::Real=-69.0, z_max_cm::Real=145.6,
        D_n_bins::Int=100,  D_max_cm::Real=100.0)
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
    into
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
