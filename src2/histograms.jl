# src2/histograms.jl — Per-event control histograms accumulated by the MC.
#
# Six histograms describe the visible-deposit structure of each photon:
#   1. SS / MS / no_cluster bucket counts
#   2. Δz from the first :active deposit to every subsequent :active deposit
#   3. Energy of the first :active deposit
#   4. Energy of each cluster (group of :active deposits with adjacent Δz < 3 mm)
#   5. Number of clusters per photon
#   6. Number of "extra" clusters (= n_clusters − 1)  [item 6 = option (a)]
#
# Only :active deposits feed into these histograms — same definition the
# main SS/MS outcome logic uses. Deposits in :skin or :inert are ignored
# here. Pair production and other events with no :active deposits go into
# the `no_cluster` bucket (and N_clusters = 0).

# ---------------------------------------------------------------------------
# Per-photon scratch buffer (reused across photons within a thread)
# ---------------------------------------------------------------------------

"""
    PhotonScratch

Per-thread scratch buffer that records (z, E_dep) for every visible
:active deposit during one photon's tracking. Reset to empty before each
photon by `track_one_photon!`. Pre-allocate once per thread.
"""
mutable struct PhotonScratch
    deposits::Vector{Tuple{Float64, Float64}}   # (z_cm, E_dep_MeV)
end
PhotonScratch() = PhotonScratch(Tuple{Float64, Float64}[])

# ---------------------------------------------------------------------------
# HistogramSet
# ---------------------------------------------------------------------------

"""
    HistogramSet

Six control histograms. All bins are uniform within their stated ranges.
Out-of-range fills are silently dropped (counted neither in nor out).

Defaults:
  - Δz: 100 bins on [0, 50] cm  (0.5 cm/bin)
  - E_first / E_cluster: 270 bins on [0, 2.7] MeV  (10 keV/bin)
  - N_clusters / N_extra: integer bins 0..N_max (default N_max = 20)
"""
struct HistogramSet
    # SS / MS / no_cluster
    ssms_counts::Vector{Int}                    # length 3

    # Δz from first to other interactions (cm)
    Δz_n_bins::Int
    Δz_max_cm::Float64
    Δz_counts::Vector{Int}

    # Energies (MeV)
    E_n_bins::Int
    E_max_MeV::Float64
    E_first_counts::Vector{Int}
    E_cluster_counts::Vector{Int}

    # Cluster multiplicities (integer bins 0 .. N_max)
    N_max::Int
    N_clusters_counts::Vector{Int}
    N_extra_counts::Vector{Int}
end

function HistogramSet(; Δz_n_bins::Int=100, Δz_max_cm::Real=50.0,
                       E_n_bins::Int=270,  E_max_MeV::Real=2.7,
                       N_max::Int=20)
    HistogramSet(
        zeros(Int, 3),
        Δz_n_bins, Float64(Δz_max_cm), zeros(Int, Δz_n_bins),
        E_n_bins,  Float64(E_max_MeV),
        zeros(Int, E_n_bins),
        zeros(Int, E_n_bins),
        N_max,
        zeros(Int, N_max + 1),
        zeros(Int, N_max + 1),
    )
end

# ---------------------------------------------------------------------------
# Filling helpers
# ---------------------------------------------------------------------------

"Bin index in 1..n for `x ∈ [lo, hi)`, else 0."
@inline function _bin_idx(x::Real, lo::Real, hi::Real, n::Int)::Int
    x < lo && return 0
    x >= hi && return 0
    return clamp(floor(Int, (x - lo) / (hi - lo) * n) + 1, 1, n)
end

@inline function fill_Δz!(h::HistogramSet, Δz_cm::Float64)
    i = _bin_idx(Δz_cm, 0.0, h.Δz_max_cm, h.Δz_n_bins)
    i > 0 && (h.Δz_counts[i] += 1)
    nothing
end

@inline function fill_E_first!(h::HistogramSet, E_MeV::Float64)
    i = _bin_idx(E_MeV, 0.0, h.E_max_MeV, h.E_n_bins)
    i > 0 && (h.E_first_counts[i] += 1)
    nothing
end

@inline function fill_E_cluster!(h::HistogramSet, E_MeV::Float64)
    i = _bin_idx(E_MeV, 0.0, h.E_max_MeV, h.E_n_bins)
    i > 0 && (h.E_cluster_counts[i] += 1)
    nothing
end

@inline function fill_N_clusters!(h::HistogramSet, n::Int)
    i = clamp(n, 0, h.N_max) + 1
    h.N_clusters_counts[i] += 1
    nothing
end

@inline function fill_N_extra!(h::HistogramSet, n::Int)
    i = clamp(n, 0, h.N_max) + 1
    h.N_extra_counts[i] += 1
    nothing
end

@inline function fill_ssms!(h::HistogramSet, kind::Symbol)
    i = kind === :SS ? 1 :
        kind === :MS ? 2 : 3   # default :no_cluster
    h.ssms_counts[i] += 1
    nothing
end

# ---------------------------------------------------------------------------
# Per-photon update from scratch buffer
# ---------------------------------------------------------------------------

"""
    update_histograms!(h, scratch, params)

Compute clusters from `scratch.deposits` (groups of consecutive entries
with Δz < `params.Δz_threshold_mm` after sorting by z) and fill all six
histograms in `h`.
"""
function update_histograms!(h::HistogramSet, scratch::PhotonScratch,
                             params::MCParams)
    n = length(scratch.deposits)
    if n == 0
        fill_ssms!(h, :no_cluster)
        fill_N_clusters!(h, 0)
        fill_N_extra!(h, 0)
        return
    end

    # Items 2 and 3: use tracking order (first deposit chronologically)
    z_first, E_first = scratch.deposits[1]
    fill_E_first!(h, E_first)
    @inbounds for i in 2:n
        Δz = abs(scratch.deposits[i][1] - z_first)
        fill_Δz!(h, Δz)
    end

    # Items 4–6: sort by z and group into clusters by Δz < threshold
    sorted = sort(scratch.deposits; by = x -> x[1])
    Δz_thresh = Δz_threshold_cm(params)

    cluster_E  = sorted[1][2]
    z_prev     = sorted[1][1]
    n_clusters = 1
    @inbounds for i in 2:n
        z_i = sorted[i][1]
        if z_i - z_prev < Δz_thresh
            cluster_E += sorted[i][2]
        else
            fill_E_cluster!(h, cluster_E)
            cluster_E = sorted[i][2]
            n_clusters += 1
        end
        z_prev = z_i
    end
    fill_E_cluster!(h, cluster_E)         # last cluster

    fill_N_clusters!(h, n_clusters)
    fill_N_extra!(h, n_clusters - 1)
    fill_ssms!(h, n_clusters == 1 ? :SS : :MS)
    nothing
end

# ---------------------------------------------------------------------------
# Merge (across threads)
# ---------------------------------------------------------------------------

"""
    merge_histograms!(into::HistogramSet, src::HistogramSet)

In-place sum of `src` into `into`. Both must have matching binning.
"""
function merge_histograms!(into::HistogramSet, src::HistogramSet)
    @assert into.Δz_n_bins == src.Δz_n_bins
    @assert into.E_n_bins  == src.E_n_bins
    @assert into.N_max     == src.N_max
    @. into.ssms_counts       += src.ssms_counts
    @. into.Δz_counts         += src.Δz_counts
    @. into.E_first_counts    += src.E_first_counts
    @. into.E_cluster_counts  += src.E_cluster_counts
    @. into.N_clusters_counts += src.N_clusters_counts
    @. into.N_extra_counts    += src.N_extra_counts
    into
end
