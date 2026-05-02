# src3/select.jl — Pure boolean predicates over event-level data.
#
# Atomic decisions used by the event-outcome classifier:
#
#   * `select_SC(clusters, params)` — "is this a single-VISIBLE-cluster
#                                       event?" (sub-visible clusters
#                                       ignored, like in `select_FV`)
#   * `select_ROI(cluster, params)` — "does the SMEARED cluster energy
#                                      fall inside the ROI window?"
#   * `select_FV(clusters, params)` — "do all clusters above the visible
#                                      threshold lie inside the FV box?"
#                                      (slow per-cluster FV check)
#   * `select_skin(stack, params)`  — "is the cumulative skin energy
#                                      below the veto threshold?"
#                                      (slow cumulative skin check)
#
# All four are pure (no RNG, no global state). `select_ROI` reads
# `cluster.es` (already smeared by `build_clusters`); it does not
# resample.
#
# These functions deliberately do NOT decide the per-event outcome
# (:SS_in_ROI / :MS_rejected / :skin_vetoed / etc.). The outcome
# classifier `classify_event` composes them with the tracker's
# termination status to produce the final symbol.
#
# Pre-refactor counterpart `classify_ss_energy` (in mc_params.jl) returns
# Symbols and resamples internally; it dies at cutover.

"""
    select_SC(clusters::Vector{Cluster}, params::MCParams) -> Bool

True iff the event has exactly one **visible** cluster, where "visible"
means `c.ec * 1000 > params.E_visible_keV`. Sub-visible clusters are
ignored (consistent with `select_FV`, which also ignores them).

Why visibility matters: `build_clusters` includes every `:TPC` row with
`edep > 0`, so a sub-visible Compton kick or below-threshold residual
deposit creates its own cluster if it is z-separated from the main
deposit. Counting those would penalize vertical photon trajectories
(where sub-visible deposits z-spread by ~mfp) far more than horizontal
ones — even though the sub-visible deposits would never be seen by the
detector. Counting only visible clusters removes that bias.
"""
@inline function select_SC(clusters::Vector{Cluster}, params::MCParams)::Bool
    thr_MeV = params.E_visible_keV / 1000.0
    n = 0
    @inbounds for c in clusters
        if c.ec > thr_MeV
            n += 1
            n > 1 && return false
        end
    end
    return n == 1
end

"""
    select_ROI(cluster::Cluster, params::MCParams) -> Bool

Test whether the cluster's smeared energy `cluster.es` lies within
`±params.ROI_halfwidth_keV` of `params.Q_betabeta_keV`. Pure: no RNG,
no smearing — `cluster.es` is set once at cluster construction
(`build_clusters`) and consumed deterministically here.
"""
@inline function select_ROI(cluster::Cluster, params::MCParams)::Bool
    return abs(cluster.es * 1000.0 - params.Q_betabeta_keV) <= params.ROI_halfwidth_keV
end

"""
    select_FV(clusters::Vector{Cluster}, params::MCParams) -> Bool

True iff every cluster with `ec * 1000 > params.E_visible_keV` has its
energy-weighted centroid inside the FV box (`in_fv`). Sub-visible
clusters are ignored — they cannot drive a fiducial-volume rejection.

Empty `clusters` returns `true` (no offending data).

Stronger than the per-event SS-only FV check: catches MS events where
any of the multiple clusters is outside the FV box.
"""
function select_FV(clusters::Vector{Cluster}, params::MCParams)::Bool
    threshold_MeV = params.E_visible_keV / 1000.0
    @inbounds for c in clusters
        if c.ec > threshold_MeV && !in_fv(c.xc, c.yc, c.zc, params)
            return false
        end
    end
    return true
end

"""
    select_skin(stack::PhotonStack, params::MCParams) -> Bool

True iff the cumulative energy deposited in the LXe skin
(Σ `edep` over rows with `region === :Skin`), expressed in keV, is
strictly less than `params.E_skin_veto_keV`.

Empty stack (or stack with no `:Skin` rows) returns `true`.
"""
function select_skin(stack::PhotonStack, params::MCParams)::Bool
    skin_E_MeV = 0.0
    @inbounds for r in stack.rows
        if r.region === :Skin
            skin_E_MeV += r.edep
        end
    end
    return skin_E_MeV * 1000.0 < params.E_skin_veto_keV
end
