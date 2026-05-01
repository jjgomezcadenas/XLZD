# src3/select.jl — Pure boolean predicates over already-built clusters.
#
# Two atomic decisions used by the upcoming event-outcome classifier:
#
#   * `select_SC(clusters)`        — "is this a single-cluster event?"
#   * `select_ROI(cluster, params)` — "does the SMEARED cluster energy
#                                      fall inside the ROI window?"
#
# Both are pure (no RNG, no global state). `select_ROI` reads
# `cluster.es` (already smeared by `build_clusters`); it does not
# resample. This makes the cluster the single source of truth for the
# measured energy and lets histograms plot `cluster.es` directly without
# re-smearing.
#
# These functions deliberately do NOT decide the per-event outcome
# (:SS_in_ROI / :MS_rejected / :skin_vetoed / etc.). The outcome
# classifier (`classify_event`, added later) composes these predicates
# with the FV check and rejection state to produce the final symbol.
#
# Pre-refactor counterpart `classify_ss_energy` (in mc_params.jl) returns
# Symbols and resamples internally; it dies at cutover.

"""
    select_SC(clusters::Vector{Cluster}) -> Bool

True iff the event has exactly one cluster. Pure cluster-count predicate.
Does not consider FV, skin overflow, or energy. Stateless.
"""
@inline select_SC(clusters::Vector{Cluster})::Bool = length(clusters) == 1

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
