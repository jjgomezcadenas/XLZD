# src3/clusters.jl — Cluster type and cluster-builder.
#
# Defines:
#   * `Cluster`              — z-cluster data type (energy-weighted
#                               centroid, true energy `ec`, smeared `es`).
#   * `build_clusters(rng,
#                      stack, params)` — group :TPC rows of a PhotonStack
#                              into z-clusters; smear each cluster's
#                              energy at construction.
#
# Note: this file does not classify the event. The event-outcome
# classifier `classify_event` lives in `src3/classify.jl`.

using Random: AbstractRNG, randn

"""
    Cluster

A z-cluster of energy depositions in the active LXe region:

  * `xc, yc, zc` — energy-weighted centroid (cm)
  * `ec`        — true (deterministic) cluster energy = Σ edep (MeV)
  * `es`        — smeared (measured) cluster energy (MeV). One realization
                  of `ec + σ_E·ξ`, where σ_E = `params.σ_E_over_E · ec`
                  and ξ ~ N(0,1). Filled by `build_clusters`.
"""
struct Cluster
    xc::Float64
    yc::Float64
    zc::Float64
    ec::Float64
    es::Float64
end

"""
    build_clusters(rng::AbstractRNG, stack::PhotonStack, params::MCParams)
        -> Vector{Cluster}

Group `:TPC` rows of `stack` into z-clusters: sort by z, group
consecutive rows with `Δz < params.Δz_threshold_mm`, compute energy-
weighted (x, y, z) centroid and total energy `ec` per cluster, and
sample a smeared energy `es = ec + σ_E·randn(rng)` with
`σ_E = params.σ_E_over_E · ec`.

Filter rule: `region === :TPC && edep > 0`. This includes
`INT_BELOW_THRESH` rows in `:TPC` (residual energy at threshold) and
excludes `:Skin`, `:Inert`, `:Gas`, and any source-region rows.

Note on naming: this function only *builds* clusters; it does not
classify the event. The event-outcome classifier is `classify_event`
(in `src3/classify.jl`), which decides :SS_in_ROI / :MS_rejected /
:skin_vetoed / etc. from clusters and rejection state.
"""
function build_clusters(rng::AbstractRNG, stack::PhotonStack,
                         params::MCParams)::Vector{Cluster}
    actives = [r for r in stack.rows if r.region === :TPC && r.edep > 0]
    n = length(actives)
    n == 0 && return Cluster[]

    sorted = sort(actives; by = r -> r.z)
    Δz_thresh = Δz_threshold_cm(params)

    @inline function smear(ec::Float64)::Float64
        ec + (params.σ_E_over_E * ec) * randn(rng)
    end

    out = Cluster[]
    cum_E  = sorted[1].edep
    cum_xE = sorted[1].x * sorted[1].edep
    cum_yE = sorted[1].y * sorted[1].edep
    cum_zE = sorted[1].z * sorted[1].edep
    z_prev = sorted[1].z
    @inbounds for i in 2:n
        z_i = sorted[i].z
        if z_i - z_prev < Δz_thresh
            cum_E  += sorted[i].edep
            cum_xE += sorted[i].x * sorted[i].edep
            cum_yE += sorted[i].y * sorted[i].edep
            cum_zE += sorted[i].z * sorted[i].edep
        else
            push!(out, Cluster(cum_xE/cum_E, cum_yE/cum_E, cum_zE/cum_E,
                               cum_E, smear(cum_E)))
            cum_E  = sorted[i].edep
            cum_xE = sorted[i].x * sorted[i].edep
            cum_yE = sorted[i].y * sorted[i].edep
            cum_zE = sorted[i].z * sorted[i].edep
        end
        z_prev = z_i
    end
    push!(out, Cluster(cum_xE/cum_E, cum_yE/cum_E, cum_zE/cum_E,
                       cum_E, smear(cum_E)))
    out
end
