# src3/classify.jl — Cluster type and the functions that produce clusters.
#
# Two cluster-building functions live here:
#
#   * `compute_clusters(deposits, params)` — operates on the legacy
#     `Vector{LXeDeposit}` produced by the current src3 tracker
#     (`PhotonScratch.deposits`). Pre-refactor; will be deleted at cutover.
#
#   * `build_clusters(stack, params)` — operates on the new `PhotonStack`
#     produced by the upcoming stack-based tracker.
#
# Both produce a `Vector{Cluster}` and use the same z-only adjacency rule
# (`Δz < params.Δz_threshold_mm`) and energy-weighted centroid.
#
# Note: neither of these classifies the event. The event-outcome classifier
# (`classify_event`, added later) decides the per-event outcome symbol.

"""
    Cluster

A z-cluster of energy depositions in the active LXe region:
energy-weighted centroid (xc, yc, zc) and total energy `ec` (MeV).
"""
struct Cluster
    xc::Float64
    yc::Float64
    zc::Float64
    ec::Float64
end

"""
    compute_clusters(deposits::Vector{LXeDeposit}, params::MCParams) -> Vector{Cluster}

Group `:active` deposits in `deposits` (any region; non-`:active` ignored)
into z-clusters: sort by z, group consecutive entries with Δz <
`params.Δz_threshold_mm` after sorting, compute energy-weighted (x, y, z)
and total E per cluster.
"""
function compute_clusters(deposits::Vector{LXeDeposit},
                           params::MCParams)::Vector{Cluster}
    actives = [d for d in deposits if d.region === :active]
    n = length(actives)
    n == 0 && return Cluster[]

    sorted = sort(actives; by = d -> d.z)
    Δz_thresh = Δz_threshold_cm(params)

    out = Cluster[]
    cum_E = sorted[1].E_dep
    cum_xE = sorted[1].x * sorted[1].E_dep
    cum_yE = sorted[1].y * sorted[1].E_dep
    cum_zE = sorted[1].z * sorted[1].E_dep
    z_prev = sorted[1].z
    @inbounds for i in 2:n
        z_i = sorted[i].z
        if z_i - z_prev < Δz_thresh
            cum_E  += sorted[i].E_dep
            cum_xE += sorted[i].x * sorted[i].E_dep
            cum_yE += sorted[i].y * sorted[i].E_dep
            cum_zE += sorted[i].z * sorted[i].E_dep
        else
            push!(out, Cluster(cum_xE/cum_E, cum_yE/cum_E, cum_zE/cum_E, cum_E))
            cum_E  = sorted[i].E_dep
            cum_xE = sorted[i].x * sorted[i].E_dep
            cum_yE = sorted[i].y * sorted[i].E_dep
            cum_zE = sorted[i].z * sorted[i].E_dep
        end
        z_prev = z_i
    end
    push!(out, Cluster(cum_xE/cum_E, cum_yE/cum_E, cum_zE/cum_E, cum_E))
    out
end

"""
    build_clusters(stack::PhotonStack, params::MCParams) -> Vector{Cluster}

Build z-clusters from the rows of `stack`. Same algorithm as
`compute_clusters`, but the input is the new `PhotonStack` produced
by the stack-based tracker. Once the stack tracker replaces the
deposit-list tracker, `compute_clusters` will be deleted.

Filter rule: `region === :TPC && edep > 0`. This includes
`INT_BELOW_THRESH` rows in `:TPC` (whose `edep` is the residual energy
dumped at threshold) and excludes `:Skin`, `:Inert`, `:Gas`, and any
source-region rows.

Grouping: sort filtered rows by z, group consecutive rows with
`Δz < params.Δz_threshold_mm`. Centroid is energy-weighted; cluster
energy is the sum of `edep`.

Note on naming: this function only *builds* clusters; it does not
classify the event. The event-outcome classifier is `classify_event`
(added in a later step), which decides :SS_in_ROI / :MS_rejected /
:skin_vetoed / etc. from clusters and rejection state.
"""
function build_clusters(stack::PhotonStack,
                         params::MCParams)::Vector{Cluster}
    actives = [r for r in stack.rows if r.region === :TPC && r.edep > 0]
    n = length(actives)
    n == 0 && return Cluster[]

    sorted = sort(actives; by = r -> r.z)
    Δz_thresh = Δz_threshold_cm(params)

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
            push!(out, Cluster(cum_xE/cum_E, cum_yE/cum_E, cum_zE/cum_E, cum_E))
            cum_E  = sorted[i].edep
            cum_xE = sorted[i].x * sorted[i].edep
            cum_yE = sorted[i].y * sorted[i].edep
            cum_zE = sorted[i].z * sorted[i].edep
        end
        z_prev = z_i
    end
    push!(out, Cluster(cum_xE/cum_E, cum_yE/cum_E, cum_zE/cum_E, cum_E))
    out
end
