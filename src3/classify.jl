# src3/classify.jl — Event-outcome classifier.
#
# `classify_event` consumes the result of one photon's transport
# (the tracker's termination status, plus the clusters produced by
# `build_clusters`) and assigns the event to one of the outcome
# categories used downstream for counting and histograms.
#
# This is multi-category classification via a chain of selections:
# skin → cluster-count → FV → ROI. The atomic predicates live in
# `select.jl`; this file composes them with the tracker's termination
# status and the FV box check.
#
# Output alphabet (`CLASSIFY_EVENT_OUTCOMES`):
#   :escaped          — no LXe entry, or LXe entry but no visible clusters
#   :skin_vetoed      — photon deposited > E_skin_veto_keV in :Skin
#                       (tracker's early-skin-reject decision)
#   :MS_rejected      — multiple clusters
#   :SS_outside_FV    — single cluster (or first deposit) outside FV box
#   :SS_outside_ROI   — single cluster, in FV, smeared E outside ROI window
#   :SS_in_ROI        — single cluster, in FV, smeared E inside ROI window
#
# Note: `:companion_vetoed` is NOT produced by this function. It is
# applied later in `run_mc` for Tl-208 events whose cascade companion
# γ produces a visible deposit (post-classification veto).
#
# Tracker termination-status alphabet (`TRACK_STATUSES`) — returned by
# `track_photon_stack` and consumed here:
#   :completed     — tracker ran to completion; use the populated stack
#   :escaped       — photon left LXe with no visible deposit
#   :vetoed_skin   — skin-energy threshold exceeded; tracker bailed out
#   :rejected_fv   — first :TPC deposit was outside FV; tracker bailed out

const CLASSIFY_EVENT_OUTCOMES = (
    :escaped, :MS_rejected, :skin_vetoed,
    :SS_outside_FV, :SS_outside_ROI, :SS_in_ROI,
)

const TRACK_STATUSES = (
    :completed, :escaped, :vetoed_skin, :rejected_fv,
)

"""
    classify_event(status::Symbol, clusters::Vector{Cluster},
                    params::MCParams) -> Symbol

Assign one of the symbols in `CLASSIFY_EVENT_OUTCOMES` (excluding
`:companion_vetoed`, which is applied later in `run_mc`) to the event.

`status` is the tracker's termination status (one of `TRACK_STATUSES`);
`clusters` is the output of `build_clusters` (smearing already done).

Decision chain (first match wins):

    vetoed_skin    → :skin_vetoed
    rejected_fv    → :SS_outside_FV
    escaped        → :escaped
    empty clusters → :escaped
    >1 cluster     → :MS_rejected
    1 cluster:
        cluster outside FV → :SS_outside_FV
        else if in ROI     → :SS_in_ROI
        else               → :SS_outside_ROI

Pure: no RNG (smearing already happened in `build_clusters`), no
detector geometry beyond what `params` carries (FV box and ROI window).
"""
function classify_event(status::Symbol,
                         clusters::Vector{Cluster},
                         params::MCParams)::Symbol
    status === :vetoed_skin   && return :skin_vetoed
    status === :rejected_fv   && return :SS_outside_FV
    status === :escaped       && return :escaped
    isempty(clusters)         && return :escaped
    if !select_SC(clusters)
        return :MS_rejected
    end
    c = clusters[1]
    if !in_fv(c.xc, c.yc, c.zc, params)
        return :SS_outside_FV
    end
    return select_ROI(c, params) ? :SS_in_ROI : :SS_outside_ROI
end
