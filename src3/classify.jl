# src3/classify.jl — Event-outcome classifier.
#
# `classify_event` consumes the result of one photon's transport
# (the tracker's termination status, the populated stack, and the
# clusters produced by `build_clusters`) and assigns the event to
# one of the outcome categories used downstream for counting and
# histograms.
#
# This is multi-category classification via a chain of selections:
#   skin → cluster-count → FV → ROI
#
# Two layers of skin/FV checks compose:
#   - FAST tracker-level (in `fast_veto`): single-deposit check on the
#     FIRST interaction, returns :vetoed_skin / :rejected_fv early.
#   - SLOW classifier-level (here): cumulative skin (`select_skin`) and
#     per-cluster FV (`select_FV`) checks across the whole event.
#
# Output alphabet (`CLASSIFY_EVENT_OUTCOMES`):
#   :escaped         — no LXe entry, or LXe entry but no visible clusters
#   :skin_vetoed     — skin energy ≥ E_skin_veto_keV (fast or cumulative)
#   :MS_rejected     — multiple clusters (and not already skin/FV-rejected)
#   :outside_FV      — any cluster with ec > E_visible_keV outside the FV
#                      (fast first-deposit OR slow per-cluster check)
#   :SS_outside_ROI  — single cluster, FV OK, smeared E outside ROI window
#   :SS_in_ROI       — single cluster, FV OK, smeared E inside ROI window
#
# Note: `:companion_vetoed` is NOT produced by this function. It is
# applied later in `run_mc` for Tl-208 events whose cascade companion
# γ produces a visible deposit (post-classification veto).
#
# Tracker termination-status alphabet (`TRACK_STATUSES`) — returned by
# `track_photon_stack` and consumed here:
#   :completed     — tracker ran to completion; use the populated stack
#   :escaped       — photon left LXe with no visible deposit
#   :vetoed_skin   — fast skin-veto fired
#   :rejected_fv   — fast FV-reject fired

const CLASSIFY_EVENT_OUTCOMES = (
    :escaped, :MS_rejected, :skin_vetoed,
    :outside_FV, :SS_outside_ROI, :SS_in_ROI,
)

const TRACK_STATUSES = (
    :completed, :escaped, :vetoed_skin, :rejected_fv,
)

"""
    classify_event(status::Symbol, stack::PhotonStack,
                    clusters::Vector{Cluster}, params::MCParams) -> Symbol

Assign one of the symbols in `CLASSIFY_EVENT_OUTCOMES` (excluding
`:companion_vetoed`, applied later in `run_mc`) to the event.

Decision chain (first match wins):

    vetoed_skin              → :skin_vetoed       (fast tracker veto)
    rejected_fv              → :outside_FV         (fast tracker reject)
    escaped                  → :escaped
    empty clusters           → :escaped
    !select_skin(stack,p)    → :skin_vetoed       (slow cumulative check)
    !select_FV(clusters,p)   → :outside_FV         (slow per-cluster check)
    >1 cluster               → :MS_rejected
    1 cluster:
        select_ROI(c,p) ? :SS_in_ROI : :SS_outside_ROI

Pure: no RNG (smearing already happened in `build_clusters`).
"""
function classify_event(status::Symbol,
                         stack::PhotonStack,
                         clusters::Vector{Cluster},
                         params::MCParams)::Symbol
    status === :vetoed_skin   && return :skin_vetoed
    status === :rejected_fv   && return :outside_FV
    status === :escaped       && return :escaped
    isempty(clusters)         && return :escaped
    select_skin(stack, params) || return :skin_vetoed
    select_FV(clusters, params) || return :outside_FV
    select_SC(clusters)        || return :MS_rejected
    return select_ROI(clusters[1], params) ? :SS_in_ROI : :SS_outside_ROI
end
