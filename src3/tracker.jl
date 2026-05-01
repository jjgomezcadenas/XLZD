# src3/tracker.jl — Stack-based per-photon tracker (NEW; replaces the
# legacy `track_one_photon!` in src3/mc.jl at cutover, step 11i).
#
# Contract (final, after step 11f):
#
#   track_photon_stack(rng, det, eff, xcom, params, stack;
#                       rej_hist=nothing) -> status::Symbol
#
# 1. Sample one γ from `eff` (entry point + direction at the ICV inner
#    surface). Mutate `stack` in place; the caller is expected to have
#    called `empty!(stack)` first if reusing.
#
# 2. Transport the photon through the LXe regions (TPC / Skin / Inert /
#    Gas) using μ_LXe(epre) for the interaction-distance sampling and
#    `path_to_next_region` for the region-boundary distance. Whichever
#    is smaller wins.
#
# 3. At each interaction, sample type via cross-section ratios
#    (photo / Compton / pair). Push one StackRow per interaction. For
#    Compton, recurse on the outgoing γ. For pair, push the vertex row,
#    spawn two 511 keV children with a back-to-back random axis, and
#    recurse on each.
#
# 4. Below-threshold termination: when a photon's energy falls below
#    `params.E_tracking_cutoff_keV` between interactions, push one
#    INT_BELOW_THRESH row at the current position with edep equal to
#    the residual energy, and return.
#
# 5. Early-reject paths (write to `rej_hist` if provided, then return):
#       - first :Skin interaction with cumulative skin energy above
#         `det.E_skin_veto_keV`  -> :vetoed_skin
#       - first :TPC interaction outside the FV box                -> :rejected_fv
#
# 6. Status return values (subset of TRACK_STATUSES):
#       :completed   — tracker ran to completion; stack populated
#       :escaped     — photon left LXe with no visible deposit
#       :vetoed_skin — early-reject due to skin overflow
#       :rejected_fv — early-reject due to FV box check
#
# Helpers `_propagate_to_first_lxe`, `_first_interaction`,
# `_track_child_photon!` will be added in steps 11b–11e as their
# callers materialize. They are deliberately NOT scaffolded here to
# avoid unused-stub noise.
#
# At cutover (step 11i) the legacy tracker is deleted from src3/mc.jl
# and `path_to_next_region`, `companion_visible!`, `companion_reach_prob`
# remain there as the only contents.

using Random: AbstractRNG

"""
    track_photon_stack(rng, det, eff, xcom, params, stack;
                        rej_hist=nothing) -> status::Symbol

Stack-based per-photon tracker (work in progress; see file header for
the final contract).

**Step 11a stub:** returns `:escaped` and leaves `stack` untouched.
The full transport, interaction, recursion, and early-reject logic
land in steps 11b through 11f. Until then `track_photon_stack` is not
called from `run_mc` — the legacy `track_one_photon!` (in `src3/mc.jl`)
remains the active tracker.

The returned symbol is a member of `TRACK_STATUSES` and is consumed by
`classify_event` to produce the per-event outcome.
"""
function track_photon_stack(rng::AbstractRNG, det::LXeDetector,
                             eff::EffectiveSource, xcom::XCOMTable,
                             params::MCParams, stack::PhotonStack;
                             rej_hist::Union{RejectionHistograms, Nothing}=nothing
                             )::Symbol
    # 11a stub: full implementation lands across steps 11b–11f.
    return :escaped
end
