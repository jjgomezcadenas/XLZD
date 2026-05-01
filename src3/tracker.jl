# src3/tracker.jl — Stack-based per-photon tracker (NEW; replaces the
# legacy `track_one_photon!` in src3/mc.jl at cutover, step 11i).
#
# Status:
#   - 11a: scaffold + stub                                            DONE
#   - 11b: source sampling + transparent advance + first interaction  DONE
#          (forced PHOTO; one row per event)
#   - 11c: cross-section sampling (photo / Compton / pair)            TODO
#   - 11d: Compton recursion (`_track_child_photon!`)                 TODO
#   - 11e: pair production children (back-to-back 511 keV)            TODO
#   - 11f: skin + FV early-reject; rej_hist filling                   TODO
#   - 11g: integrate into run_mc behind kwarg                         TODO
#   - 11i: delete legacy tracker from src3/mc.jl                      TODO
#
# Region-class semantics:
#   :TPC, :Skin, :Inert          -> interactive LXe; sample with μ_LXe(e)
#   :FC, :Gas                    -> transparent (μ ≡ 0); advance through
#   :Outside                     -> photon left LXe; tracker returns :escaped
#
# At cutover (step 11i) the legacy tracker is deleted from src3/mc.jl
# and `path_to_next_region`, `companion_visible!`, `companion_reach_prob`
# remain there as the only contents.

using Random: AbstractRNG

"""
    track_photon_stack(rng, det, eff, xcom, params, stack;
                        rej_hist=nothing) -> status::Symbol

Stack-based per-photon tracker (work in progress; see file header for
the step-by-step status).

**Step 11b:** samples one γ from `eff`, walks through LXe regions
until either the photon exits (`:Outside`) or an interaction is
sampled. Interactions are currently FORCED to `INT_PHOTO` (full
deposit; cross-section sampling lands in 11c). The function pushes at
most one row to `stack` and returns `:completed` (deposit happened) or
`:escaped` (photon left LXe before interacting).

The returned symbol is a member of `TRACK_STATUSES` and is consumed by
`classify_event` to produce the per-event outcome.

`rej_hist` is accepted but unused until 11f.
"""
function track_photon_stack(rng::AbstractRNG, det::LXeDetector,
                             eff::EffectiveSource, xcom::XCOMTable,
                             params::MCParams, stack::PhotonStack;
                             rej_hist::Union{RejectionHistograms, Nothing}=nothing
                             )::Symbol
    # 1. Sample initial photon state from the source.
    x, y, z, dx, dy, dz = sample_entry(rng, det, eff)
    e             = eff.E_MeV
    parent_region = eff.region        # :barrel / :endcap_top / :endcap_bottom

    # 2. Transport loop: alternate between transparent advances (across
    #    :FC / :Gas regions) and interactions (μ-step in interactive LXe).
    ε = 1.0e-7   # boundary nudge; larger than path_to_next_region's EPS=1e-9
    while true
        region = region_at(det, x, y, z)

        # Photon left the LXe envelope -> :escaped.
        region === :Outside && return :escaped

        # Defensive below-threshold dump (never triggers in 11b — source γ
        # are all > 500 keV and there are no Compton recursions yet — but
        # the symmetric handling for child photons lands here in 11d).
        if e * 1000.0 < params.E_tracking_cutoff_keV
            push_row!(stack; nm=0, parent_region=parent_region,
                              region=region, interaction=INT_BELOW_THRESH,
                              x=x, y=y, z=z, epre=e, edep=e)
            return :completed
        end

        # :FC and :Gas are transparent — advance to the next boundary
        # without sampling an interaction. If no boundary lies ahead
        # (`Inf`), the photon escapes.
        if region === :FC || region === :Gas
            d = path_to_next_region(x, y, z, dx, dy, dz, det)
            d == Inf && return :escaped
            x += dx * (d + ε)
            y += dy * (d + ε)
            z += dz * (d + ε)
            continue
        end

        # Interactive LXe region: sample μ-step and compare with boundary.
        μ     = μ_LXe(det, e)
        dint  = -log(rand(rng)) / μ
        dnext = path_to_next_region(x, y, z, dx, dy, dz, det)

        if dnext < dint
            # Transparent advance to boundary; loop continues.
            x += dx * (dnext + ε)
            y += dy * (dnext + ε)
            z += dz * (dnext + ε)
        else
            # Interaction. Forced PHOTO until 11c.
            # Tag the row with `region` (the start region), not a
            # post-advance reclassification — dint was sampled from this
            # region's μ, so the interaction belongs to it. This also
            # avoids floating-point edge cases at region boundaries when
            # dint ≈ dnext.
            x += dx * dint
            y += dy * dint
            z += dz * dint
            push_row!(stack; nm=0, parent_region=parent_region,
                              region=region,
                              interaction=INT_PHOTO,
                              x=x, y=y, z=z, epre=e, edep=e)
            return :completed
        end
    end
end
