# src3/tracker.jl — Stack-based per-photon tracker (NEW; replaces the
# legacy `track_one_photon!` in src3/mc.jl at cutover, step 11i).
#
# Status:
#   - 11a: scaffold + stub                                            DONE
#   - 11b: source sampling + transparent advance + first interaction  DONE
#          (forced PHOTO; one row per event)
#   - 11c: cross-section sampling (photo / Compton / pair)            DONE
#          (one row per event; no recursion or children)
#   - 11d: Compton recursion (`_track_child_photon!`)                 DONE
#          (Compton outgoing γ tracked; pair children still TODO)
#   - 11e: pair production children (back-to-back 511 keV)            DONE
#   - 11f: fast_veto + slow checks in classify_event                  DONE
#          (fast_veto here; select_skin / select_FV in select.jl;
#           classify_event composes them; not yet called from run_mc)
#   - 11g: integrate fast_veto + track_photon_stack into run_mc       TODO
#          via per-thread RNG snapshot/restore
#   - 11i: delete legacy tracker from src3/mc.jl                      TODO
#
# Region-class semantics:
#   :TPC, :Skin, :Inert          -> interactive LXe; sample with μ_LXe(e)
#   :FC, :Gas                    -> transparent (μ ≡ 0); advance through
#   :Outside                     -> photon left LXe; tracker returns :escaped
#
# Mother-chain semantics (parent_region, parent_ng):
#   - For the source γ's first interaction: parent_region = eff.region
#     (e.g. :barrel / :endcap_top / :endcap_bottom), parent_ng = 0.
#   - For a child γ (Compton outgoing): parent_region = region of the
#     creating Compton vertex, parent_ng = ng of that vertex row.
#   - Updated tail-loop-style inside `_track_child_photon!` so a single
#     call frame handles arbitrarily long Compton chains (no actual
#     recursion -> no risk of stack overflow).
#
# At cutover (step 11i) the legacy tracker is deleted from src3/mc.jl
# and `path_to_next_region`, `companion_visible!`, `companion_reach_prob`
# remain there as the only contents.

using Random: AbstractRNG

# Boundary nudge for transparent advances; larger than path_to_next_region's
# internal EPS=1e-9 to ensure we cross the boundary on every step.
const _TRACKER_BOUNDARY_NUDGE = 1.0e-7

"""
    track_photon_stack(rng, det, eff, xcom, params, stack;
                        rej_hist=nothing) -> status::Symbol

Stack-based per-photon tracker (work in progress; see file header for
the step-by-step status).

**Step 11e:** samples one γ from `eff` and delegates to
`_track_child_photon!` for transport + interactions. Compton-outgoing
γ are tracked through the cascade (tail-loop, no real recursion). At
a pair vertex, two 511 keV annihilation γ are spawned at the vertex
position with back-to-back random axis; each is tracked by its own
`_track_child_photon!` call (one level of real recursion). Pair γ's
cannot pair-produce again (511 keV << 1.022 MeV threshold), so any
mother-chain walk contains at most one INT_PAIR row.

The status return is derived from whether any row was pushed:
  `length(stack) == 0` -> :escaped (source photon left LXe before any
                                     interaction).
  `length(stack) >= 1` -> :completed.

`rej_hist` is accepted but unused until 11f.
"""
function track_photon_stack(rng::AbstractRNG, det::LXeDetector,
                             eff::EffectiveSource, xcom::XCOMTable,
                             params::MCParams, stack::PhotonStack;
                             rej_hist::Union{RejectionHistograms, Nothing}=nothing
                             )::Symbol
    x, y, z, dx, dy, dz = sample_entry(rng, det, eff)
    _track_child_photon!(stack, rng, det, xcom, params,
                          x, y, z, dx, dy, dz, eff.E_MeV,
                          eff.region, 0)
    return length(stack) > 0 ? :completed : :escaped
end

# ---------------------------------------------------------------------------
# Internal: transport + interaction loop for one photon (and, via tail-loop,
# its Compton-outgoing γ's). Mutates `stack` in place.
# ---------------------------------------------------------------------------
"""
    _track_child_photon!(stack, rng, det, xcom, params,
                          x, y, z, dx, dy, dz, e,
                          parent_region, parent_ng)

Transport one photon from `(x, y, z)` along `(dx, dy, dz)` with energy
`e` (MeV) until it terminates (PHOTO absorption, PAIR vertex, escape,
or below-threshold). On Compton, tail-loop into the outgoing γ with
updated `(parent_region, parent_ng)` linking back to the Compton vertex.
"""
function _track_child_photon!(stack::PhotonStack,
                               rng::AbstractRNG,
                               det::LXeDetector,
                               xcom::XCOMTable,
                               params::MCParams,
                               x::Float64, y::Float64, z::Float64,
                               dx::Float64, dy::Float64, dz::Float64,
                               e::Float64,
                               parent_region::Symbol,
                               parent_ng::Int)
    ε = _TRACKER_BOUNDARY_NUDGE
    while true
        region = region_at(det, x, y, z)

        # Photon left LXe envelope -> done.
        region === :Outside && return

        # Below-threshold: dump residual energy locally, terminate.
        if e * 1000.0 < params.E_tracking_cutoff_keV
            push_row!(stack; nm=parent_ng, parent_region=parent_region,
                              region=region, interaction=INT_BELOW_THRESH,
                              x=x, y=y, z=z, epre=e, edep=e)
            return
        end

        # :FC and :Gas are transparent — advance to next boundary or escape.
        if region === :FC || region === :Gas
            d = path_to_next_region(x, y, z, dx, dy, dz, det)
            d == Inf && return
            x += dx * (d + ε)
            y += dy * (d + ε)
            z += dz * (d + ε)
            continue
        end

        # Interactive LXe region: μ-step vs region boundary.
        μ     = μ_LXe(det, e)
        dint  = -log(rand(rng)) / μ
        dnext = path_to_next_region(x, y, z, dx, dy, dz, det)

        if dnext < dint
            # Transparent advance to boundary; loop continues.
            x += dx * (dnext + ε)
            y += dy * (dnext + ε)
            z += dz * (dnext + ε)
            continue
        end

        # Interaction at distance dint inside `region`. Row tagged with
        # `region` (start region) — dint was sampled from its μ.
        x += dx * dint
        y += dy * dint
        z += dz * dint

        sP = σ_photo(xcom, e)
        sC = σ_Compton(xcom, e)
        sX = σ_pair(xcom, e)        # 0 below 1.022 MeV
        sT = sP + sC + sX
        r  = rand(rng) * sT

        if r < sP
            # Photoelectric: full absorption — terminal.
            push_row!(stack; nm=parent_ng, parent_region=parent_region,
                              region=region, interaction=INT_PHOTO,
                              x=x, y=y, z=z, epre=e, edep=e)
            return
        elseif r < sP + sC
            # Compton: push vertex row, then tail-loop on outgoing γ.
            E_prime, cos_θ = sample_klein_nishina(rng, e)
            edep_e = e - E_prime
            ng_self = push_row!(stack; nm=parent_ng,
                                        parent_region=parent_region,
                                        region=region,
                                        interaction=INT_COMPTON,
                                        x=x, y=y, z=z,
                                        epre=e, edep=edep_e)
            # Sample azimuth around the input direction (KN is φ-symmetric)
            # and rotate to obtain the scattered γ's new direction.
            φ = 2.0π * rand(rng)
            ndx, ndy, ndz = rotate_direction(dx, dy, dz, cos_θ, φ)

            # Tail-loop: track the outgoing γ from the same (x,y,z),
            # with new direction, energy, and parent linkage.
            dx, dy, dz    = ndx, ndy, ndz
            e             = E_prime
            parent_region = region
            parent_ng     = ng_self
            continue
        else
            # Pair: push the vertex row, then spawn two 511 keV
            # annihilation γ at the vertex with back-to-back direction.
            # Each child is tracked by its own _track_child_photon! call.
            edep_pair = e - 2.0 * ME_C2_MEV
            ng_self   = push_row!(stack; nm=parent_ng,
                                          parent_region=parent_region,
                                          region=region,
                                          interaction=INT_PAIR,
                                          x=x, y=y, z=z,
                                          epre=e, edep=edep_pair)

            # Sample one isotropic direction in the lab frame; the second
            # γ is the antipode. Standard approximation — positron is
            # assumed to thermalize in negligible distance/time, so the
            # CM frame is approximately at rest in the lab.
            cos_θ_iso = 2.0 * rand(rng) - 1.0
            sin_θ_iso = sqrt(max(0.0, 1.0 - cos_θ_iso * cos_θ_iso))
            φ_iso     = 2.0π * rand(rng)
            g1x = sin_θ_iso * cos(φ_iso)
            g1y = sin_θ_iso * sin(φ_iso)
            g1z = cos_θ_iso

            # Track first 511 keV γ.
            _track_child_photon!(stack, rng, det, xcom, params,
                                  x, y, z,  g1x,  g1y,  g1z,
                                  ME_C2_MEV, region, ng_self)
            # Track second 511 keV γ (antipodal direction).
            _track_child_photon!(stack, rng, det, xcom, params,
                                  x, y, z, -g1x, -g1y, -g1z,
                                  ME_C2_MEV, region, ng_self)
            return
        end
    end
end

# ===========================================================================
# fast_veto — pre-tracking single-deposit reject screen
# ===========================================================================
"""
    fast_veto(rng, det, eff, xcom, params; rej_hist=nothing) -> Symbol

Quick first-interaction screen returning one of:
  :pass         — first interaction does NOT trigger reject; caller
                  should restore the RNG snapshot and run
                  `track_photon_stack` for full tracking.
  :vetoed_skin  — first interaction in `:Skin` with edep > E_skin_veto.
  :rejected_fv  — first interaction in `:TPC` outside FV with
                  edep > E_visible.

Self-contained: samples source γ, transports through transparent
regions, samples the FIRST interaction (μ-step + interaction-type)
and applies the threshold check. Does NOT push to any stack and does
NOT call `track_photon_stack` or `_track_child_photon!`.

Caller (`run_mc` in 11g) should:
    copy!(rng_snapshot, rng)
    fv = fast_veto(rng, det, eff, xcom, params; rej_hist)
    if fv === :pass
        copy!(rng, rng_snapshot)
        track_photon_stack(rng, det, eff, xcom, params, stack)
        # ...build_clusters / classify_event...
    end
so that `track_photon_stack` re-samples the SAME first interaction
(deterministic given equal RNG state) and continues the cascade.

Threshold sources: `params.E_visible_keV`, `params.E_skin_veto_keV`.
"""
function fast_veto(rng::AbstractRNG, det::LXeDetector,
                    eff::EffectiveSource, xcom::XCOMTable,
                    params::MCParams;
                    rej_hist::Union{RejectionHistograms, Nothing}=nothing
                    )::Symbol
    x, y, z, dx, dy, dz = sample_entry(rng, det, eff)
    e   = eff.E_MeV
    ε   = _TRACKER_BOUNDARY_NUDGE
    E_visible_MeV   = params.E_visible_keV   / 1000.0
    E_skin_veto_MeV = params.E_skin_veto_keV / 1000.0

    while true
        region = region_at(det, x, y, z)

        # Photon left LXe envelope before any interaction → no reject possible.
        region === :Outside && return :pass

        # Sub-cutoff photon (defensive; source γ are all > 500 keV).
        e * 1000.0 < params.E_tracking_cutoff_keV && return :pass

        # :FC and :Gas are transparent — advance to next boundary or escape.
        if region === :FC || region === :Gas
            d = path_to_next_region(x, y, z, dx, dy, dz, det)
            d == Inf && return :pass
            x += dx * (d + ε); y += dy * (d + ε); z += dz * (d + ε)
            continue
        end

        # Interactive LXe region: μ-step vs region boundary.
        μ     = μ_LXe(det, e)
        dint  = -log(rand(rng)) / μ
        dnext = path_to_next_region(x, y, z, dx, dy, dz, det)

        if dnext < dint
            x += dx * (dnext + ε); y += dy * (dnext + ε); z += dz * (dnext + ε)
            continue
        end

        # First interaction inside `region` at distance dint.
        x += dx * dint
        y += dy * dint
        z += dz * dint

        # Sample interaction type and compute edep.
        sP = σ_photo(xcom, e)
        sC = σ_Compton(xcom, e)
        sX = σ_pair(xcom, e)
        sT = sP + sC + sX
        r  = rand(rng) * sT

        edep = if r < sP
            e
        elseif r < sP + sC
            E_prime, _cosθ = sample_klein_nishina(rng, e)
            e - E_prime
        else
            e - 2.0 * ME_C2_MEV
        end

        # Apply the two fast-reject thresholds.
        if region === :TPC && !in_fv(x, y, z, params) && edep > E_visible_MeV
            rej_hist !== nothing && fill_rejected_fv!(rej_hist, x, y, z, edep)
            return :rejected_fv
        end
        if region === :Skin && edep > E_skin_veto_MeV
            rej_hist !== nothing && fill_rejected_skin!(rej_hist, x, y, z, edep)
            return :vetoed_skin
        end

        # First interaction did not trigger fast reject; defer to slow checks.
        return :pass
    end
end
