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

**Step 11d:** samples one γ from `eff` and delegates to
`_track_child_photon!` for transport + interactions. Compton-outgoing
γ are tracked through the cascade (multiple rows per event possible,
linked via the `nm` field). Pair-vertex rows are still terminal —
the two 511 keV annihilation γ are not yet spawned (lands in 11e).

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
            # Pair: vertex row only in 11d. The two 511 keV annihilation γ
            # are tracked as children in 11e.
            edep_pair = e - 2.0 * ME_C2_MEV
            push_row!(stack; nm=parent_ng, parent_region=parent_region,
                              region=region, interaction=INT_PAIR,
                              x=x, y=y, z=z, epre=e, edep=edep_pair)
            return
        end
    end
end
