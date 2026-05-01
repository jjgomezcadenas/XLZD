# src3/mc.jl — Per-photon Monte Carlo tracker.
#
# `track_one_photon!` shoots one γ from an EffectiveSource into the LXe
# detector and returns one of six outcome Symbols:
#
#   :escaped         — photon left without depositing visible energy
#   :MS_rejected     — multiple sites separated by ≥ Δz_threshold (or pair production)
#   :skin_vetoed     — accumulated skin deposit ≥ E_skin_veto
#   :SS_outside_FV   — single site, but cluster outside the FV box
#   :SS_outside_ROI  — single site in FV, energy outside the ROI window
#   :SS_in_ROI       — single site in FV, smeared energy in the ROI (the signal)
#
# The MC uses `region_at` to classify every interaction into one of
# {:active, :skin, :inert, :gas, :fc_region, :outside_lxe}. The two
# transparent regions (:gas, :fc_region) are advanced through without
# interaction; LXe regions interact with μ_LXe(E) and have region-specific
# deposit bookkeeping. Companion-γ veto for Tl-208 is NOT in this commit
# (handled in MC4).

# ---------------------------------------------------------------------------
# Boundary finder
# ---------------------------------------------------------------------------

"""
    path_to_next_region(x, y, z, dx, dy, dz, det) -> Float64

Smallest positive distance along the ray from (x, y, z) in direction
(dx, dy, dz) to ANY of the detector's region-defining surfaces:
  - cylindrical: r = R_FC_inner, R_FC_outer, R_ICV_inner
  - planar: z = z_gate, z_RFR_bottom, z_LXe_bottom, z_cathode

Returns `Inf` if the ray hits no boundary in finite distance.
"""
function path_to_next_region(x::Float64, y::Float64, z::Float64,
                              dx::Float64, dy::Float64, dz::Float64,
                              det::LXeDetector)::Float64
    EPS = 1.0e-9
    best = Inf

    # Cylindrical surfaces
    a = dx*dx + dy*dy
    b = 2.0 * (x*dx + y*dy)
    if a > EPS
        for R in (det.R_FC_inner, det.R_FC_outer, det.R_ICV_inner)
            c = x*x + y*y - R*R
            disc = b*b - 4.0*a*c
            disc < 0 && continue
            sq = sqrt(disc)
            for t in ((-b - sq)/(2.0*a), (-b + sq)/(2.0*a))
                if t > EPS && t < best
                    best = t
                end
            end
        end
    end

    # Axial planes
    if abs(dz) > EPS
        for z_plane in (det.z_gate, det.z_RFR_bottom,
                        det.z_LXe_bottom, det.z_cathode)
            t = (z_plane - z) / dz
            if t > EPS && t < best
                best = t
            end
        end
    end

    return best
end

# ---------------------------------------------------------------------------
# Photon state (mutable scratch struct kept short-lived)
# ---------------------------------------------------------------------------

mutable struct PhotonState
    cluster_started::Bool
    x_cluster::Float64
    y_cluster::Float64
    z_cluster::Float64
    E_cluster::Float64       # sum of deposits within Δz of first deposit
    E_skin_total::Float64
    nskin::Int               # count of :skin interactions
    skin_overflow::Bool      # E_skin_total ≥ veto threshold
    n_int::Int
    total_path::Float64
    outcome::Symbol
end

PhotonState() = PhotonState(false, NaN, NaN, NaN, 0.0, 0.0, 0, false, 0, 0.0, :in_progress)

# ---------------------------------------------------------------------------
# Per-region deposit handling
# ---------------------------------------------------------------------------

"""
    handle_deposit!(state, det, params, x, y, z, E_dep_MeV)

Apply the deposit bookkeeping for the given interaction point. May set
`state.outcome` to `:MS_rejected` (Δz exceeded) or `:skin_vetoed`
(skin total ≥ veto). Does nothing for `:inert`, `:gas`, `:fc_region`,
or `:outside_lxe`.
"""
function handle_deposit!(state::PhotonState, det::LXeDetector,
                          params::MCParams,
                          x::Float64, y::Float64, z::Float64,
                          E_dep::Float64,
                          scratch::Union{PhotonScratch, Nothing}=nothing)
    reg = region_at(det, x, y, z)
    E_visible_MeV  = det.E_visible_keV   * 1.0e-3
    E_skin_veto_MeV = det.E_skin_veto_keV * 1.0e-3
    Δz_thresh = Δz_threshold_cm(params)

    if reg === :active
        E_dep < E_visible_MeV && return
        # Record for control histograms / cluster computation
        scratch !== nothing && push!(scratch.deposits,
                                      LXeDeposit(x, y, z, E_dep, :active))
        if !state.cluster_started
            state.x_cluster = x
            state.y_cluster = y
            state.z_cluster = z
            state.E_cluster = E_dep
            state.cluster_started = true
        else
            # Sum deposits within Δz of first; do NOT MS-reject during tracking.
            # Multi-cluster classification is done post-tracking.
            if abs(z - state.z_cluster) < Δz_thresh
                state.E_cluster += E_dep
            end
        end
    elseif reg === :skin
        scratch !== nothing && push!(scratch.deposits,
                                      LXeDeposit(x, y, z, E_dep, :skin))
        state.E_skin_total += E_dep
        state.nskin        += 1
        if state.E_skin_total >= E_skin_veto_MeV
            state.skin_overflow = true
        end
    end
    # :inert / :gas / :fc_region / :outside_lxe → no signal
    nothing
end

# ---------------------------------------------------------------------------
# Finalisation: classify cluster -> SS_in_ROI / SS_outside_ROI / SS_outside_FV / escaped
# ---------------------------------------------------------------------------

"""
    finalize_outcome!(state, rng, params, scratch=nothing) -> Symbol

Post-tracking classification. Called once after `_track_photon_segment!`
returns. Sets `state.outcome` based on:
  - `state.skin_overflow` → :skin_vetoed
  - no cluster started   → :escaped
  - first :active deposit outside FV → :SS_outside_FV
  - else: cluster grouping from `scratch.deposits`. ng > 1 → :MS_rejected;
    ng = 1 → smear E and apply ROI cut.

The cluster-grouping pass uses `compute_clusters` (defined in
`histograms.jl`) so SC/MS classification is consistent with the control
histograms.
"""
function finalize_outcome!(state::PhotonState, rng::AbstractRNG,
                            params::MCParams,
                            scratch::Union{PhotonScratch, Nothing}=nothing)
    state.outcome === :in_progress || return

    # Skin veto if not already early-rejected
    if state.skin_overflow
        state.outcome = :skin_vetoed
        return
    end

    if !state.cluster_started
        state.outcome = :escaped
        return
    end

    # First :active deposit outside FV?
    if !in_fv(state.x_cluster, state.y_cluster, state.z_cluster, params)
        state.outcome = :SS_outside_FV
        return
    end

    # SS / MS classification from cluster grouping
    clusters = scratch === nothing ?
               Cluster[Cluster(state.x_cluster, state.y_cluster,
                               state.z_cluster, state.E_cluster)] :
               compute_clusters(scratch.deposits, params)
    if length(clusters) > 1
        state.outcome = :MS_rejected
        return
    end
    # SS in FV — smear and apply ROI cut
    ec = clusters[1].ec
    state.outcome = classify_ss_energy(rng, ec, params)
    nothing
end

# ---------------------------------------------------------------------------
# Photon tracking segment (recursive helper)
# ---------------------------------------------------------------------------

"""
    _track_photon_segment!(rng, det, params, xcom, state, scratch,
                            x, y, z, dx, dy, dz, E)

Track one photon segment through the LXe until it is absorbed, escapes
the detector, or `state.outcome` becomes terminal. Mutates `state` and
(if non-`nothing`) pushes deposits to `scratch`. Returns nothing.

Pair production triggers two recursive calls — one per 511 keV
annihilation γ emitted back-to-back from the conversion vertex.
Recursion depth is bounded at 2: 511 keV is below the pair-production
threshold, so secondaries cannot themselves pair-produce.
"""
function _track_photon_segment!(rng::AbstractRNG, det::LXeDetector,
                                 params::MCParams, xcom::XCOMTable,
                                 state::PhotonState,
                                 scratch::Union{PhotonScratch, Nothing},
                                 rej_hist::Union{RejectionHistograms, Nothing},
                                 early_skin_reject::Bool,
                                 early_fv_reject::Bool,
                                 x::Float64, y::Float64, z::Float64,
                                 dx::Float64, dy::Float64, dz::Float64,
                                 E::Float64)
    ρ = det.material.density
    E_cutoff = E_tracking_cutoff_MeV(params)

    while state.outcome === :in_progress
        reg = region_at(det, x, y, z)

        if reg === :outside_lxe
            return
        end

        if reg === :gas || reg === :fc_region
            d = path_to_next_region(x, y, z, dx, dy, dz, det)
            d == Inf && return
            d_step = d + 1.0e-7
            x += d_step * dx; y += d_step * dy; z += d_step * dz
            state.total_path += d_step
            continue
        end

        # LXe region: μ-step vs region boundary
        μ = μ_total_lin(xcom, E, ρ)
        d_int = -log(rand(rng)) / μ
        d_bnd = path_to_next_region(x, y, z, dx, dy, dz, det)

        if d_int >= d_bnd
            d_step = d_bnd + 1.0e-7
            x += d_step * dx; y += d_step * dy; z += d_step * dz
            state.total_path += d_step
            continue
        end

        # Interact
        x += d_int * dx; y += d_int * dy; z += d_int * dz
        state.total_path += d_int
        state.n_int += 1

        sP = σ_photo(xcom, E)
        sC = σ_Compton(xcom, E)
        sX = σ_pair(xcom, E)
        sT = sP + sC + sX
        r_branch = rand(rng)

        # Helper: post-deposit checks for early rejection (skin / FV).
        # Captures the current (x, y, z, E_dep_at_this_call) at each call
        # site by closure over `was_cluster_started_before` and the
        # most-recent E_dep variable.
        function _post_deposit_checks!(E_dep_here::Float64,
                                       prev_cluster_started::Bool)
            # Skin overflow → write rejection histogram + early-reject if requested
            if state.skin_overflow && state.outcome === :in_progress
                if early_skin_reject
                    if rej_hist !== nothing
                        fill_rejected_skin!(rej_hist, x, y, z, E_dep_here)
                    end
                    state.outcome = :skin_vetoed
                end
            end
            # First :active deposit just happened? Check FV
            if state.outcome === :in_progress &&
               !prev_cluster_started && state.cluster_started
                if !in_fv(state.x_cluster, state.y_cluster, state.z_cluster, params)
                    if early_fv_reject
                        if rej_hist !== nothing
                            fill_rejected_fv!(rej_hist,
                                              state.x_cluster, state.y_cluster,
                                              state.z_cluster, state.E_cluster)
                        end
                        state.outcome = :SS_outside_FV
                    end
                end
            end
        end

        if r_branch < sX / sT
            # Pair production: e⁺e⁻ kinetic energy deposited at the vertex,
            # then two 511 keV annihilation γ's back-to-back.
            E_pair_kin = E - 2.0 * ME_C2_MEV
            prev_cs = state.cluster_started
            handle_deposit!(state, det, params, x, y, z, E_pair_kin, scratch)
            _post_deposit_checks!(E_pair_kin, prev_cs)
            state.outcome === :in_progress || return

            # Isotropic direction for the first; second is opposite.
            u_aux = 2.0 * rand(rng) - 1.0
            φ_aux = 2π * rand(rng)
            s_aux = sqrt(max(0.0, 1.0 - u_aux*u_aux))
            d1x = s_aux * cos(φ_aux)
            d1y = s_aux * sin(φ_aux)
            d1z = u_aux

            _track_photon_segment!(rng, det, params, xcom, state, scratch,
                                    rej_hist, early_skin_reject, early_fv_reject,
                                    x, y, z,  d1x,  d1y,  d1z, ME_C2_MEV)
            state.outcome === :in_progress || return
            _track_photon_segment!(rng, det, params, xcom, state, scratch,
                                    rej_hist, early_skin_reject, early_fv_reject,
                                    x, y, z, -d1x, -d1y, -d1z, ME_C2_MEV)
            return

        elseif r_branch < (sX + sP) / sT
            # Photoelectric absorption: full energy locally
            prev_cs = state.cluster_started
            handle_deposit!(state, det, params, x, y, z, E, scratch)
            _post_deposit_checks!(E, prev_cs)
            return

        else
            # Compton: KN scatter
            E_scatt, cos_θ = sample_klein_nishina(rng, E)
            E_dep = E - E_scatt
            prev_cs = state.cluster_started
            handle_deposit!(state, det, params, x, y, z, E_dep, scratch)
            _post_deposit_checks!(E_dep, prev_cs)
            state.outcome === :in_progress || return

            φ_az = 2π * rand(rng)
            dx, dy, dz = rotate_direction(dx, dy, dz, cos_θ, φ_az)
            E = E_scatt

            # Energy cutoff: deposit remainder locally
            if E < E_cutoff
                prev_cs = state.cluster_started
                handle_deposit!(state, det, params, x, y, z, E, scratch)
                _post_deposit_checks!(E, prev_cs)
                return
            end
        end
    end
    return
end

# ---------------------------------------------------------------------------
# Main tracker
# ---------------------------------------------------------------------------

"""
    track_one_photon!(rng, det, eff, xcom, params) -> Symbol

Shoot one γ from `eff` into `det`. Pair-production interactions are
fully tracked: the e⁺e⁻ kinetic energy is deposited at the conversion
vertex and the two 511 keV annihilation γ's are propagated through the
LXe with full physics.

Returns one of `:escaped`, `:MS_rejected`, `:skin_vetoed`,
`:SS_outside_FV`, `:SS_outside_ROI`, `:SS_in_ROI`.

Companion-γ veto is applied separately in `run_mc`.
"""
function track_one_photon!(rng::AbstractRNG,
                            det::LXeDetector,
                            eff::EffectiveSource,
                            xcom::XCOMTable,
                            params::MCParams;
                            hist::Union{HistogramSet, Nothing}=nothing,
                            scratch::Union{PhotonScratch, Nothing}=nothing,
                            rej_hist::Union{RejectionHistograms, Nothing}=nothing,
                            early_skin_reject::Bool=true,
                            early_fv_reject::Bool=true)::Symbol
    scratch !== nothing && empty!(scratch.deposits)

    x, y, z, dx, dy, dz = sample_entry(rng, det, eff)
    state = PhotonState()
    _track_photon_segment!(rng, det, params, xcom, state, scratch,
                            rej_hist, early_skin_reject, early_fv_reject,
                            x, y, z, dx, dy, dz, eff.E_MeV)

    # Post-tracking classification (handles MS/SS, skin overflow, FV cut).
    if state.outcome === :in_progress
        finalize_outcome!(state, rng, params, scratch)
    end

    if hist !== nothing && scratch !== nothing
        update_histograms!(hist, scratch, params)
    end

    return state.outcome
end

# ---------------------------------------------------------------------------
# Companion veto (Tl-208)
# ---------------------------------------------------------------------------

"""
    companion_visible!(rng, det, comp_eff, xcom, params) -> Bool

Track one companion γ from `comp_eff` (lumped at 583 keV) into the
LXe and return `true` iff it produces a visible deposit anywhere —
i.e. any single deposit ≥ `det.E_visible_keV` in `:active` LXe, OR the
accumulated skin deposit ≥ `det.E_skin_veto_keV`. Pair production is
treated as visible (the two annihilation γ guarantee a visible deposit).

This is a stripped-down version of `track_one_photon!`: no SS/MS
clustering, no FV cut, no ROI smearing — just a Bool answer for the
veto step. Intended to be called only when the main γ has already
been tagged `:SS_in_ROI`.
"""
function companion_visible!(rng::AbstractRNG, det::LXeDetector,
                             comp_eff::EffectiveSource,
                             xcom::XCOMTable,
                             params::MCParams)::Bool
    x, y, z, dx, dy, dz = sample_entry(rng, det, comp_eff)
    E = comp_eff.E_MeV
    E_visible_MeV   = det.E_visible_keV   * 1.0e-3
    E_skin_veto_MeV = det.E_skin_veto_keV * 1.0e-3
    E_cutoff = E_tracking_cutoff_MeV(params)
    ρ = det.material.density
    E_skin_total = 0.0

    while true
        reg = region_at(det, x, y, z)
        if reg === :outside_lxe
            return false
        end

        if reg === :gas || reg === :fc_region
            d = path_to_next_region(x, y, z, dx, dy, dz, det)
            d == Inf && return false
            d_step = d + 1.0e-7
            x += d_step * dx; y += d_step * dy; z += d_step * dz
            continue
        end

        # LXe region: sample step length and compare to next boundary
        μ = μ_total_lin(xcom, E, ρ)
        d_int = -log(rand(rng)) / μ
        d_bnd = path_to_next_region(x, y, z, dx, dy, dz, det)
        if d_int >= d_bnd
            d_step = d_bnd + 1.0e-7
            x += d_step * dx; y += d_step * dy; z += d_step * dz
            continue
        end

        # Interact
        x += d_int * dx; y += d_int * dy; z += d_int * dz

        sP = σ_photo(xcom, E)
        sC = σ_Compton(xcom, E)
        sX = σ_pair(xcom, E)
        sT = sP + sC + sX
        r_branch = rand(rng)

        # Determine deposit + termination
        E_dep = 0.0
        terminate = false
        if r_branch < sX / sT
            # Pair production: 511 keV γ pair certainly deposit visibly nearby
            return true
        elseif r_branch < (sX + sP) / sT
            E_dep = E
            terminate = true
        else
            # Compton
            E_scatt, cos_θ = sample_klein_nishina(rng, E)
            E_dep = E - E_scatt
            φ_az = 2π * rand(rng)
            dx, dy, dz = rotate_direction(dx, dy, dz, cos_θ, φ_az)
            E = E_scatt
            if E < E_cutoff
                E_dep += E
                terminate = true
            end
        end

        # Visibility check on this deposit
        reg_dep = region_at(det, x, y, z)
        if reg_dep === :active
            E_dep >= E_visible_MeV && return true
        elseif reg_dep === :skin
            E_skin_total += E_dep
            E_skin_total >= E_skin_veto_MeV && return true
        end

        terminate && return false
    end
end

# ---------------------------------------------------------------------------
# Tl-208 helper: per-decay probability that the companion reaches the LXe
# ---------------------------------------------------------------------------

"""
    companion_reach_prob(comp_eff::EffectiveSource) -> Float64

Probability that, given a Tl-208 decay in the contributing physical
volumes, the cascade companion γ reaches the LXe inner exit surface.
Equal to `comp_eff.total_per_yr / decay_rate`, where
`decay_rate = (Σ contrib.source.produced_per_yr) / BR_TL208_COMPANION`.
"""
function companion_reach_prob(comp_eff::EffectiveSource)::Float64
    total_produced = sum(c.source.produced_per_yr for c in comp_eff.contributions)
    return (comp_eff.total_per_yr / total_produced) * BR_TL208_COMPANION
end

