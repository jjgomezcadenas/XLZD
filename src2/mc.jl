# src2/mc.jl — Per-photon Monte Carlo tracker.
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
    E_cluster::Float64
    E_skin_total::Float64
    n_int::Int
    total_path::Float64
    outcome::Symbol
end

PhotonState() = PhotonState(false, NaN, NaN, NaN, 0.0, 0.0, 0, 0.0, :in_progress)

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
                          E_dep::Float64)
    reg = region_at(det, x, y, z)
    E_visible_MeV  = det.E_visible_keV   * 1.0e-3
    E_skin_veto_MeV = det.E_skin_veto_keV * 1.0e-3
    Δz_thresh = Δz_threshold_cm(params)

    if reg === :active
        E_dep < E_visible_MeV && return
        if !state.cluster_started
            state.x_cluster = x
            state.y_cluster = y
            state.z_cluster = z
            state.E_cluster = E_dep
            state.cluster_started = true
        else
            if abs(z - state.z_cluster) < Δz_thresh
                state.E_cluster += E_dep
            else
                state.outcome = :MS_rejected
            end
        end
    elseif reg === :skin
        state.E_skin_total += E_dep
        if state.E_skin_total >= E_skin_veto_MeV
            state.outcome = :skin_vetoed
        end
    end
    # :inert / :gas / :fc_region / :outside_lxe → no signal
    nothing
end

# ---------------------------------------------------------------------------
# Finalisation: classify cluster -> SS_in_ROI / SS_outside_ROI / SS_outside_FV / escaped
# ---------------------------------------------------------------------------

function finalize_outcome!(state::PhotonState, rng::AbstractRNG,
                            params::MCParams)
    state.outcome === :in_progress || return
    if state.cluster_started
        if in_fv(state.x_cluster, state.y_cluster, state.z_cluster, params)
            state.outcome = classify_ss_energy(rng, state.E_cluster, params)
        else
            state.outcome = :SS_outside_FV
        end
    else
        state.outcome = :escaped
    end
    nothing
end

# ---------------------------------------------------------------------------
# Main tracker
# ---------------------------------------------------------------------------

"""
    track_one_photon!(rng, det, eff, xcom, params) -> Symbol

Shoot one γ from `eff` into `det`, sampling step lengths via μ_LXe and
sampling interaction types from `xcom`'s per-component cross sections.
Returns one of `:escaped`, `:MS_rejected`, `:skin_vetoed`,
`:SS_outside_FV`, `:SS_outside_ROI`, `:SS_in_ROI`.

Companion-γ veto is NOT applied here — it will be added in MC4.
"""
function track_one_photon!(rng::AbstractRNG,
                            det::LXeDetector,
                            eff::EffectiveSource,
                            xcom::XCOMTable,
                            params::MCParams)::Symbol
    # Sample entry point + initial direction
    x, y, z, dx, dy, dz = sample_entry(rng, det, eff)
    E = eff.E_MeV
    state = PhotonState()
    ρ = det.material.density
    E_cutoff = E_tracking_cutoff_MeV(params)

    # Tracking loop. The photon may pass through transparent regions
    # (:gas, :fc_region) without interaction; in LXe regions it samples
    # an exponential step length and either interacts or hits a region
    # boundary first (in which case we re-classify).
    while state.outcome === :in_progress
        reg = region_at(det, x, y, z)

        if reg === :outside_lxe
            # Photon left the LXe envelope (entered the Ti vessel or beyond).
            finalize_outcome!(state, rng, params)
            break
        end

        if reg === :gas || reg === :fc_region
            d = path_to_next_region(x, y, z, dx, dy, dz, det)
            if d == Inf
                # No further boundary — photon escapes
                finalize_outcome!(state, rng, params)
                break
            end
            # Step slightly past the boundary so re-classification lands
            # in the next region rather than on the surface.
            d_step = d + 1.0e-7
            x += d_step * dx
            y += d_step * dy
            z += d_step * dz
            state.total_path += d_step
            continue
        end

        # reg ∈ (:active, :skin, :inert) — LXe with full attenuation
        μ = μ_total_lin(xcom, E, ρ)
        d_int = -log(rand(rng)) / μ
        d_bnd = path_to_next_region(x, y, z, dx, dy, dz, det)

        if d_int >= d_bnd
            # Cross region boundary before interacting
            d_step = d_bnd + 1.0e-7
            x += d_step * dx
            y += d_step * dy
            z += d_step * dz
            state.total_path += d_step
            continue
        end

        # Interact
        x += d_int * dx
        y += d_int * dy
        z += d_int * dz
        state.total_path += d_int
        state.n_int += 1

        # Sample interaction type by σ-fractions
        sP = σ_photo(xcom, E)
        sC = σ_Compton(xcom, E)
        sX = σ_pair(xcom, E)
        sT = sP + sC + sX
        r_branch = rand(rng)

        if r_branch < sX / sT
            # Pair production → MS rejected (the two 511 keV annihilation
            # γ travel ~6 cm in LXe, virtually always producing spatially
            # separated deposits).
            state.outcome = :MS_rejected
            break

        elseif r_branch < (sX + sP) / sT
            # Photoelectric absorption: full energy locally
            handle_deposit!(state, det, params, x, y, z, E)
            if state.outcome === :in_progress
                finalize_outcome!(state, rng, params)
            end
            break

        else
            # Compton: KN scatter
            E_scatt, cos_θ = sample_klein_nishina(rng, E)
            E_dep = E - E_scatt
            handle_deposit!(state, det, params, x, y, z, E_dep)
            state.outcome === :in_progress || break

            # Update direction & energy
            φ_az = 2π * rand(rng)
            dx, dy, dz = rotate_direction(dx, dy, dz, cos_θ, φ_az)
            E = E_scatt

            # Energy cutoff: deposit remainder locally
            if E < E_cutoff
                handle_deposit!(state, det, params, x, y, z, E)
                if state.outcome === :in_progress
                    finalize_outcome!(state, rng, params)
                end
                break
            end
        end
    end

    return state.outcome
end
