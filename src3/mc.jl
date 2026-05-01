# src3/mc.jl — Per-photon MC helper functions.
#
# After cutover (step 11i) this file contains ONLY:
#   * `path_to_next_region`     — geometric ray/boundary distance
#   * `companion_visible!`      — Tl-208 cascade companion visibility check
#   * `companion_reach_prob`    — companion source per-decay reach probability
#
# The legacy tracker (`track_one_photon!` and helpers — `PhotonState`,
# `handle_deposit!`, `_track_photon_segment!`, `finalize_outcome!`)
# was removed in step 11i. Photon tracking is now in `src3/tracker.jl`
# (`track_photon_stack`, `_track_child_photon!`, `fast_veto`).

using Random: AbstractRNG

# ---------------------------------------------------------------------------
# Boundary finder
# ---------------------------------------------------------------------------

"""
    path_to_next_region(x, y, z, dx, dy, dz, det) -> Float64

Smallest positive distance along the ray from `(x, y, z)` in direction
`(dx, dy, dz)` to ANY of the detector's region-defining surfaces:
  - cylindrical: `r = R_FC_inner, R_FC_outer, R_ICV_inner`
  - planar:      `z = z_gate, z_RFR_bottom, z_LXe_bottom, z_cathode`

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
# Companion veto (Tl-208)
# ---------------------------------------------------------------------------

"""
    companion_visible!(rng, det, comp_eff, xcom, params) -> Bool

Track one companion γ from `comp_eff` (lumped at 583 keV) into the
LXe and return `true` iff it produces a visible deposit anywhere —
i.e. any single deposit ≥ `params.E_visible_keV` in `:TPC` LXe, OR the
accumulated skin deposit ≥ `params.E_skin_veto_keV`. Pair production is
treated as visible (the two annihilation γ guarantee a visible deposit).

Stripped-down: no SS/MS clustering, no FV cut, no ROI smearing — just
a Bool answer for the veto step. Intended to be called only when the
main γ has already been tagged `:SS_in_ROI`.
"""
function companion_visible!(rng::AbstractRNG, det::LXeDetector,
                             comp_eff::EffectiveSource,
                             xcom::XCOMTable,
                             params::MCParams)::Bool
    x, y, z, dx, dy, dz = sample_entry(rng, det, comp_eff)
    E = comp_eff.E_MeV
    E_visible_MeV   = params.E_visible_keV   * 1.0e-3
    E_skin_veto_MeV = params.E_skin_veto_keV * 1.0e-3
    E_cutoff = E_tracking_cutoff_MeV(params)
    ρ = det.material.density
    E_skin_total = 0.0

    while true
        reg = region_at(det, x, y, z)
        if reg === :Outside
            return false
        end

        if reg === :Gas || reg === :FC
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
        if reg_dep === :TPC
            E_dep >= E_visible_MeV && return true
        elseif reg_dep === :Skin
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
