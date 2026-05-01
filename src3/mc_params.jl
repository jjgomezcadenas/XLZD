# src3/mc_params.jl — Monte Carlo control parameters: Q-value, energy
# resolution, ROI window, FV box, tracking cutoff, SS/MS clustering
# threshold. Detector visibility thresholds (E_visible, E_skin_veto) live
# on LXeDetector, not here.

"""
    MCParams

Per-photon-MC control knobs.

Fields:
  * Physics — `Q_betabeta_keV`, `σ_E_over_E`, `ROI_halfwidth_keV`
  * Tracking — `E_tracking_cutoff_keV`, `Δz_threshold_mm`
  * Fiducial volume box — `fv_z_min_cm`, `fv_z_max_cm`, `fv_r2_max_cm2`

Defaults reproduce the LZ 0νββ analysis (`Q_ββ = 2458 keV`,
σ/E = 0.7 %, ±1σ ROI = ±17.2 keV, FV = 39 cm × [26, 96] cm).
"""
Base.@kwdef struct MCParams
    Q_betabeta_keV::Float64        = 2458.0
    σ_E_over_E::Float64             = 0.007
    ROI_halfwidth_keV::Float64      = 17.2
    E_tracking_cutoff_keV::Float64  = 40.0
    Δz_threshold_mm::Float64        = 3.0
    fv_z_min_cm::Float64            = 26.0
    fv_z_max_cm::Float64            = 96.0
    fv_r2_max_cm2::Float64          = 1521.0       # = 39²
end

"Tracking cutoff (MeV)."
E_tracking_cutoff_MeV(p::MCParams) = p.E_tracking_cutoff_keV / 1000.0

"SS/MS clustering threshold (cm)."
Δz_threshold_cm(p::MCParams) = p.Δz_threshold_mm / 10.0

"Test if a cluster center (x, y, z) lies inside the FV box."
@inline in_fv(x::Real, y::Real, z::Real, p::MCParams)::Bool =
    z >= p.fv_z_min_cm && z <= p.fv_z_max_cm &&
    (x*x + y*y) <= p.fv_r2_max_cm2

"""
    classify_ss_energy(rng, E_cluster_MeV, p::MCParams) -> Symbol

Apply Gaussian energy smearing (σ/E = `p.σ_E_over_E`) to the cluster
energy and check the ±`p.ROI_halfwidth_keV` window around `p.Q_betabeta_keV`.
Returns `:SS_in_ROI` or `:SS_outside_ROI`.
"""
function classify_ss_energy(rng::AbstractRNG, E_cluster_MeV::Float64,
                             p::MCParams)::Symbol
    σ_E = p.σ_E_over_E * E_cluster_MeV
    E_smeared_keV = (E_cluster_MeV + σ_E * randn(rng)) * 1000.0
    if abs(E_smeared_keV - p.Q_betabeta_keV) <= p.ROI_halfwidth_keV
        return :SS_in_ROI
    else
        return :SS_outside_ROI
    end
end
