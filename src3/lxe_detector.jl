# src3/lxe_detector.jl — LXe detector geometry, region classifier.
#
# Encapsulates the LXe-side of the detector model:
#   :TPC      active TPC LXe (visible)
#   :Skin     LXe skin (PMT-readout veto layer)
#   :Inert    inert LXe (RFR + dome — invisible)
#   :FC       field-cage annulus (1.5 cm PTFE+Ti — treated as vacuum)
#   :Gas      gas (above the liquid surface — no interactions)
#   :Outside  outside the LXe envelope
#
# `region_at(det, x, y, z)` returns one of these as a Symbol. The
# visibility thresholds (E_visible_keV, E_skin_veto_keV) are physics
# policy and live on MCParams, not here.

"""
    LXeDetector

LZ LXe detector geometry + the LXe `Material`. Lengths in cm. The
visibility / veto thresholds live on `MCParams`.
"""
struct LXeDetector
    # Geometry (cm)
    R_FC_inner::Float64        # active TPC outer radius
    R_FC_outer::Float64        # field cage outer = skin inner
    R_ICV_inner::Float64       # ICV inner = skin outer
    z_cathode::Float64
    z_gate::Float64            # = LXe surface
    z_RFR_bottom::Float64      # bottom of the field-cage structure
    z_LXe_bottom::Float64      # ICV bottom apex
    z_ICV_top::Float64         # ICV top apex (upper gas bound)
    # Physics
    material::Material         # LXe; carries ρ and μ_lin(E)
end

# ---------------------------------------------------------------------------
# Loader
# ---------------------------------------------------------------------------

"""
    build_lxe_detector(csv_path, mat_LXe) -> LXeDetector

Read the single-row CSV at `csv_path` and combine with `mat_LXe` (the
`Material` returned by `load_material("LXe", ρ, "data/nist_lxe.csv")`).
The CSV's `E_visible_keV` and `E_skin_veto_keV` columns are ignored
(thresholds now live on MCParams); the columns may remain in the file
for historical reference.
"""
function build_lxe_detector(csv_path::AbstractString,
                            mat_LXe::Material)::LXeDetector
    rows = _read_csv_rows(csv_path)
    length(rows) == 1 ||
        error("build_lxe_detector: expected 1 row, got $(length(rows))")
    r = rows[1]
    LXeDetector(
        Float64(r["R_FC_inner_cm"]),
        Float64(r["R_FC_outer_cm"]),
        Float64(r["R_ICV_inner_cm"]),
        Float64(r["z_cathode_cm"]),
        Float64(r["z_gate_cm"]),
        Float64(r["z_RFR_bottom_cm"]),
        Float64(r["z_LXe_bottom_cm"]),
        Float64(r["z_ICV_top_cm"]),
        mat_LXe,
    )
end

# ---------------------------------------------------------------------------
# Region classifier
# ---------------------------------------------------------------------------

"""
    region_at(det, x, y, z) -> Symbol

Classify a point as one of `:Gas`, `:Outside`, `:Skin`, `:FC`,
`:TPC`, `:Inert`. The FC annulus (`R_FC_inner < r < R_FC_outer`) is
treated as transparent (vacuum) per project convention; the per-photon
MC handles `:FC` exactly like `:Gas`.
"""
function region_at(det::LXeDetector, x::Real, y::Real, z::Real)::Symbol
    r2 = x*x + y*y
    z > det.z_gate                                              && return :Gas
    r2 > det.R_ICV_inner^2                                       && return :Outside

    if z >= det.z_RFR_bottom
        # In the FC-bearing z range
        r2 > det.R_FC_outer^2  && return :Skin
        r2 > det.R_FC_inner^2  && return :FC
        # r ≤ R_FC_inner
        z >= det.z_cathode      && return :TPC   # drift region
        return :Inert                               # RFR (below cathode)
    else
        # Below z_RFR_bottom — dome LXe (no FC structure here)
        return :Inert
    end
end

# ---------------------------------------------------------------------------
# Convenience accessors / sanity helpers
# ---------------------------------------------------------------------------

"LXe linear attenuation coefficient (cm⁻¹) at `E_MeV`."
μ_LXe(det::LXeDetector, E_MeV::Real) = det.material.μ_lin(E_MeV)

"Active TPC volume π·R_FC_inner²·(z_gate − z_cathode) (cm³)."
active_volume_cm3(det::LXeDetector) =
    π * det.R_FC_inner^2 * (det.z_gate - det.z_cathode)

"Active TPC LXe mass (kg)."
active_mass_kg(det::LXeDetector) =
    active_volume_cm3(det) * det.material.density / 1000.0

"LXe skin volume π·(R_ICV_inner² − R_FC_outer²)·(z_gate − z_RFR_bottom) (cm³)."
skin_volume_cm3(det::LXeDetector) =
    π * (det.R_ICV_inner^2 - det.R_FC_outer^2) *
        (det.z_gate - det.z_RFR_bottom)

"LXe skin mass (kg)."
skin_mass_kg(det::LXeDetector) =
    skin_volume_cm3(det) * det.material.density / 1000.0
