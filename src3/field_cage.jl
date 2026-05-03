# src3/field_cage.jl — LZ field cage + grid-holder components.
#
# Six effective sources, all sitting *inside* the LXe (no cryostat-Ti
# downstream), parallel to the cryostat barrel/heads built in cryostat.jl:
#
#   FCRN  — field-cage Ti rings           (barrel,  R=74.3, Ti)
#   FCRS  — field-cage resistors          (barrel,  R=74.3, SS proxy)
#   FCSE  — TPC sensors                   (barrel,  R=74.3, SS proxy)
#   FCPT  — PTFE reflector walls          (barrel,  R=72.8, PTFE)
#   FCTG  — top grids+holders (anode+gate),  endcap-top (annular slab at z_gate, +H above)
#   FCBG  — bottom grid+holder (cathode only), endcap-bot (annular slab at z=0, +H below)
#
# Masses and specific activities from the LZ bb0nu paper Table I:
#   Field-cage rings        93.0 kg   Bi=0.35  Tl=0.24
#   Field-cage resistors     0.06 kg  Bi=1350  Tl=2010
#   TPC sensors              5.02 kg  Bi=5.82  Tl=1.88
#   PTFE walls               184  kg  Bi=0.04  Tl=0.01
#   Field grids and holders  89.1 kg  Bi=2.63  Tl=1.46  (4-grid lump)
#       → FCTG = 2/4 of 89.1 = 44.55 kg  (anode + gate at top)
#       → FCBG = 1/4 of 89.1 = 22.28 kg  (cathode only;
#                                          bottom-shield grid dropped)
#
# The "barrel" components are smeared cylindrical shells over the full
# field-cage axial extent (drift + RFR). Effective wall thickness is
# back-derived from mass/(ρ·area) so the GCyl carries a physically
# meaningful slab depth for the self-shielding integral.
#
# Grid-holder geometry: thick-walled annulus from R_in=72.8 (TPC inner)
# to R_out=80.3 (skin outer), axial slab thickness H back-derived from
# the SS density. The slab extends *away* from the active LXe (above the
# gate plane for FCTG; below the cathode plane for FCBG); the LXe-facing
# face lies at the gate / cathode plane respectively.

# ---------------------------------------------------------------------------
# Geometry constants
# ---------------------------------------------------------------------------

# TPC anchors (from docs/LZ_detector_summary.md §2)
const FC_R_TPC_INNER_CM   = 72.8         # PTFE inner = TPC inner
const FC_R_RING_CM        = 74.3         # rings, resistors, sensors at FC outer
const FC_R_HOLDER_OUT_CM  = 80.3         # grid-holder outer (= skin outer)
const FC_BARREL_Z_BOT_CM  = -13.75       # bottom of RFR
const FC_BARREL_Z_TOP_CM  =  145.6       # gate plane (top of drift)
const FC_BARREL_HEIGHT_CM = FC_BARREL_Z_TOP_CM - FC_BARREL_Z_BOT_CM  # 159.35

const FC_Z_GATE_CM        = FC_BARREL_Z_TOP_CM   # top grid emission face
const FC_Z_CATHODE_CM     = 0.0                  # bottom grid emission face

# Material densities (g/cm³)
const ρ_Ti_GCM3           = 4.510        # field-cage rings
const ρ_SS_GCM3           = 7.930        # SS-304 (Fe XCOM proxy)
const ρ_PTFE_GCM3         = 2.200        # PTFE reflectors

# bb0nu Table I — masses (kg)
const FC_M_RINGS_KG       = 93.0
const FC_M_RESISTORS_KG   = 0.06
const FC_M_SENSORS_KG     = 5.02
const FC_M_PTFE_KG        = 184.0
const FC_M_GRIDS_TOTAL_KG = 89.1
const FC_M_FCTG_KG        = FC_M_GRIDS_TOTAL_KG / 2          # 44.55
const FC_M_FCBG_KG        = FC_M_GRIDS_TOTAL_KG / 4          # 22.275

# bb0nu Table I — late-chain specific activities (mBq/kg)
const FC_A_RINGS_BI214      = 0.35
const FC_A_RINGS_TL208      = 0.24
const FC_A_RESISTORS_BI214  = 1350.0
const FC_A_RESISTORS_TL208  = 2010.0
const FC_A_SENSORS_BI214    = 5.82
const FC_A_SENSORS_TL208    = 1.88
const FC_A_PTFE_BI214       = 0.04
const FC_A_PTFE_TL208       = 0.01
const FC_A_GRIDS_BI214      = 2.63
const FC_A_GRIDS_TL208      = 1.46

# ---------------------------------------------------------------------------
# Composite type
# ---------------------------------------------------------------------------

"""
    FieldCage

Container for the six field-cage PObjects. Built once at startup from
`build_field_cage(mat_Ti, mat_SS, mat_PTFE)`; passed read-only to the
effective-source builder.

Fields:
  * `FCRN` — Ti rings, barrel, PCyl
  * `FCRS` — resistors, barrel, PCyl (SS proxy)
  * `FCSE` — TPC sensors, barrel, PCyl (SS proxy)
  * `FCPT` — PTFE walls, barrel, PCyl
  * `FCTG` — top grids+holders, PAnnularDisk (single emission face = gate plane)
  * `FCBG` — cathode grid+holder, PAnnularDisk (single emission face = cathode plane)
"""
struct FieldCage
    FCRN::PCyl
    FCRS::PCyl
    FCSE::PCyl
    FCPT::PCyl
    FCTG::PAnnularDisk
    FCBG::PAnnularDisk
end

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

"""
    _back_derive_wall_thickness_cm(mass_kg, ρ_g_cm3, R_in_cm, height_cm) -> cm

Back-solve for the cylindrical-shell wall thickness that reproduces a
target mass at fixed inner radius and axial height. This smears the
component's physical metal over the full FC barrel envelope; the slab
self-shielding integral then sees the correct optical depth (μ·t)
without needing an explicit per-component geometry.
"""
function _back_derive_wall_thickness_cm(mass_kg::Real, ρ_g_cm3::Real,
                                        R_in_cm::Real, height_cm::Real)::Float64
    # mass(GCyl) = π·(R_o² − R_i²)·H · ρ / 1000.0 = kg
    # Solve: (R_in + t)² − R_in² = mass·1000 / (π·H·ρ)
    target_diff = Float64(mass_kg) * 1000.0 / (π * height_cm * ρ_g_cm3)
    R_out = sqrt(R_in_cm^2 + target_diff)
    R_out - R_in_cm
end

"""
    _fcbarrel_pcyl(name, mass_kg, R_in_cm, mat, aBi, aTh) -> PCyl

Build a barrel field-cage PCyl whose wall thickness is back-derived
from the target mass at the supplied inner radius and the FC barrel
axial extent.
"""
function _fcbarrel_pcyl(name::AbstractString, mass_kg::Real, R_in_cm::Real,
                        mat::Material, aBi::Real, aTh::Real)::PCyl
    t = _back_derive_wall_thickness_cm(mass_kg, mat.density,
                                        R_in_cm, FC_BARREL_HEIGHT_CM)
    g = GCyl(R_in_cm, t, FC_BARREL_Z_BOT_CM, FC_BARREL_Z_TOP_CM)
    PCyl(g, mat, aBi, aTh; count=1, name=String(name))
end

"""
    _annular_holder(name, mass_kg, z_face_cm, normal_sign, mat, aBi, aTh) -> PAnnularDisk

Build an annular-slab grid holder (FCTG / FCBG) whose axial thickness
H is back-derived from the target mass, given the SS density and the
fixed annular footprint π·(R_out² − R_in²).
"""
function _annular_holder(name::AbstractString, mass_kg::Real,
                         z_face_cm::Real, normal_sign::Integer,
                         mat::Material, aBi::Real, aTh::Real)::PAnnularDisk
    H = mass_kg * 1000.0 /
        (mat.density * π * (FC_R_HOLDER_OUT_CM^2 - FC_R_TPC_INNER_CM^2))
    PAnnularDisk(FC_R_TPC_INNER_CM, FC_R_HOLDER_OUT_CM,
                 Float64(z_face_cm), H, Int(normal_sign),
                 mat, Float64(aBi), Float64(aTh);
                 count=1, name=String(name))
end

# ---------------------------------------------------------------------------
# Builder
# ---------------------------------------------------------------------------

"""
    build_field_cage(mat_Ti, mat_SS, mat_PTFE) -> FieldCage

Construct the six field-cage components with the bb0nu Table I masses
and activities. Materials must be provided by the caller (so the script
controls the NIST table source); typically:

    mat_Ti   = load_material("Ti",   4.510, "data/nist_ti.csv")
    mat_SS   = load_material("SS",   7.930, "data/nist_fe.csv")
    mat_PTFE = load_material("PTFE", 2.200, "data/nist_teflon.csv")
"""
function build_field_cage(mat_Ti::Material, mat_SS::Material,
                          mat_PTFE::Material)::FieldCage
    fcrn = _fcbarrel_pcyl("FCRN", FC_M_RINGS_KG,     FC_R_RING_CM,
                          mat_Ti,   FC_A_RINGS_BI214,     FC_A_RINGS_TL208)
    fcrs = _fcbarrel_pcyl("FCRS", FC_M_RESISTORS_KG, FC_R_RING_CM,
                          mat_SS,   FC_A_RESISTORS_BI214, FC_A_RESISTORS_TL208)
    fcse = _fcbarrel_pcyl("FCSE", FC_M_SENSORS_KG,   FC_R_RING_CM,
                          mat_SS,   FC_A_SENSORS_BI214,   FC_A_SENSORS_TL208)
    fcpt = _fcbarrel_pcyl("FCPT", FC_M_PTFE_KG,      FC_R_TPC_INNER_CM,
                          mat_PTFE, FC_A_PTFE_BI214,      FC_A_PTFE_TL208)

    # FCTG: anode+gate collapsed at the gate plane. Emission face = z_gate
    # (LXe-facing side of the slab); the slab itself extends UP into the
    # gas+anode region (irrelevant for active-LXe accounting).
    fctg = _annular_holder("FCTG", FC_M_FCTG_KG, FC_Z_GATE_CM, -1,
                           mat_SS, FC_A_GRIDS_BI214, FC_A_GRIDS_TL208)
    # FCBG: cathode only. Emission face = z=0 (LXe-facing top of slab);
    # slab extends DOWN into the RFR (also irrelevant).
    fcbg = _annular_holder("FCBG", FC_M_FCBG_KG, FC_Z_CATHODE_CM, +1,
                           mat_SS, FC_A_GRIDS_BI214, FC_A_GRIDS_TL208)

    FieldCage(fcrn, fcrs, fcse, fcpt, fctg, fcbg)
end

# ---------------------------------------------------------------------------
# Iteration / lookup helpers
# ---------------------------------------------------------------------------

"""
    fc_components(fc::FieldCage) -> Vector{PObject}

Return the six FC PObjects in canonical order:
    [FCRN, FCRS, FCSE, FCPT, FCTG, FCBG]
"""
fc_components(fc::FieldCage)::Vector{PObject} =
    PObject[fc.FCRN, fc.FCRS, fc.FCSE, fc.FCPT, fc.FCTG, fc.FCBG]

"""
    fc_total_mass(fc::FieldCage) -> Float64

Total field-cage mass (kg) summed over the six components. For sanity
checks vs. bb0nu Table I (≈ 350.6 kg with FCBG at 1/4 mass).
"""
fc_total_mass(fc::FieldCage)::Float64 = sum(mass(p) for p in fc_components(fc))
