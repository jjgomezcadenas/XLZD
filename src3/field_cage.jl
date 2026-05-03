# src3/field_cage.jl — LZ field cage + grid-holder components.
#
# Six effective sources, all sitting *inside* the LXe (no cryostat-Ti
# downstream), parallel to the cryostat barrel/heads built in cryostat.jl:
#
#   FCRN  — field-cage Ti rings           (barrel,  R=74.3, Ti)
#   FCRS  — field-cage resistors          (barrel,  R=74.3, SS proxy)
#   FCSE  — TPC sensors                   (barrel,  R=74.3, SS proxy)
#   FCPT  — PTFE reflector walls          (barrel,  R=72.8, PTFE)
#   FCTG  — top grids+holders (anode+gate),  endcap-top (annular slab at z_gate)
#   FCBG  — bottom grid+holder (cathode only), endcap-bot (annular slab at z=0)
#
# Data lives in two CSVs (single source of truth):
#   data/lz_fc_barrels.csv  — 4 PCyl rows  (FCRN, FCRS, FCSE, FCPT)
#   data/lz_fc_grids.csv    — 2 PAnnularDisk rows (FCTG, FCBG)
#
# CSV schema is self-contained: every per-component dimension, mass and
# specific activity sits in the row. Materials (Ti, SS, PTFE) are named
# as strings; the loader looks each up in a `materials` Dict supplied by
# the caller, so NIST table paths stay in code (alongside the LXe + Ti
# loaders the cryostat already uses).
#
# Effective wall thickness (PCyl) and slab height (PAnnularDisk) are
# back-derived from mass + density at load time so the slab self-shielding
# integral receives the correct optical depth without ever needing the
# thickness as an explicit CSV column.

# ---------------------------------------------------------------------------
# Composite type
# ---------------------------------------------------------------------------

"""
    FieldCage

Container for the six field-cage PObjects, ordered as they appear in
the source CSVs:

  * `barrels` — Vector{PCyl}        (FCRN, FCRS, FCSE, FCPT)
  * `grids`   — Vector{PAnnularDisk} (FCTG, FCBG)

Built once at startup from `build_field_cage(barrels_csv, grids_csv,
materials)`; passed read-only to the effective-source builder.
"""
struct FieldCage
    barrels::Vector{PCyl}
    grids::Vector{PAnnularDisk}
end

# ---------------------------------------------------------------------------
# CSV → PObject
# ---------------------------------------------------------------------------

"""
    _back_derive_wall_thickness_cm(mass_kg, ρ_g_cm3, R_in_cm, height_cm) -> cm

Back-solve for the cylindrical-shell wall thickness that reproduces a
target mass at fixed inner radius and axial height. Smears the
component's physical metal over the full FC barrel envelope; the slab
self-shielding integral then sees the correct optical depth (μ·t)
without needing an explicit per-component geometry.
"""
function _back_derive_wall_thickness_cm(mass_kg::Real, ρ_g_cm3::Real,
                                        R_in_cm::Real, height_cm::Real)::Float64
    # mass(GCyl) = π·(R_o² − R_i²)·H · ρ / 1000.0 = kg
    # Solve:    (R_in + t)² − R_in² = mass·1000 / (π·H·ρ)
    target_diff = Float64(mass_kg) * 1000.0 / (π * height_cm * ρ_g_cm3)
    R_out = sqrt(R_in_cm^2 + target_diff)
    R_out - R_in_cm
end

"""
    _lookup_material(materials, name) -> Material

Resolve a material name (a CSV cell) to the `Material` object the
caller supplied. Errors with a clear message if the name is unknown,
listing the available keys for fast diagnosis of typos in CSVs.
"""
function _lookup_material(materials::Dict{String,Material},
                          name::AbstractString)::Material
    haskey(materials, name) ||
        error("FieldCage CSV: unknown material '$name'. " *
              "Known: $(sort(collect(keys(materials))))")
    materials[name]
end

"""
    _build_barrel_pcyl(row, materials) -> PCyl

Build one barrel-source PCyl from a parsed CSV row. Wall thickness is
back-derived from mass.
"""
function _build_barrel_pcyl(row::Dict{String,Any},
                            materials::Dict{String,Material})::PCyl
    name    = strip(string(row["name"]))
    mat     = _lookup_material(materials, strip(string(row["material"])))
    R_in    = Float64(row["R_inner_cm"])
    z_min   = Float64(row["z_min_cm"])
    z_max   = Float64(row["z_max_cm"])
    mass_kg = Float64(row["mass_kg"])
    aBi     = Float64(row["act_U238_late_mBqkg"])
    aTh     = Float64(row["act_Th232_late_mBqkg"])

    t = _back_derive_wall_thickness_cm(mass_kg, mat.density, R_in, z_max - z_min)
    g = GCyl(R_in, t, z_min, z_max)
    PCyl(g, mat, aBi, aTh; count=1, name=name)
end

"""
    _build_grid_annulus(row, materials) -> PAnnularDisk

Build one grid-holder PAnnularDisk from a parsed CSV row. Slab height
H is back-derived from mass and density.
"""
function _build_grid_annulus(row::Dict{String,Any},
                             materials::Dict{String,Material})::PAnnularDisk
    name        = strip(string(row["name"]))
    mat         = _lookup_material(materials, strip(string(row["material"])))
    R_in        = Float64(row["R_in_cm"])
    R_out       = Float64(row["R_out_cm"])
    z_face      = Float64(row["z_face_cm"])
    normal_sign = Int(row["normal_sign"])
    mass_kg     = Float64(row["mass_kg"])
    aBi         = Float64(row["act_U238_late_mBqkg"])
    aTh         = Float64(row["act_Th232_late_mBqkg"])

    H = mass_kg * 1000.0 / (mat.density * π * (R_out^2 - R_in^2))
    PAnnularDisk(R_in, R_out, z_face, H, normal_sign, mat, aBi, aTh;
                 count=1, name=name)
end

# ---------------------------------------------------------------------------
# Builder
# ---------------------------------------------------------------------------

"""
    build_field_cage(barrels_csv, grids_csv, materials) -> FieldCage

Construct the six field-cage components by reading the two CSVs:

  * `barrels_csv` (default `data/lz_fc_barrels.csv`) — 4 PCyl rows
  * `grids_csv`   (default `data/lz_fc_grids.csv`)   — 2 PAnnularDisk rows

`materials` is a `Dict{String,Material}` that maps every material name
appearing in the CSVs to a fully-loaded `Material`. Typical wiring:

    mat_Ti   = load_material("Ti",   4.510, "data/nist_ti.csv")
    mat_SS   = load_material("SS",   7.930, "data/nist_fe.csv")
    mat_PTFE = load_material("PTFE", 2.200, "data/nist_teflon.csv")
    fc = build_field_cage("data/lz_fc_barrels.csv",
                           "data/lz_fc_grids.csv",
                           Dict("Ti" => mat_Ti, "SS" => mat_SS,
                                "PTFE" => mat_PTFE))
"""
function build_field_cage(barrels_csv::AbstractString,
                          grids_csv::AbstractString,
                          materials::Dict{String,Material})::FieldCage
    barrels = PCyl[_build_barrel_pcyl(row, materials)
                   for row in _read_csv_rows(barrels_csv)]
    grids   = PAnnularDisk[_build_grid_annulus(row, materials)
                           for row in _read_csv_rows(grids_csv)]
    FieldCage(barrels, grids)
end

# ---------------------------------------------------------------------------
# Iteration / lookup helpers
# ---------------------------------------------------------------------------

"""
    fc_components(fc::FieldCage) -> Vector{PObject}

Return the six FC PObjects in canonical CSV order:
    [FCRN, FCRS, FCSE, FCPT, FCTG, FCBG]
(barrels in lz_fc_barrels.csv order, then grids in lz_fc_grids.csv order)
"""
fc_components(fc::FieldCage)::Vector{PObject} =
    PObject[fc.barrels..., fc.grids...]

"""
    fc_total_mass(fc::FieldCage) -> Float64

Total field-cage mass (kg) summed over the six components.
"""
fc_total_mass(fc::FieldCage)::Float64 = sum(mass(p) for p in fc_components(fc))
