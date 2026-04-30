# geometry.jl — Detector geometry, source configuration, entry-point sampling

"""
    Params

All user-tunable inputs to the XLZD background MC.
Immutable. Constructed once at startup from CLI args or defaults.
Fields are prefix-grouped: geom_*, phys_*, cut_*, fv_*, mc_*, hist_*, out_*.

Ring-specific fields (geom_pitch, act_Cu, etc.) have been removed.
Source configuration is now handled by SourceConfig.
"""
Base.@kwdef struct Params
    # ─── Geometry (geom_*) — LZ TPC dimensions ─────────────────────────
    geom_R_lxe::Float64        = 72.8     # TPC inner radius (cm), 1456 mm diameter
    geom_L_lxe::Float64        = 145.6    # TPC drift region height (cm), cathode to gate
    geom_ρ_LXe::Float64        = 2.953    # LXe density at saturation, ~165 K (g/cm³)

    # ─── Extended tracking volume (for PMT sources) ──────────────────
    # The tracking cylinder extends below the cathode to include the RFR
    # and the space down to the bottom PMTs. Interactions below the
    # cathode (z < 0 in drift coords) are invisible (no S2 signal).
    # Top PMTs are in gas above the liquid surface — no LXe above gate.
    geom_z_bottom_pmt::Float64 = -15.75   # cm, bottom PMT faces (2 cm below shield grid)
    geom_z_top_pmt::Float64    = 152.6    # cm, top PMT faces (7 cm above gate, in gas)
    geom_z_cathode::Float64    = 0.0      # cm, cathode position (z=0 = bottom of drift)
    geom_z_rfr_bottom::Float64 = -13.75   # cm, bottom of RFR (= bottom shield grid)

    # ─── Physics constants and tables (phys_*) ──────────────────────────
    phys_Q_betabeta_keV::Float64   = 2458.0  # ¹³⁶Xe Q value
    phys_σ_E_over_E::Float64       = 0.007   # energy resolution σ/E at Q_ββ (0.7%)
    phys_xcom_data_path::String    = "data/nist.csv"

    # ─── Event-selection cuts (cut_*) ───────────────────────────────────
    # ROI: ±1σ around Q_ββ. At σ/E=0.7%: σ_E = 17.2 keV → halfwidth = 17.2 keV
    cut_ROI_halfwidth_keV::Float64       = 17.2      # ±1σ around Q_ββ
    cut_Δz_threshold_mm::Float64         = 3.0       # SS/MS threshold in z
    cut_E_visible_threshold_keV::Float64 = 5.0       # below this, deposit invisible
    cut_E_tracking_cutoff_keV::Float64   = 40.0      # below this, photon stops

    # ─── Companion gamma veto (for ²⁰⁸Tl two-gamma events) ───────────
    cut_companion_veto_keV::Float64      = 5.0       # min deposit in active LXe to tag companion
    cut_skin_veto_keV::Float64           = 100.0     # min deposit in LXe skin to veto event

    # ─── Skin geometry (for companion veto) ───────────────────────────
    # Cylindrical shell, 6 cm average thickness, same height as TPC
    geom_R_skin_inner::Float64   = 74.3    # cm, FC outer wall = skin inner boundary
    geom_R_skin_outer::Float64   = 80.3    # cm, 74.3 + 6.0 cm skin

    # ─── Fiducial volume (fv_*) — from LZ 0νββ paper (967 kg inner vol) ─
    fv_z_min_cm::Float64             = 26.0     # cm
    fv_z_max_cm::Float64             = 96.0     # cm
    fv_r2_max_cm2::Float64           = 1521.0   # = 39² cm²

    # ─── Monte Carlo control (mc_*) ─────────────────────────────────────
    mc_N_samples::Int                = 100_000_000
    mc_seed::Int                     = 1234
    mc_n_traj_per_outcome::Int       = 100

    # ─── Histogram binning (hist_*) ─────────────────────────────────────
    hist_z_min_cm::Float64         = 0.0
    hist_z_max_cm::Float64         = 145.6        # = L_lxe
    hist_r2_min_cm2::Float64       = 0.0
    hist_r2_max_cm2::Float64       = 5300.0       # ≈ 72.8² cm²
    hist_n_z_bins::Int             = 100
    hist_n_r2_bins::Int            = 100
    hist_n_E_bins::Int             = 250
    hist_E_max_keV::Float64        = 2700.0       # cover Tl-208 line at 2615 keV
    hist_r2_inner_max_cm2::Float64 = 2500.0  # = 50² cm², zoomed view around FV (r≤39)

    # ─── Output (out_*) ─────────────────────────────────────────────────
    out_dir::String                = "output/last_run/"
    out_viewer_template::String    = "viewer/viewer.html"
end

"""
    copy_params(p::Params; kwargs...) -> Params

Create a copy of Params with selected fields overridden.
Julia @kwdef structs don't support copy-with-modification natively.
"""
function copy_params(p::Params; kwargs...)
    fields = Dict{Symbol, Any}()
    for f in fieldnames(Params)
        fields[f] = getfield(p, f)
    end
    for (k, v) in kwargs
        fields[k] = v
    end
    Params(; fields...)
end

"""
    BFV

Buffer Fiducial Volume — FV expanded by 1 mm on all sides.
Events whose first interaction falls outside the BFV are immediately
rejected (no tracking needed). Events inside the BFV are tracked, but
only count as signal if the final cluster is inside the FV.
"""
struct BFV
    z_min::Float64
    z_max::Float64
    r2_max::Float64
end

"""
    build_bfv(params::Params) -> BFV

Construct the BFV from FV parameters + 1 mm buffer.
"""
function build_bfv(params::Params)::BFV
    buffer = 0.1  # cm = 1 mm
    BFV(
        params.fv_z_min_cm - buffer,
        params.fv_z_max_cm + buffer,
        (sqrt(params.fv_r2_max_cm2) + buffer)^2
    )
end

"""
    in_bfv(x, y, z, bfv::BFV) -> Bool

Test if a point is inside the Buffer Fiducial Volume.
"""
@inline function in_bfv(x::Float64, y::Float64, z::Float64, bfv::BFV)::Bool
    z >= bfv.z_min && z <= bfv.z_max && (x*x + y*y) <= bfv.r2_max
end

"""
    in_fv(x, y, z, params::Params) -> Bool

Test if a point is inside the Fiducial Volume.
"""
@inline function in_fv(x::Float64, y::Float64, z::Float64, params::Params)::Bool
    z >= params.fv_z_min_cm && z <= params.fv_z_max_cm &&
        (x*x + y*y) <= params.fv_r2_max_cm2
end

"""
    in_skin(x, y, z, params::Params, L_lxe::Float64) -> Bool

Test if a point is in the LXe skin region: inside the LXe cylinder
but outside the field cage (R_skin_inner < r < R_skin_outer).
Also includes the dome region below the TPC.
"""
@inline function in_skin(x::Float64, y::Float64, z::Float64,
                          params::Params, L_lxe::Float64)::Bool
    r2 = x*x + y*y
    R_inner2 = params.geom_R_skin_inner^2
    R_outer2 = params.geom_R_skin_outer^2
    # Barrel skin: between FC outer wall and ICV inner wall
    if r2 > R_inner2 && r2 <= R_outer2 && z >= 0.0 && z <= L_lxe
        return true
    end
    # Dome region: below the TPC (z < ~some threshold), inside FC radius
    # For simplicity, everything inside the LXe cylinder but outside
    # the active region at r < R_skin_inner is considered "active LXe"
    return false
end

"""
    in_active_lxe(x, y, z, params::Params, L_lxe::Float64) -> Bool

Test if a point is in the active LXe (inside the field cage, not skin).
"""
@inline function in_active_lxe(x::Float64, y::Float64, z::Float64,
                                params::Params, L_lxe::Float64)::Bool
    r2 = x*x + y*y
    R_inner2 = params.geom_R_skin_inner^2
    return r2 <= R_inner2 && z >= 0.0 && z <= L_lxe
end

"""
    Geometry

Derived detector geometry, computed once from `Params` at startup.
LXe cylinder dimensions and fiducial volume masses.
"""
struct Geometry
    R_lxe::Float64         # TPC inner radius (cm) — primary tracking boundary
    L_lxe::Float64         # TPC drift height (cm), cathode to gate
    R_lxe_full::Float64    # Full LXe cylinder radius (cm) — includes skin, for companion
    L_lxe_full::Float64    # Full LXe cylinder height (cm) — same as TPC (skin is radial)
    ρ_LXe::Float64
    LXe_mass_kg::Float64
    FV_LXe_mass_kg::Float64
    bfv::BFV
    # Extended tracking volume (for PMT sources)
    z_bottom::Float64      # bottom of extended tracking volume (bottom PMT faces)
    z_top::Float64         # top of extended tracking volume (gate = liquid surface)
    z_cathode::Float64     # cathode z position (below = invisible)
    L_extended::Float64    # total height of extended tracking: z_top - z_bottom
end

"""
    build_geometry(params::Params) -> Geometry

Derive LXe cylinder, fiducial volume, and BFV from `Params`.
"""
function build_geometry(params::Params)::Geometry
    R     = params.geom_R_lxe
    L     = params.geom_L_lxe
    R_full = params.geom_R_skin_outer
    L_full = L
    ρ     = params.geom_ρ_LXe

    # TPC active volume
    TPC_volume_cm3 = π * R^2 * L
    TPC_mass_kg    = TPC_volume_cm3 * ρ / 1000.0

    # Full LXe (TPC + skin)
    full_volume_cm3 = π * R_full^2 * L_full
    full_mass_kg    = full_volume_cm3 * ρ / 1000.0

    # Skin mass
    skin_mass_kg = full_mass_kg - TPC_mass_kg

    r_fv = sqrt(params.fv_r2_max_cm2)
    FV_volume_cm3  = π * r_fv^2 * (params.fv_z_max_cm - params.fv_z_min_cm)
    FV_LXe_mass_kg = FV_volume_cm3 * ρ / 1000.0

    bfv = build_bfv(params)

    # Extended tracking volume for PMT sources
    z_bottom = params.geom_z_bottom_pmt
    z_top = params.geom_L_lxe       # gate = liquid surface = top of LXe
    z_cathode = params.geom_z_cathode
    L_extended = z_top - z_bottom

    @printf("── Geometry summary ──\n")
    @printf("  TPC:  R = %.1f cm, H = %.1f cm, mass = %.1f t\n", R, L, TPC_mass_kg / 1000.0)
    @printf("  Skin: R_outer = %.1f cm, thickness = %.1f cm, mass = %.1f t\n",
            R_full, R_full - params.geom_R_skin_inner, skin_mass_kg / 1000.0)
    @printf("  Full LXe: R = %.1f cm, H = %.1f cm, mass = %.1f t\n",
            R_full, L_full, full_mass_kg / 1000.0)
    @printf("  FV:  z ∈ [%.1f, %.1f] cm, r ≤ %.1f cm, mass = %.2f t\n",
            params.fv_z_min_cm, params.fv_z_max_cm, r_fv, FV_LXe_mass_kg / 1000.0)
    @printf("  BFV: z ∈ [%.1f, %.1f] cm, r ≤ %.1f cm\n",
            bfv.z_min, bfv.z_max, sqrt(bfv.r2_max))
    @printf("  PMT top:    z = %.1f cm (in gas, no LXe above gate)\n", params.geom_z_top_pmt)
    @printf("  PMT bottom: z = %.1f cm (below RFR, in LXe)\n", z_bottom)
    @printf("  Tracking:   z ∈ [%.1f, %.1f] cm (%.1f cm total)\n",
            z_bottom, z_top, L_extended)
    @printf("  Cathode:    z = %.1f cm (below = invisible, no S2)\n", z_cathode)
    @printf("──────────────────────\n")

    Geometry(R, L, R_full, L_full, ρ, full_mass_kg, FV_LXe_mass_kg, bfv,
             z_bottom, z_top, z_cathode, L_extended)
end

# =========================================================================
# Source Configuration
# =========================================================================

"""
    SourceConfig

Defines a gamma source for the MC. Specifies the entry surface,
initial energy, and angular distribution of gammas entering the LXe.

Fields:
  label        — human-readable name (e.g., "cryo_barrel_Bi214")
  E_MeV        — initial gamma energy (2.448 for Bi-214, 2.615 for Tl-208)
  entry        — :barrel or :endcap
  R_entry      — radius of the entry surface (cm)
  H_entry      — height of barrel entry surface (cm); ignored for endcap
  z_min_entry  — z start of barrel surface; ignored for endcap
  angular      — :flat or :shaped
  a, b         — linear fit coefficients for shaped: p(u) = (a + b*u)/norm
  u_min        — lower cutoff for shaped sampling (default 0.3)
  norm         — normalization of linear fit over [u_min, 1.0]
  gammas_per_yr — total gamma rate entering active volume (for normalization)
"""
# ²⁰⁸Tl decay companion gamma branching ratios.
# The 2.6145 MeV gamma is always emitted (99.75%).
# Companion gamma energies and probabilities:
#   583 keV (85%), 860 keV (12%), 763 keV (2%), none (1%)
# The remaining 1% has no companion (or negligible energy).
const TL208_COMPANION_ENERGIES = [0.583, 0.860, 0.763]  # MeV
const TL208_COMPANION_PROBS    = [0.85,  0.12,  0.02]   # sum = 0.99

Base.@kwdef struct SourceConfig
    label::String
    E_MeV::Float64
    entry::Symbol              # :barrel, :endcap, :endcap_top, :endcap_bottom
    R_entry::Float64           # cm
    H_entry::Float64 = 0.0     # cm (barrel height; 0 for endcap)
    z_min_entry::Float64 = 0.0 # cm (barrel z start)
    z_entry::Float64 = 0.0     # cm (fixed z for PMT endcaps)
    angular::Symbol = :flat    # :flat or :shaped
    a::Float64 = 0.0           # linear fit intercept
    b::Float64 = 0.0           # linear fit slope
    u_min::Float64 = 0.3       # shaped lower cutoff
    norm::Float64 = 1.0        # linear fit normalization
    gammas_per_yr::Float64 = 0.0
    is_Tl208::Bool = false     # if true, emit companion gamma
    use_extended_volume::Bool = false  # if true, track in extended cylinder (PMT sources)
end

"""
    sample_companion_energy(rng) -> Float64

Sample the companion gamma energy for a ²⁰⁸Tl decay.
Returns 0.0 if no companion is emitted (1% of decays).
"""
function sample_companion_energy(rng::AbstractRNG)::Float64
    r = rand(rng)
    cumul = 0.0
    for (E, p) in zip(TL208_COMPANION_ENERGIES, TL208_COMPANION_PROBS)
        cumul += p
        if r < cumul
            return E
        end
    end
    return 0.0  # no companion (1% case)
end

# =========================================================================
# Entry-point and direction sampling
# =========================================================================

"""
    sample_entry_point(rng, source::SourceConfig, L_lxe::Float64) -> (x, y, z)

Sample one entry point uniformly on the source's entry surface.

Entry types:
  :barrel        — cylindrical shell at R_entry, z ∈ [z_min, z_min + H]
  :endcap        — disc at R_entry, 50/50 top (z=L_lxe) or bottom (z=0)
  :endcap_top    — disc at R_entry, fixed z = z_entry (PMT top, in gas)
  :endcap_bottom — disc at R_entry, fixed z = z_entry (PMT bottom, below RFR)
"""
function sample_entry_point(rng::AbstractRNG, source::SourceConfig, L_lxe::Float64)
    φ = 2π * rand(rng)

    if source.entry == :barrel
        r = source.R_entry
        z = source.z_min_entry + source.H_entry * rand(rng)
        x = r * cos(φ)
        y = r * sin(φ)
        return (x, y, z)
    elseif source.entry == :endcap_top || source.entry == :endcap_bottom
        # PMT endcap: uniform on disc at fixed z
        r = source.R_entry * sqrt(rand(rng))
        x = r * cos(φ)
        y = r * sin(φ)
        z = source.z_entry
        return (x, y, z)
    else  # :endcap (generic, 50/50 top/bottom)
        r = source.R_entry * sqrt(rand(rng))
        x = r * cos(φ)
        y = r * sin(φ)
        z = rand(rng) < 0.5 ? 0.0 : L_lxe
        return (x, y, z)
    end
end

"""
    sample_u(rng, source::SourceConfig) -> Float64

Sample u = cos θ from the angular distribution.

Flat: u ~ Uniform(0, 1).
Shaped: inverse CDF of the linear fit p(u) = (a + b*u)/norm over [u_min, 1.0].
"""
function sample_u(rng::AbstractRNG, source::SourceConfig)::Float64
    if source.angular == :flat
        return rand(rng)
    else
        # Inverse CDF: solve F(u) = r for the linear PDF p(u) = (a + b*u)/norm
        # C = r * norm + a * u_min + b/2 * u_min^2
        # u = (-a + sqrt(a^2 + 2*b*C)) / b
        r = rand(rng)
        a = source.a
        b = source.b
        u_min = source.u_min
        C = r * source.norm + a * u_min + b / 2.0 * u_min^2
        if abs(b) < 1e-15
            return u_min + r * (1.0 - u_min)
        end
        return (-a + sqrt(a^2 + 2.0 * b * C)) / b
    end
end

"""
    sample_entry_direction(rng, source::SourceConfig, x, y, z, L_lxe) -> (dx, dy, dz)

Sample an inward-directed unit vector at the entry point.

u = cos θ is sampled from the source's angular distribution (flat or shaped).
θ is measured from the inward normal. φ_dir is uniform in [0, 2π).

For barrel: inward normal is -r̂ = (-cos φ_pos, -sin φ_pos, 0).
For endcap: inward normal is -ẑ (top) or +ẑ (bottom).

The direction is constructed by rotating the normal by θ and φ_dir.
"""
function sample_entry_direction(rng::AbstractRNG, source::SourceConfig,
                                 x::Float64, y::Float64, z::Float64,
                                 L_lxe::Float64)
    u = sample_u(rng, source)
    sin_theta = sqrt(max(0.0, 1.0 - u * u))
    φ_dir = 2π * rand(rng)
    cos_φ = cos(φ_dir)
    sin_φ = sin(φ_dir)

    if source.entry == :barrel
        # Inward normal = -r̂ at the entry point
        r = sqrt(x * x + y * y)
        if r < 1e-10
            # Degenerate case (shouldn't happen for barrel)
            return (0.0, 0.0, -1.0)
        end
        # Normal pointing inward: n = (-x/r, -y/r, 0)
        nx = -x / r;  ny = -y / r;  nz = 0.0

        # Build orthogonal basis: n is in xy-plane, so ẑ is perpendicular
        # u1 = ẑ = (0, 0, 1)
        # u2 = n × u1 = (-ny, nx, 0) (tangential direction)
        u1x = 0.0;  u1y = 0.0;  u1z = 1.0
        u2x = -ny;  u2y = nx;   u2z = 0.0

        dx = u * nx + sin_theta * cos_φ * u2x + sin_theta * sin_φ * u1x
        dy = u * ny + sin_theta * cos_φ * u2y + sin_theta * sin_φ * u1y
        dz = u * nz + sin_theta * cos_φ * u2z + sin_theta * sin_φ * u1z

    elseif source.entry == :endcap_top
        # PMT top: inward = -ẑ (downward into LXe)
        dx = sin_theta * cos_φ
        dy = sin_theta * sin_φ
        dz = -u

    elseif source.entry == :endcap_bottom
        # PMT bottom: inward = +ẑ (upward into LXe)
        dx = sin_theta * cos_φ
        dy = sin_theta * sin_φ
        dz = u

    else  # :endcap (generic, 50/50)
        if z < 1.0  # bottom endcap
            dx = sin_theta * cos_φ
            dy = sin_theta * sin_φ
            dz = u
        else  # top endcap
            dx = sin_theta * cos_φ
            dy = sin_theta * sin_φ
            dz = -u
        end
    end

    # Normalize for safety
    norm = sqrt(dx * dx + dy * dy + dz * dz)
    return (dx / norm, dy / norm, dz / norm)
end

"""
    path_to_cylinder_exit(x, y, z, dx, dy, dz, R, z_min, z_max) -> Float64

Ray–cylinder intersection: path length until the photon exits
via side wall, top (z_max), or bottom (z_min), whichever comes first.
Returns 0.0 if the ray misses or starts outside and heading away.
"""
function path_to_cylinder_exit(x::Float64, y::Float64, z::Float64,
                               dx::Float64, dy::Float64, dz::Float64,
                               R::Float64, z_min::Float64, z_max::Float64)::Float64
    EPS = 1.0e-6

    # --- Radial (infinite cylinder) intersection ---
    a = dx * dx + dy * dy
    b = 2.0 * (x * dx + y * dy)
    c = x * x + y * y - R * R

    t_radial_exit = Inf
    if a > EPS * EPS
        disc = b * b - 4.0 * a * c
        if disc < 0.0
            return 0.0
        end
        sq = sqrt(disc)
        t1 = (-b - sq) / (2.0 * a)
        t2 = (-b + sq) / (2.0 * a)
        t_radial_exit = max(t2, EPS)
    else
        if c > EPS
            return 0.0
        end
        t_radial_exit = Inf
    end

    # --- Axial caps ---
    t_bot = Inf
    t_top = Inf
    if dz < -EPS
        t_bot = (z_min - z) / dz
    end
    if dz > EPS
        t_top = (z_max - z) / dz
    end

    d_boundary = min(t_radial_exit, t_bot, t_top)
    return max(d_boundary, 0.0)
end
