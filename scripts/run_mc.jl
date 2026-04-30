using XLZD, ArgParse, DelimitedFiles, Printf

# =========================================================================
# Source definitions — hard-wired from Python CSV output files
# =========================================================================
# These are the 8 distinct MC runs needed (see design/lz_cryostat_model.md):
#   4 shaped (cryostat): barrel/endcap × Bi214/Tl208
#   4 flat (field cage + PMTs): barrel/endcap × Bi214/Tl208
#
# The gammas_per_yr values come from:
#   data/lz_cryo_linear_fit.csv  (cryostat, shaped)
#   data/lz_fc_sampling.csv      (field cage, flat)
#   data/lz_pmt_sampling.csv     (PMTs, flat)
#
# For flat sources, the total rate is the SUM of field cage + PMT rates
# for the same (isotope, region), since they share the same entry geometry
# and angular distribution.
# =========================================================================

# Entry surface parameters (LZ geometry)
const R_TPC        = 72.8    # cm — TPC inner radius
const H_TPC        = 145.6   # cm — TPC drift region height
const Z_MIN_BARREL = 0.0     # cm — barrel starts at cathode (z=0)
const L_LXE        = 145.6   # cm — TPC height
# PMT positions (from TDR)
const Z_PMT_TOP    = 152.6   # cm — top PMT faces (7 cm above gate, in gas)
const Z_PMT_BOTTOM = -15.75  # cm — bottom PMT faces (2 cm below shield grid)
const Z_HEAD_TOP   = 190.0   # cm — ICV top head (2:1 ellipsoidal, in gas)
const Z_HEAD_BOT   = -44.0   # cm — ICV bottom head (3:1 ellipsoidal, in LXe)
const R_FC_OUTER   = 74.3    # cm — FC outer wall (rings, resistors)
const R_ICV_INNER  = 82.1    # cm — ICV inner wall (83.0 - 0.9 cm wall)

function build_sources()
    # ─── Cryostat (shaped, from linear fit) ─────────────────────────
    # Values from data/lz_cryo_linear_fit.csv (updated to bbonu paper)
    # Cryostat gammas enter at the ICV inner wall (R = 80.9 cm).
    # They have been attenuated through Ti only (no LXe skin in Python).
    # The Julia MC tracks them through the LXe skin (6.6 cm) + TPC.
    # use_extended_volume = true so tracking includes the skin region.
    # Cryostat with real wall thicknesses (9/7 mm).
    # Extra 1160 kg assigned to bottom endcap (heavily shielded).
    # The endcap PDF includes both heads + extra mass combined.
    cryo_sources = [
        SourceConfig(
            label = "cryo_barrel_Bi214",
            E_MeV = 2.448, entry = :barrel, angular = :shaped,
            R_entry = R_ICV_INNER, H_entry = H_TPC, z_min_entry = Z_MIN_BARREL,
            a = 0.881389, b = 0.507266, u_min = 0.3, norm = 0.847778,
            gammas_per_yr = 3.382267e4,
            use_extended_volume = true,
        ),
        # Endcap: top and bottom heads at actual positions.
        # Top gets half the head γ/yr, bottom gets half + extra mass γ/yr.
        # From Python: total endcap = 2.986e4, heads = 7.39e3, extra = 2.25e4
        # Top head ≈ 3.69e3, bottom head + extra ≈ 2.62e4
        SourceConfig(
            label = "cryo_head_top_Bi214",
            E_MeV = 2.448, entry = :endcap_top, angular = :shaped,
            R_entry = R_ICV_INNER, z_entry = Z_HEAD_TOP,
            a = 0.959574, b = 0.156725, u_min = 0.3, norm = 0.743012,
            gammas_per_yr = 3.69e3,
            use_extended_volume = true,
        ),
        SourceConfig(
            label = "cryo_head_bot_Bi214",
            E_MeV = 2.448, entry = :endcap_bottom, angular = :shaped,
            R_entry = R_ICV_INNER, z_entry = Z_HEAD_BOT,
            a = 0.959574, b = 0.156725, u_min = 0.3, norm = 0.743012,
            gammas_per_yr = 2.62e4,
            use_extended_volume = true,
        ),
        SourceConfig(
            label = "cryo_barrel_Tl208",
            E_MeV = 2.615, entry = :barrel, angular = :shaped,
            R_entry = R_ICV_INNER, H_entry = H_TPC, z_min_entry = Z_MIN_BARREL,
            a = 0.893791, b = 0.465098, u_min = 0.3, norm = 0.837273,
            gammas_per_yr = 1.046081e6,
            is_Tl208 = true,
            use_extended_volume = true,
        ),
        SourceConfig(
            label = "cryo_head_top_Tl208",
            E_MeV = 2.615, entry = :endcap_top, angular = :shaped,
            R_entry = R_ICV_INNER, z_entry = Z_HEAD_TOP,
            a = 0.961187, b = 0.153947, u_min = 0.3, norm = 0.742877,
            gammas_per_yr = 2.40e5,
            is_Tl208 = true,
            use_extended_volume = true,
        ),
        SourceConfig(
            label = "cryo_head_bot_Tl208",
            E_MeV = 2.615, entry = :endcap_bottom, angular = :shaped,
            R_entry = R_ICV_INNER, z_entry = Z_HEAD_BOT,
            a = 0.961187, b = 0.153947, u_min = 0.3, norm = 0.742877,
            gammas_per_yr = 1.67e6,
            is_Tl208 = true,
            use_extended_volume = true,
        ),
    ]

    # ─── Field cage barrel: split into 3 sub-sources ───────────────
    # Rates from data/lz_fc_gammas.csv (gammas_entering_per_yr)
    # 1. Rings (Ti, 93 kg) at R = 74.3 cm (FC outer wall)
    # 2. PTFE reflectors (184 kg) at R = 72.8 cm (TPC inner wall)
    # 3. Resistors + sensors (0.06 + 5.02 kg) at R = 74.3 cm

    fc_barrel_sources = [
        # Rings
        SourceConfig(
            label = "fc_rings_Bi214",
            E_MeV = 2.448, entry = :barrel, angular = :flat,
            R_entry = R_FC_OUTER, H_entry = H_TPC, z_min_entry = Z_MIN_BARREL,
            gammas_per_yr = 7.96e3,
        ),
        SourceConfig(
            label = "fc_rings_Tl208",
            E_MeV = 2.615, entry = :barrel, angular = :flat,
            R_entry = R_FC_OUTER, H_entry = H_TPC, z_min_entry = Z_MIN_BARREL,
            gammas_per_yr = 1.26e5,
            is_Tl208 = true,
        ),
        # PTFE reflectors (at TPC inner wall)
        SourceConfig(
            label = "fc_ptfe_Bi214",
            E_MeV = 2.448, entry = :barrel, angular = :flat,
            R_entry = R_TPC, H_entry = H_TPC, z_min_entry = Z_MIN_BARREL,
            gammas_per_yr = 1.80e3,
        ),
        SourceConfig(
            label = "fc_ptfe_Tl208",
            E_MeV = 2.615, entry = :barrel, angular = :flat,
            R_entry = R_TPC, H_entry = H_TPC, z_min_entry = Z_MIN_BARREL,
            gammas_per_yr = 1.04e4,
            is_Tl208 = true,
        ),
        # Resistors + sensors
        SourceConfig(
            label = "fc_ressens_Bi214",
            E_MeV = 2.448, entry = :barrel, angular = :flat,
            R_entry = R_FC_OUTER, H_entry = H_TPC, z_min_entry = Z_MIN_BARREL,
            gammas_per_yr = 2.70e4,  # 1.98e4 (resistors) + 7.15e3 (sensors)
        ),
        SourceConfig(
            label = "fc_ressens_Tl208",
            E_MeV = 2.615, entry = :barrel, angular = :flat,
            R_entry = R_FC_OUTER, H_entry = H_TPC, z_min_entry = Z_MIN_BARREL,
            gammas_per_yr = 7.37e5,  # 6.83e5 (resistors) + 5.35e4 (sensors)
            is_Tl208 = true,
        ),
    ]

    # ─── Grid holder rings (at TPC perimeter, short barrel sources) ──
    # The holder rings (88.3 kg total, 44.15 each top/bottom) are ring-shaped
    # frames at the outer edge of the TPC. Modeled as short barrel sources
    # at R = 72.8 cm, H = 2.5 cm, centered at z=0 (cathode) and z=145.6 (gate).
    # Wire mesh mass (0.8 kg) is negligible and dropped.
    # Rates: 44.15 kg × 2.63 mBq/kg × BR × sec/yr × 1/2
    fc_holder_sources = [
        SourceConfig(
            label = "fc_holder_top_Bi214",
            E_MeV = 2.448, entry = :barrel, angular = :flat,
            R_entry = R_TPC, H_entry = 2.5, z_min_entry = 143.1,
            gammas_per_yr = 2.84e4,
        ),
        SourceConfig(
            label = "fc_holder_top_Tl208",
            E_MeV = 2.615, entry = :barrel, angular = :flat,
            R_entry = R_TPC, H_entry = 2.5, z_min_entry = 143.1,
            gammas_per_yr = 3.65e5,
            is_Tl208 = true,
        ),
        SourceConfig(
            label = "fc_holder_bot_Bi214",
            E_MeV = 2.448, entry = :barrel, angular = :flat,
            R_entry = R_TPC, H_entry = 2.5, z_min_entry = -1.25,
            gammas_per_yr = 2.84e4,
            use_extended_volume = true,
        ),
        SourceConfig(
            label = "fc_holder_bot_Tl208",
            E_MeV = 2.615, entry = :barrel, angular = :flat,
            R_entry = R_TPC, H_entry = 2.5, z_min_entry = -1.25,
            gammas_per_yr = 3.65e5,
            is_Tl208 = true,
            use_extended_volume = true,
        ),
    ]

    # ─── PMT top (endcap_top, in gas above gate) ─────────────────
    # Top PMTs are at z = 152.6 cm (in gas, 7 cm above gate).
    # Gas is transparent — no attenuation between PMTs and LXe surface.
    # The tracking volume must extend above z_top (145.6) to include
    # the PMT position, so the gamma can propagate through gas (which
    # we approximate as vacuum — no interactions) down to the LXe surface.
    # Rates from data/lz_pmt_sampling.csv (endcap_top)
    pmt_top_sources = [
        SourceConfig(
            label = "pmt_top_Bi214",
            E_MeV = 2.448, entry = :endcap_top, angular = :flat,
            R_entry = R_TPC, z_entry = Z_PMT_TOP,
            gammas_per_yr = 9.62e4,
            use_extended_volume = true,
        ),
        SourceConfig(
            label = "pmt_top_Tl208",
            E_MeV = 2.615, entry = :endcap_top, angular = :flat,
            R_entry = R_TPC, z_entry = Z_PMT_TOP,
            gammas_per_yr = 1.20e6,
            is_Tl208 = true,
            use_extended_volume = true,
        ),
    ]

    # ─── PMT bottom (endcap_bottom, below RFR in LXe) ────────────
    # Gammas must traverse RFR LXe (invisible) before reaching drift.
    # Source at z = Z_PMT_BOTTOM.
    # Rates from data/lz_pmt_sampling.csv (endcap_bottom)
    pmt_bottom_sources = [
        SourceConfig(
            label = "pmt_bottom_Bi214",
            E_MeV = 2.448, entry = :endcap_bottom, angular = :flat,
            R_entry = R_TPC, z_entry = Z_PMT_BOTTOM,
            gammas_per_yr = 9.63e4,
            use_extended_volume = true,
        ),
        SourceConfig(
            label = "pmt_bottom_Tl208",
            E_MeV = 2.615, entry = :endcap_bottom, angular = :flat,
            R_entry = R_TPC, z_entry = Z_PMT_BOTTOM,
            gammas_per_yr = 1.19e6,
            is_Tl208 = true,
            use_extended_volume = true,
        ),
    ]

    # ─── PMT barrel (cables, side skin — same entry as FC barrel) ─
    # These are at the FC outer wall, same geometry as FC barrel.
    pmt_barrel_sources = [
        SourceConfig(
            label = "pmt_barrel_Bi214",
            E_MeV = 2.448, entry = :barrel, angular = :flat,
            R_entry = R_TPC, H_entry = H_TPC, z_min_entry = Z_MIN_BARREL,
            gammas_per_yr = 1.66e5,
        ),
        SourceConfig(
            label = "pmt_barrel_Tl208",
            E_MeV = 2.615, entry = :barrel, angular = :flat,
            R_entry = R_TPC, H_entry = H_TPC, z_min_entry = Z_MIN_BARREL,
            gammas_per_yr = 9.53e5,
            is_Tl208 = true,
        ),
    ]

    return vcat(cryo_sources, fc_barrel_sources, fc_holder_sources,
                pmt_top_sources, pmt_bottom_sources, pmt_barrel_sources)
end

# =========================================================================
# CLI
# =========================================================================

function parse_cli()
    s = ArgParseSettings(description="XLZD Background MC — Multi-source")
    @add_arg_table! s begin
        "--n-samples"
            arg_type = Int
            default  = 100_000_000
            help     = "Number of gammas per source run"
        "--seed"
            arg_type = Int
            default  = 1234
            help     = "Master RNG seed"
        "--output"
            default  = "last_run"
            help     = "Run name (outputs go to output/<name>/)"
        "--xcom"
            default  = "data/nist.csv"
            help     = "Path to NIST XCOM table"
        "--n-traj"
            arg_type = Int
            default  = 100
            help     = "Trajectories to save per outcome category"
        "--source-index"
            arg_type = Int
            default  = 0
            help     = "Source index (1-8) to run, or 0 for all"
        # ─── FV parameters (CLI) — defaults from LZ 0νββ paper (967 kg) ─
        "--fv-z-min"
            arg_type = Float64
            default  = 26.0
            help     = "FV z_min (cm)"
        "--fv-z-max"
            arg_type = Float64
            default  = 96.0
            help     = "FV z_max (cm)"
        "--fv-r-max"
            arg_type = Float64
            default  = 39.0
            help     = "FV r_max (cm), converted to r²"
        # ─── Cut parameters (CLI) ────────────────────────────────────
        "--roi-halfwidth"
            arg_type = Float64
            default  = 17.2
            help     = "ROI half-width in keV (±1σ at σ/E=0.7%)"
        "--dz-threshold"
            arg_type = Float64
            default  = 3.0
            help     = "SS/MS z-threshold in mm"
        "--e-visible"
            arg_type = Float64
            default  = 5.0
            help     = "Visible deposit threshold in keV"
        "--e-cutoff"
            arg_type = Float64
            default  = 40.0
            help     = "Photon tracking energy cutoff in keV"
        # ─── Companion veto (²⁰⁸Tl) ─────────────────────────────────
        "--companion-veto"
            arg_type = Float64
            default  = 5.0
            help     = "Companion gamma veto threshold in active LXe (keV)"
        "--skin-veto"
            arg_type = Float64
            default  = 100.0
            help     = "Companion gamma veto threshold in LXe skin (keV)"
    end
    parse_args(s)
end

# =========================================================================
# Main
# =========================================================================

function main()
    args = parse_cli()
    sources = build_sources()

    params = Params(
        mc_N_samples                = args["n-samples"],
        mc_seed                     = args["seed"],
        phys_xcom_data_path         = args["xcom"],
        mc_n_traj_per_outcome       = args["n-traj"],
        fv_z_min_cm                 = args["fv-z-min"],
        fv_z_max_cm                 = args["fv-z-max"],
        fv_r2_max_cm2               = args["fv-r-max"]^2,
        cut_ROI_halfwidth_keV       = args["roi-halfwidth"],
        cut_Δz_threshold_mm         = args["dz-threshold"],
        cut_E_visible_threshold_keV = args["e-visible"],
        cut_E_tracking_cutoff_keV   = args["e-cutoff"],
        cut_companion_veto_keV      = args["companion-veto"],
        cut_skin_veto_keV           = args["skin-veto"],
    )

    src_idx = args["source-index"]

    if src_idx == 0
        # Run all sources sequentially
        println("Running all $(length(sources)) sources...\n")
        total_bg = 0.0
        for (i, source) in enumerate(sources)
            @printf("━━━ Source %d / %d: %s ━━━\n", i, length(sources), source.label)
            p = copy_params(params;
                       out_dir = joinpath("output", args["output"], source.label),
                       mc_seed = params.mc_seed + (i - 1) * 100)
            result = run_mc(p, source)
            write_outputs(result)
            n_roi = result.counts[:SS_in_ROI]
            f = n_roi / p.mc_N_samples
            bg = f * source.gammas_per_yr
            total_bg += bg
            @printf("  → f_ROI = %.4e,  bg = %.4e events/yr\n\n", f, bg)
        end
        @printf("\n═══ TOTAL BACKGROUND FROM ALL SOURCES ═══\n")
        @printf("  %.4e events/yr\n", total_bg)
        @printf("═════════════════════════════════════════\n")
    else
        # Run a single source
        if src_idx < 1 || src_idx > length(sources)
            error("Invalid source index $src_idx. Valid range: 1-$(length(sources))")
        end
        source = sources[src_idx]
        @printf("Running source %d: %s\n", src_idx, source.label)
        p = copy_params(params;
                   out_dir = joinpath("output", args["output"], source.label))
        result = run_mc(p, source)
        write_outputs(result)
    end
end

main()
