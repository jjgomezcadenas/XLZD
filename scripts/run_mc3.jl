# scripts/run_mc3.jl — End-to-end driver for the src3/ Monte Carlo.
#
# Runs the per-photon MC over all 6 main effective sources (Bi-214 and
# Tl-208 main × 3 regions) with the companion veto applied for the
# Tl-208 sources. Prints a summary table to stdout and writes a CSV.
#
# Run from project root:
#     julia --project=. -t 8 scripts/run_mc3.jl              # 1e6 default
#     julia --project=. -t 8 scripts/run_mc3.jl --n-samples 10000000
#     julia --project=. -t 8 scripts/run_mc3.jl --source-index 1

using ArgParse
using Printf

include(joinpath(@__DIR__, "..", "src3", "XLZD3.jl"))
using .XLZD3

const PROJECT_ROOT = abspath(joinpath(@__DIR__, ".."))
const TI_PATH      = joinpath(PROJECT_ROOT, "data", "nist_ti.csv")
const LXE_NIST     = joinpath(PROJECT_ROOT, "data", "nist_lxe.csv")
const XCOM_PATH    = joinpath(PROJECT_ROOT, "data", "nist.csv")
const LXE_CSV      = joinpath(PROJECT_ROOT, "data", "lxe_detector.csv")
const GEOM_CSV     = joinpath(PROJECT_ROOT, "data", "lz_cryo_geometry.csv")
const EXTRAS_CSV   = joinpath(PROJECT_ROOT, "data", "lz_cryo_extras.csv")
const SURFACES_CSV = joinpath(PROJECT_ROOT, "data", "lz_cryo_surface_sources.csv")

function parse_cli()
    s = ArgParseSettings(description="src3/ MC driver — cryostat-only")
    @add_arg_table! s begin
        "--n-samples"
            arg_type = Int
            default  = 1_000_000
        "--seed"
            arg_type = Int
            default  = 1234
        "--run-name"
            default  = "last_run_v2"
            help     = "Descriptive run name. Each source-job writes to " *
                       "output/<run-name>/<source>/. The Python summarizer " *
                       "globs **/summary.csv under output/<run-name>/."
        "--source-index"
            arg_type = Int
            default  = 0           # 0 = all 6 main sources
        "--sigma-e"
            arg_type = Float64
            default  = 0.007       # σ_E / E (LZ default 0.7 %)
            help     = "Energy resolution σ/E. bb0nu uses 0.010."
        "--roi-halfwidth"
            arg_type = Float64
            default  = 17.2        # keV; ±1σ at σ/E = 0.7 % and Q = 2458 keV
            help     = "ROI half-width in keV (±1σ around Q_ββ)."
        "--fv-z-min"
            arg_type = Float64
            default  = 26.0
            help     = "FV z minimum (cm). Default = bb0nu inner volume."
        "--fv-z-max"
            arg_type = Float64
            default  = 96.0
            help     = "FV z maximum (cm)."
        "--fv-r-max"
            arg_type = Float64
            default  = 39.0
            help     = "FV r maximum (cm). r² = (fv-r-max)² is what MCParams stores."
        "--histos"
            arg_type = Bool
            default  = true
            help     = "Write all histogram CSVs (cuts + diagnostics). " *
                       "Pass `--histos false` to skip them entirely."
        "--sample-stack"
            arg_type = Int
            default  = 0
            help     = "Stack-dump CSV: 0=off, -1=all events, n>0 ≈ uniform sample of n events per source."
    end
    parse_args(s)
end

function main()
    args = parse_cli()
    N         = args["n-samples"]
    seed      = args["seed"]
    run_name  = args["run-name"]
    src_i     = args["source-index"]
    sigma_e   = args["sigma-e"]
    roi_half  = args["roi-halfwidth"]
    fv_z_min  = args["fv-z-min"]
    fv_z_max  = args["fv-z-max"]
    fv_r_max  = args["fv-r-max"]
    do_histos         = args["histos"]
    sample_stack_n    = args["sample-stack"]

    # Construct output directory up front so per-source paths can be
    # built before the run_mc calls (needed for sample-stack output).
    run_dir = joinpath(PROJECT_ROOT, "output", run_name)
    mkpath(run_dir)

    println("\n── src3/ MC driver — Cryostat backgrounds ──")
    @printf("  N samples per source : %d\n", N)
    @printf("  Master seed          : %d\n", seed)
    @printf("  Threads              : %d\n", Threads.nthreads())
    @printf("  σ_E / E              : %.4f\n", sigma_e)
    @printf("  ROI half-width       : %.2f keV\n", roi_half)
    @printf("  FV z range           : [%.1f, %.1f] cm\n", fv_z_min, fv_z_max)
    @printf("  FV r max             : %.1f cm\n", fv_r_max)

    # Load
    mat_LXe = load_material("LXe", 2.953, LXE_NIST)
    mat_Ti  = load_material("Ti",  4.510, TI_PATH)
    det     = build_lxe_detector(LXE_CSV, mat_LXe)
    cryo    = build_cryostat(GEOM_CSV, EXTRAS_CSV, SURFACES_CSV)
    indiv   = build_individual_sources(cryo, mat_Ti)
    effs    = build_effective_sources(indiv, cryo, mat_Ti)
    xcom    = load_xcom(XCOM_PATH)
    params  = MCParams(σ_E_over_E       = sigma_e,
                       ROI_halfwidth_keV = roi_half,
                       fv_z_min_cm       = fv_z_min,
                       fv_z_max_cm       = fv_z_max,
                       fv_r2_max_cm2     = fv_r_max^2)

    # Run
    main_names = ["CB_Bi214", "CTH_Bi214", "CBH_Bi214",
                  "CB_Tl208", "CTH_Tl208", "CBH_Tl208"]
    by_name = Dict(e.name => e for e in effs)
    results = MCResult[]
    if src_i == 0
        results = run_mc_all(det, effs, xcom, params, N;
                              mc_seed=seed, verbose=true,
                              with_stack_histograms=do_histos,
                              with_cluster_histograms=do_histos,
                              with_cut_histograms=do_histos,
                              sample_stack=sample_stack_n,
                              sample_stack_dir=run_dir)
    else
        1 <= src_i <= length(main_names) ||
            error("source-index must be 1..$(length(main_names))")
        eff = by_name[main_names[src_i]]
        comp_eff = if eff.isotope === :Tl208
            get(by_name, replace(eff.name, "_Tl208" => "_Tl208c"), nothing)
        else
            nothing
        end
        @printf("\n── Running %s ──\n", eff.name)
        sample_path = sample_stack_n != 0 ?
                      joinpath(run_dir, eff.name, "stack_sample.csv") :
                      nothing
        push!(results, run_mc(det, eff, comp_eff, xcom, params, N;
                               mc_seed=seed, verbose=true,
                               with_stack_histograms=do_histos,
                               with_cluster_histograms=do_histos,
                               with_cut_histograms=do_histos,
                               sample_stack=sample_stack_n,
                               sample_stack_path=sample_path))
    end

    # Print summary
    println("\n── Per-source results ──")
    @printf("  %-12s %-7s %12s %10s %12s %12s\n",
            "source", "iso", "γ/yr_in", "f_SS_ROI", "bg γ/yr", "runtime (s)")
    println("  ", "─"^77)
    total_bg = 0.0
    for r in results
        @printf("  %-12s %-7s %12.3e %10.3e %12.3e %12.2f\n",
                r.name, string(r.isotope),
                r.γ_per_yr_total, r.f_SS_in_ROI, r.bg_per_yr, r.runtime_s)
        total_bg += r.bg_per_yr
    end
    println("  ", "─"^77)
    @printf("  TOTAL background (cryostat sources) : %.4e events/yr\n", total_bg)
    println()

    # Run-wide analysis cuts (same for all sources). Printed together so the
    # numbers above can be read against the configuration that produced them.
    let
        roi_lo_keV = params.Q_betabeta_keV - params.ROI_halfwidth_keV
        roi_hi_keV = params.Q_betabeta_keV + params.ROI_halfwidth_keV
        println("── Analysis cuts ──")
        @printf("  Q_ββ                : %.1f keV\n", params.Q_betabeta_keV)
        @printf("  σ_E / E             : %.4f\n",     params.σ_E_over_E)
        @printf("  ROI window          : Q_ββ ± %.2f keV  →  [%.2f, %.2f] keV\n",
                params.ROI_halfwidth_keV, roi_lo_keV, roi_hi_keV)
        @printf("  FV box              : z ∈ [%.1f, %.1f] cm,  r ≤ %.1f cm\n",
                params.fv_z_min_cm, params.fv_z_max_cm,
                sqrt(params.fv_r2_max_cm2))
        @printf("  Visible threshold   : %.1f keV  (per-cluster)\n", params.E_visible_keV)
        @printf("  Skin-veto threshold : %.1f keV  (cumulative on :Skin)\n", params.E_skin_veto_keV)
        @printf("  Tracking cutoff     : %.1f keV  (residual γ energy)\n", params.E_tracking_cutoff_keV)
        println()
    end

    # Per-source funnel — cumulative survivors at each selection stage,
    # so the killer cut is visible at a glance. Each stage shows:
    #   N_i        absolute count surviving up to this stage
    #   cum=...    N_i / N_total            (cumulative survival)
    #   acc=...    N_i / N_{i-1}            (per-stage acceptance)
    println("── Per-source funnel ──")
    for r in results
        n0 = r.n_total
        n1 = n0 - r.counts[:escaped]                       # into detector
        n2 = n1 - r.counts[:skin_vetoed]                   # pass skin veto
        n3 = n2 - r.counts[:outside_FV]                    # inside FV
        n4 = n3 - r.counts[:MS_rejected]                   # accepted as SS
        n5 = n4 - r.counts[:SS_outside_ROI]                # in ROI (pre-companion)
        n6 = n5 - r.counts[:companion_vetoed]              # after companion veto

        @printf("\n  ── %s ──\n", r.name)
        @printf("    %-24s N  = %12d  cum=%.3e\n",
                "Total events run", n0, 1.0)
        rows = (
            ("Into detector",         n1),
            ("Pass skin veto",        n2),
            ("Inside FV",             n3),
            ("Accepted as SS",        n4),
            ("In ROI (pre-companion)", n5),
            ("After companion veto",  n6),
        )
        prev = n0
        for (i, (label, ni)) in enumerate(rows)
            cum = ni / n0
            acc_str = prev == 0 ? "  ---  " : @sprintf("%.3e", ni / prev)
            @printf("    %-24s N%-1d = %12d  cum=%.3e  acc=%s\n",
                    label, i, ni, cum, acc_str)
            prev = ni
        end
    end
    println()

    # Per-source output: one directory per source under run_dir, holding
    # a single-row summary.csv plus all histogram CSVs. Multiple sources
    # in one process (the all-sources mode) write multiple dirs; the
    # Python summarizer aggregates by globbing **/summary.csv under run_dir.
    fv_r_max_cm = sqrt(params.fv_r2_max_cm2)
    for r in results
        src_dir = joinpath(run_dir, r.name)
        mkpath(src_dir)

        csv_path = joinpath(src_dir, "summary.csv")
        open(csv_path, "w") do f
            println(f, "source,isotope,n_total,gamma_per_yr_total,f_SS_in_ROI,",
                       "bg_per_yr,r_comp,runtime_s,",
                       "n_escaped,n_MS,n_skin_vetoed,n_SS_outside_FV,",
                       "n_SS_outside_ROI,n_SS_in_ROI,n_companion_vetoed,",
                       "Q_betabeta_keV,sigma_E_over_E,ROI_halfwidth_keV,",
                       "fv_z_min_cm,fv_z_max_cm,fv_r_max_cm,",
                       "E_visible_keV,E_skin_veto_keV,",
                       "R_LXe_outer_cm")
            println(f, join([
                r.name, r.isotope, r.n_total,
                @sprintf("%.6e", r.γ_per_yr_total),
                @sprintf("%.6e", r.f_SS_in_ROI),
                @sprintf("%.6e", r.bg_per_yr),
                @sprintf("%.6e", r.r_comp),
                @sprintf("%.3f",  r.runtime_s),
                r.counts[:escaped], r.counts[:MS_rejected],
                r.counts[:skin_vetoed], r.counts[:outside_FV],
                r.counts[:SS_outside_ROI], r.counts[:SS_in_ROI],
                r.counts[:companion_vetoed],
                @sprintf("%.3f", params.Q_betabeta_keV),
                @sprintf("%.6f", params.σ_E_over_E),
                @sprintf("%.3f", params.ROI_halfwidth_keV),
                @sprintf("%.3f", params.fv_z_min_cm),
                @sprintf("%.3f", params.fv_z_max_cm),
                @sprintf("%.3f", fv_r_max_cm),
                @sprintf("%.3f", params.E_visible_keV),
                @sprintf("%.3f", params.E_skin_veto_keV),
                @sprintf("%.3f", det.R_ICV_inner),
            ], ","))
        end

        if do_histos
            r.stack_hists   !== nothing && _save_stack_csvs(src_dir,   r.stack_hists)
            r.cluster_hists !== nothing && _save_cluster_csvs(src_dir, r.cluster_hists)
            r.cut_hists     !== nothing && _save_cut_csvs(src_dir,     r.cut_hists)
        end

        println("  → wrote $src_dir")
    end
    println()
end

# ---------------------------------------------------------------------------
# Generic CSV writers (shared by stack and cluster sets)
# ---------------------------------------------------------------------------

"Write integer-bucket histogram (bins 0..length(V)-1) to a 2-column CSV."
function _write_int_bins(path::String, V::Vector{Int}; label::String="n")
    open(path, "w") do f
        println(f, "$label,count")
        for i in 1:length(V)
            println(f, "$(i-1),$(V[i])")
        end
    end
end

"Write a uniform 1-D histogram to a 3-column CSV with bin edges + count."
function _write_1d(path::String, V::Vector{Int};
                    lo::Float64, hi::Float64,
                    left_label::String, right_label::String)
    n = length(V)
    bw = (hi - lo) / n
    open(path, "w") do f
        println(f, "$left_label,$right_label,count")
        for i in 1:n
            println(f, "$(lo + (i-1)*bw),$(lo + i*bw),$(V[i])")
        end
    end
end

"Write a uniform 2-D histogram (i along x, j along y) to a 5-column CSV."
function _write_2d(path::String, M::Matrix{Int};
                    x_lo::Float64, x_hi::Float64, x_n::Int, x_left::String, x_right::String,
                    y_lo::Float64, y_hi::Float64, y_n::Int, y_left::String, y_right::String)
    x_bw = (x_hi - x_lo) / x_n
    y_bw = (y_hi - y_lo) / y_n
    open(path, "w") do f
        println(f, "$x_left,$x_right,$y_left,$y_right,count")
        for i in 1:x_n, j in 1:y_n
            xl = x_lo + (i-1)*x_bw; xr = x_lo + i*x_bw
            yl = y_lo + (j-1)*y_bw; yr = y_lo + j*y_bw
            println(f, "$xl,$xr,$yl,$yr,$(M[i, j])")
        end
    end
end

# ---------------------------------------------------------------------------
# StackHistogramSet → CSVs
# ---------------------------------------------------------------------------

const _INTERACTION_LABELS = ("PHOTO", "COMPTON", "PAIR", "BELOW_THRESH")
const _STACK_REGION_LABELS = ("TPC", "Skin", "Inert")

function _save_stack_csvs(dir::String, sh::StackHistogramSet)
    # Interaction-type frequency (one entry per stack row, normalisable to
    # the cross-section ratios at the source energy).
    open(joinpath(dir, "diag_interaction_type_freq.csv"), "w") do f
        println(f, "interaction,count")
        for i in 1:length(_INTERACTION_LABELS)
            println(f, "$(_INTERACTION_LABELS[i]),$(sh.interaction_type_freq[i])")
        end
    end

    # Per-photon LXe path length (1D bar).
    _write_1d(joinpath(dir, "diag_path_length_LXe.csv"),
              sh.path_length_LXe_counts;
              lo=0.0, hi=sh.path_length_max_cm,
              left_label="bin_left_cm", right_label="bin_right_cm")

    # Region × interaction matrix (3×4).
    open(joinpath(dir, "diag_region_interaction.csv"), "w") do f
        println(f, "region,interaction,count")
        for ri in 1:3, ci in 1:4
            println(f, "$(_STACK_REGION_LABELS[ri]),$(_INTERACTION_LABELS[ci]),$(sh.region_interaction_counts[ri, ci])")
        end
    end
end

# ---------------------------------------------------------------------------
# ClusterHistogramSet → CSVs
# ---------------------------------------------------------------------------

function _save_cluster_csvs(dir::String, ch::ClusterHistogramSet)
    _write_1d(joinpath(dir, "diag_cluster_Ec.csv"), ch.Ec_counts;
              lo=0.0, hi=ch.E_max_MeV,
              left_label="bin_left_MeV", right_label="bin_right_MeV")
end

# ---------------------------------------------------------------------------
# CutHistograms → CSVs (cut-flow set, one CSV per histogram)
# ---------------------------------------------------------------------------

function _save_cut_csvs(dir::String, ch::CutHistograms)
    # Cut 1
    _write_1d(joinpath(dir, "cut1_h_u_sampled.csv"),
              ch.h_u_sampled;
              lo=0.0, hi=1.0,
              left_label="bin_u_left", right_label="bin_u_right")
    # Cut 2
    _write_2d(joinpath(dir, "cut2_first_interaction_r_z.csv"),
              ch.first_interaction_r_z;
              x_lo=0.0,           x_hi=ch.r_max_cm,   x_n=ch.r_n_bins,
              x_left="bin_r_left_cm",   x_right="bin_r_right_cm",
              y_lo=ch.z_min_cm,   y_hi=ch.z_max_cm,   y_n=ch.z_n_bins,
              y_left="bin_z_left_cm",   y_right="bin_z_right_cm")
    # Cut 3
    _write_1d(joinpath(dir, "cut3_dz_inclusive.csv"),
              ch.dz_inclusive;
              lo=0.0, hi=ch.dz_max_cm,
              left_label="bin_left_cm", right_label="bin_right_cm")
    _write_int_bins(joinpath(dir, "cut3_n_visible.csv"),
                    ch.n_visible; label="n_visible")
    _write_1d(joinpath(dir, "cut3_E_total.csv"),
              ch.E_total;
              lo=0.0, hi=ch.E_max_MeV,
              left_label="bin_left_MeV", right_label="bin_right_MeV")
    # Cut 4
    _write_1d(joinpath(dir, "cut4_ss_ec_pre_roi.csv"),
              ch.ss_ec_pre_roi;
              lo=0.0, hi=ch.E_max_MeV,
              left_label="bin_left_MeV", right_label="bin_right_MeV")
    _write_1d(joinpath(dir, "cut4_ss_es_pre_roi.csv"),
              ch.ss_es_pre_roi;
              lo=0.0, hi=ch.E_max_MeV,
              left_label="bin_left_MeV", right_label="bin_right_MeV")
    _write_2d(joinpath(dir, "cut4_ss_r_z.csv"),
              ch.ss_r_z;
              x_lo=0.0,           x_hi=ch.r_max_cm,   x_n=ch.r_n_bins,
              x_left="bin_r_left_cm",   x_right="bin_r_right_cm",
              y_lo=ch.z_min_cm,   y_hi=ch.z_max_cm,   y_n=ch.z_n_bins,
              y_left="bin_z_left_cm",   y_right="bin_z_right_cm")
end

main()
