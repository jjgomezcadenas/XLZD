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
        "--output"
            default  = "last_run_v2"
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
        "--no-rejection-histograms"
            action   = :store_true
            help     = "Skip writing the rejected_skin / rejected_fv histograms."
        "--no-stack-histos"
            action   = :store_true
            help     = "Skip accumulating + writing StackHistogramSet CSVs."
        "--no-cluster-histos"
            action   = :store_true
            help     = "Skip accumulating + writing ClusterHistogramSet CSVs."
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
    outname   = args["output"]
    src_i     = args["source-index"]
    sigma_e   = args["sigma-e"]
    roi_half  = args["roi-halfwidth"]
    fv_z_min  = args["fv-z-min"]
    fv_z_max  = args["fv-z-max"]
    fv_r_max  = args["fv-r-max"]
    do_rej_hist       = !args["no-rejection-histograms"]
    do_stack_hist     = !args["no-stack-histos"]
    do_cluster_hist   = !args["no-cluster-histos"]
    sample_stack_n    = args["sample-stack"]

    # Construct output directory up front so per-source paths can be
    # built before the run_mc calls (needed for sample-stack output).
    out_dir = joinpath(PROJECT_ROOT, "output", outname)
    mkpath(out_dir)

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
                              with_stack_histograms=do_stack_hist,
                              with_cluster_histograms=do_cluster_hist,
                              with_rejection_histograms=do_rej_hist,
                              sample_stack=sample_stack_n,
                              sample_stack_dir=out_dir)
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
                      joinpath(out_dir, "hist_$(eff.name)", "stack_sample.csv") :
                      nothing
        push!(results, run_mc(det, eff, comp_eff, xcom, params, N;
                               mc_seed=seed, verbose=true,
                               with_stack_histograms=do_stack_hist,
                               with_cluster_histograms=do_cluster_hist,
                               with_rejection_histograms=do_rej_hist,
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

    # Per-outcome breakdown — counts and % of n_total per source
    println("── Per-outcome breakdown ──")
    @printf("  %-12s %12s %12s %12s %12s %12s %12s %12s\n",
            "source", "esc", "MS", "skin", "outFV",
            "outROI", "in_ROI", "comp_veto")
    println("  ", "─"^114)
    for r in results
        n  = r.n_total
        pc = (k::Symbol) -> 100.0 * r.counts[k] / n
        cell = (k::Symbol) -> @sprintf("%7d (%4.1f%%)", r.counts[k], pc(k))
        @printf("  %-12s %12s %12s %12s %12s %12s %12s %12s\n",
                r.name,
                cell(:escaped),    cell(:MS_rejected),
                cell(:skin_vetoed), cell(:outside_FV),
                cell(:SS_outside_ROI), cell(:SS_in_ROI),
                cell(:companion_vetoed))
    end
    println()

    # Save CSV
    csv_path = joinpath(out_dir, "summary.csv")
    open(csv_path, "w") do f
        println(f, "source,isotope,n_total,gamma_per_yr_total,f_SS_in_ROI,",
                   "bg_per_yr,r_comp,runtime_s,",
                   "n_escaped,n_MS,n_skin_vetoed,n_SS_outside_FV,",
                   "n_SS_outside_ROI,n_SS_in_ROI,n_companion_vetoed")
        for r in results
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
            ], ","))
        end
    end
    println("  → wrote $csv_path")

    # ───────── Histogram CSVs (one folder per source) ─────────
    if do_rej_hist || do_stack_hist || do_cluster_hist
        for r in results
            h_dir = joinpath(out_dir, "hist_$(r.name)")
            wrote = false
            if do_rej_hist && r.rej_hist !== nothing
                mkpath(h_dir); _save_rejection_csvs(h_dir, r.rej_hist); wrote = true
            end
            if do_stack_hist && r.stack_hists !== nothing
                mkpath(h_dir); _save_stack_csvs(h_dir, r.stack_hists); wrote = true
            end
            if do_cluster_hist && r.cluster_hists !== nothing
                mkpath(h_dir); _save_cluster_csvs(h_dir, r.cluster_hists); wrote = true
            end
            wrote && println("  → wrote $h_dir")
        end
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
    _write_int_bins(joinpath(dir, "stack_ng_max.csv"),
                    sh.ng_max_counts; label="ng_max")
    _write_int_bins(joinpath(dir, "stack_n_photo.csv"),
                    sh.n_photo_counts; label="n_photo")
    _write_int_bins(joinpath(dir, "stack_n_compton.csv"),
                    sh.n_compton_counts; label="n_compton")
    _write_int_bins(joinpath(dir, "stack_n_pair.csv"),
                    sh.n_pair_counts; label="n_pair")
    _write_int_bins(joinpath(dir, "stack_n_below_thresh.csv"),
                    sh.n_below_thresh_counts; label="n_below_thresh")

    # First-interaction type (4 bins, labelled rather than indexed)
    open(joinpath(dir, "stack_first_interaction.csv"), "w") do f
        println(f, "interaction,count")
        for i in 1:4
            println(f, "$(_INTERACTION_LABELS[i]),$(sh.first_interaction_counts[i])")
        end
    end

    _write_1d(joinpath(dir, "stack_inclusive_edep.csv"),
              sh.inclusive_edep_counts;
              lo=0.0, hi=sh.E_max_MeV,
              left_label="bin_left_MeV", right_label="bin_right_MeV")
    _write_1d(joinpath(dir, "stack_E_first.csv"),
              sh.E_first_counts;
              lo=0.0, hi=sh.E_max_MeV,
              left_label="bin_left_MeV", right_label="bin_right_MeV")
    _write_1d(joinpath(dir, "stack_delta_z.csv"),
              sh.Δz_counts;
              lo=0.0, hi=sh.Δz_max_cm,
              left_label="bin_left_cm", right_label="bin_right_cm")

    # Region × interaction (3×4 matrix)
    open(joinpath(dir, "stack_region_interaction.csv"), "w") do f
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
    # 1D energies
    for (name, V) in (("cluster_Ec.csv",  ch.Ec_counts),
                       ("cluster_Emax.csv", ch.Emax_counts),
                       ("cluster_Emin.csv", ch.Emin_counts),
                       ("cluster_Einc.csv", ch.Einc_counts))
        _write_1d(joinpath(dir, name), V;
                  lo=0.0, hi=ch.E_max_MeV,
                  left_label="bin_left_MeV", right_label="bin_right_MeV")
    end

    # 1D pair distances
    _write_1d(joinpath(dir, "cluster_closest_D3.csv"),  ch.closest_D3_counts;
              lo=0.0, hi=ch.D3_max_cm,
              left_label="bin_left_cm", right_label="bin_right_cm")
    _write_1d(joinpath(dir, "cluster_furthest_D3.csv"), ch.furthest_D3_counts;
              lo=0.0, hi=ch.D3_max_cm,
              left_label="bin_left_cm", right_label="bin_right_cm")
    _write_1d(joinpath(dir, "cluster_closest_dz.csv"),  ch.closest_dz_counts;
              lo=0.0, hi=ch.dz_max_cm,
              left_label="bin_left_cm", right_label="bin_right_cm")
    _write_1d(joinpath(dir, "cluster_furthest_dz.csv"), ch.furthest_dz_counts;
              lo=0.0, hi=ch.dz_max_cm,
              left_label="bin_left_cm", right_label="bin_right_cm")

    # Cluster multiplicity (integer bucket)
    _write_int_bins(joinpath(dir, "cluster_N_clusters.csv"),
                    ch.N_clusters_counts; label="n_clusters")

    # 2D heatmaps
    _write_2d(joinpath(dir, "cluster_r2_vs_z.csv"), ch.r2_vs_z_2d_counts;
              x_lo=0.0,           x_hi=ch.r2_max_cm2, x_n=ch.r2_n_bins,
              x_left="bin_r2_left_cm2",  x_right="bin_r2_right_cm2",
              y_lo=ch.z_min_cm,   y_hi=ch.z_max_cm,   y_n=ch.z_n_bins,
              y_left="bin_z_left_cm",    y_right="bin_z_right_cm")
    _write_2d(joinpath(dir, "cluster_D_vs_z.csv"), ch.D_vs_z_2d_counts;
              x_lo=0.0,           x_hi=ch.D_max_cm,   x_n=ch.D_n_bins,
              x_left="bin_D_left_cm",    x_right="bin_D_right_cm",
              y_lo=ch.z_min_cm,   y_hi=ch.z_max_cm,   y_n=ch.z_n_bins,
              y_left="bin_z_left_cm",    y_right="bin_z_right_cm")
end

function _save_rejection_csvs(dir::String, rh::RejectionHistograms)
    r2_bw = rh.r2_max_cm2 / rh.r2_n_bins
    z_bw  = (rh.z_max_cm - rh.z_min_cm) / rh.z_n_bins
    e_bw  = rh.E_max_MeV / rh.E_n_bins
    function _write_r2z(path, M)
        open(path, "w") do f
            println(f, "bin_r2_left_cm2,bin_r2_right_cm2,bin_z_left_cm,bin_z_right_cm,count")
            for i in 1:rh.r2_n_bins, j in 1:rh.z_n_bins
                r2_l = (i-1) * r2_bw
                r2_r = i      * r2_bw
                z_l  = rh.z_min_cm + (j-1) * z_bw
                z_r  = rh.z_min_cm + j     * z_bw
                println(f, "$r2_l,$r2_r,$z_l,$z_r,$(M[i, j])")
            end
        end
    end
    function _write_E(path, V)
        open(path, "w") do f
            println(f, "bin_left_MeV,bin_right_MeV,count")
            for i in 1:rh.E_n_bins
                println(f, "$((i-1)*e_bw),$(i*e_bw),$(V[i])")
            end
        end
    end
    _write_r2z(joinpath(dir, "rejected_skin_r2z.csv"), rh.skin_r2z_counts)
    _write_E(  joinpath(dir, "rejected_skin_E.csv"),   rh.skin_E_counts)
    _write_r2z(joinpath(dir, "rejected_fv_r2z.csv"),   rh.fv_r2z_counts)
    _write_E(  joinpath(dir, "rejected_fv_E.csv"),     rh.fv_E_counts)
end

main()
