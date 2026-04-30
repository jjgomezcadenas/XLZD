# scripts/run_mc2.jl — End-to-end driver for the src2/ Monte Carlo.
#
# Runs the per-photon MC over all 6 main effective sources (Bi-214 and
# Tl-208 main × 3 regions) with the companion veto applied for the
# Tl-208 sources. Prints a summary table to stdout and writes a CSV.
#
# Run from project root:
#     julia --project=. -t 8 scripts/run_mc2.jl              # 1e6 default
#     julia --project=. -t 8 scripts/run_mc2.jl --n-samples 10000000
#     julia --project=. -t 8 scripts/run_mc2.jl --source-index 1

using ArgParse
using Printf

include(joinpath(@__DIR__, "..", "src2", "XLZD2.jl"))
using .XLZD2

const PROJECT_ROOT = abspath(joinpath(@__DIR__, ".."))
const TI_PATH      = joinpath(PROJECT_ROOT, "data", "nist_ti.csv")
const LXE_NIST     = joinpath(PROJECT_ROOT, "data", "nist_lxe.csv")
const XCOM_PATH    = joinpath(PROJECT_ROOT, "data", "nist.csv")
const LXE_CSV      = joinpath(PROJECT_ROOT, "data", "lxe_detector.csv")
const GEOM_CSV     = joinpath(PROJECT_ROOT, "data", "lz_cryo_geometry.csv")
const EXTRAS_CSV   = joinpath(PROJECT_ROOT, "data", "lz_cryo_extras.csv")
const SURFACES_CSV = joinpath(PROJECT_ROOT, "data", "lz_cryo_surface_sources.csv")

function parse_cli()
    s = ArgParseSettings(description="src2/ MC driver — cryostat-only")
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
        "--histograms"
            action   = :store_true
            help     = "Accumulate and save the 6 control histograms."
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
    do_hist   = args["histograms"]

    println("\n── src2/ MC driver — Cryostat backgrounds ──")
    @printf("  N samples per source : %d\n", N)
    @printf("  Master seed          : %d\n", seed)
    @printf("  Threads              : %d\n", Threads.nthreads())
    @printf("  σ_E / E              : %.4f\n", sigma_e)
    @printf("  ROI half-width       : %.2f keV\n", roi_half)

    # Load
    mat_LXe = load_material("LXe", 2.953, LXE_NIST)
    mat_Ti  = load_material("Ti",  4.510, TI_PATH)
    det     = build_lxe_detector(LXE_CSV, mat_LXe)
    cryo    = build_cryostat(GEOM_CSV, EXTRAS_CSV, SURFACES_CSV)
    indiv   = build_individual_sources(cryo, mat_Ti)
    effs    = build_effective_sources(indiv, cryo, mat_Ti)
    xcom    = load_xcom(XCOM_PATH)
    params  = MCParams(σ_E_over_E       = sigma_e,
                       ROI_halfwidth_keV = roi_half)

    # Run
    main_names = ["CB_Bi214", "CTH_Bi214", "CBH_Bi214",
                  "CB_Tl208", "CTH_Tl208", "CBH_Tl208"]
    by_name = Dict(e.name => e for e in effs)
    results = MCResult[]
    if src_i == 0
        results = run_mc_all(det, effs, xcom, params, N;
                              mc_seed=seed, verbose=true,
                              with_histograms=do_hist)
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
        push!(results, run_mc(det, eff, comp_eff, xcom, params, N;
                               mc_seed=seed, verbose=true,
                               with_histograms=do_hist))
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

    # Per-outcome breakdown
    println("── Per-outcome breakdown ──")
    @printf("  %-12s %10s %10s %10s %10s %10s %10s %10s\n",
            "source", "esc", "MS", "skin", "outFV",
            "outROI", "in_ROI", "comp_veto")
    println("  ", "─"^96)
    for r in results
        @printf("  %-12s %10d %10d %10d %10d %10d %10d %10d\n",
                r.name, r.counts[:escaped], r.counts[:MS_rejected],
                r.counts[:skin_vetoed], r.counts[:SS_outside_FV],
                r.counts[:SS_outside_ROI], r.counts[:SS_in_ROI],
                r.counts[:companion_vetoed])
    end
    println()

    # Save CSV
    out_dir = joinpath(PROJECT_ROOT, "output", outname)
    mkpath(out_dir)
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
                r.counts[:skin_vetoed], r.counts[:SS_outside_FV],
                r.counts[:SS_outside_ROI], r.counts[:SS_in_ROI],
                r.counts[:companion_vetoed],
            ], ","))
        end
    end
    println("  → wrote $csv_path")

    # ───────── Histogram CSVs (one folder per source) ─────────
    if do_hist
        for r in results
            r.histograms === nothing && continue
            h_dir = joinpath(out_dir, "hist_$(r.name)")
            mkpath(h_dir)
            _save_histogram_csvs(h_dir, r.histograms)
            println("  → wrote $h_dir")
        end
    end
    println()
end

function _save_histogram_csvs(dir::String, h::HistogramSet)
    open(joinpath(dir, "ssms_counts.csv"), "w") do f
        println(f, "bucket,count")
        println(f, "SS,$(h.ssms_counts[1])")
        println(f, "MS,$(h.ssms_counts[2])")
        println(f, "no_cluster,$(h.ssms_counts[3])")
    end
    open(joinpath(dir, "delta_z.csv"), "w") do f
        println(f, "bin_left_cm,bin_right_cm,count")
        bw = h.Δz_max_cm / h.Δz_n_bins
        for i in 1:h.Δz_n_bins
            println(f, "$((i-1)*bw),$(i*bw),$(h.Δz_counts[i])")
        end
    end
    open(joinpath(dir, "e_first.csv"), "w") do f
        println(f, "bin_left_MeV,bin_right_MeV,count")
        bw = h.E_max_MeV / h.E_n_bins
        for i in 1:h.E_n_bins
            println(f, "$((i-1)*bw),$(i*bw),$(h.E_first_counts[i])")
        end
    end
    open(joinpath(dir, "e_cluster.csv"), "w") do f
        println(f, "bin_left_MeV,bin_right_MeV,count")
        bw = h.E_max_MeV / h.E_n_bins
        for i in 1:h.E_n_bins
            println(f, "$((i-1)*bw),$(i*bw),$(h.E_cluster_counts[i])")
        end
    end
    open(joinpath(dir, "n_clusters.csv"), "w") do f
        println(f, "n,count")
        for i in 1:length(h.N_clusters_counts)
            println(f, "$(i-1),$(h.N_clusters_counts[i])")
        end
    end
    open(joinpath(dir, "n_extra_clusters.csv"), "w") do f
        println(f, "n,count")
        for i in 1:length(h.N_extra_counts)
            println(f, "$(i-1),$(h.N_extra_counts[i])")
        end
    end
end

main()
