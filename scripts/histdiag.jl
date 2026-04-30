"""
Histogram diagnostic: detailed event-level distributions for Bi-214 or Tl-208.

Runs N events with NO BFV cut. Produces 5 histogram panels:
  1. SS vs MS bar chart
  2. SS n_interactions distribution
  3. SS (n>1) inclusive deposit energies + |Δz| from cluster center
  4. MS inclusive deposit energies + |Δz| from first deposit
  5. SS energy spectrum (true + smeared) with ROI band

Also prints summary fractions: SS/total, SS_in_ROI/total, SS/MS vs cross sections.

Usage:
  julia --project=. scripts/histdiag.jl --isotope Bi214 --n-events 500000
  julia --project=. scripts/histdiag.jl --isotope Tl208 --n-events 2000000
  julia --project=. scripts/histdiag.jl --isotope Bi214 --n-events 500000 --output output/histdiag
"""

using XLZD, Random, Printf, Plots, Statistics
gr()

# ─── Deposit record ──────────────────────────────────────────────────

struct Deposit
    x::Float64
    y::Float64
    z::Float64
    E_keV::Float64      # energy deposited (keV)
    interaction::Symbol  # :photo, :compton, :pair, :cutoff
end

# ─── Event record ────────────────────────────────────────────────────

struct EventRecord
    deposits::Vector{Deposit}
    outcome::Symbol          # :SS_in_ROI, :SS_outside_ROI, :MS_rejected, :escaped
    E_cluster_true::Float64  # keV, true cluster energy (SS only)
    E_cluster_smeared::Float64  # keV, smeared (SS only)
    n_interactions::Int
end

# ─── Verbose tracker that records all deposits ───────────────────────

function track_with_deposits(rng, source, geom, xcom, params)
    E_visible = params.cut_E_visible_threshold_keV / 1000.0
    E_cutoff  = params.cut_E_tracking_cutoff_keV / 1000.0
    Δz_thresh = params.cut_Δz_threshold_mm / 10.0
    Q_bb      = params.phys_Q_betabeta_keV
    ROI_hw    = params.cut_ROI_halfwidth_keV

    x, y, z = sample_entry_point(rng, source, geom.L_lxe)
    dx, dy, dz = sample_entry_direction(rng, source, x, y, z, geom.L_lxe)
    E = source.E_MeV

    deposits = Deposit[]
    E_cluster = 0.0
    z_cluster = NaN
    cluster_started = false
    n_int = 0
    outcome = :in_progress

    while outcome == :in_progress
        μ_t = μ_total_lin(xcom, E, geom.ρ_LXe)
        d_int = -log(rand(rng)) / μ_t
        d_boundary = path_to_cylinder_exit(x, y, z, dx, dy, dz,
                                           geom.R_lxe, geom.L_lxe)

        if d_int >= d_boundary
            if cluster_started
                σ_E = params.phys_σ_E_over_E * E_cluster
                E_smeared_keV = (E_cluster + σ_E * randn(rng)) * 1000.0
                if abs(E_smeared_keV - Q_bb) ≤ ROI_hw
                    outcome = :SS_in_ROI
                else
                    outcome = :SS_outside_ROI
                end
                return EventRecord(deposits, outcome, E_cluster * 1000.0,
                                   E_smeared_keV, n_int)
            else
                return EventRecord(deposits, :escaped, 0.0, 0.0, n_int)
            end
        end

        x += d_int * dx; y += d_int * dy; z += d_int * dz
        n_int += 1

        σ_pe = σ_photo(xcom, E)
        σ_co = σ_Compton(xcom, E)
        σ_pp = σ_pair(xcom, E)
        σ_tot = σ_pe + σ_co + σ_pp
        r_branch = rand(rng)

        if r_branch < σ_pp / σ_tot
            push!(deposits, Deposit(x, y, z, E * 1000.0, :pair))
            return EventRecord(deposits, :MS_rejected, 0.0, 0.0, n_int)

        elseif r_branch < (σ_pp + σ_pe) / σ_tot
            E_dep = E
            push!(deposits, Deposit(x, y, z, E_dep * 1000.0, :photo))

            if E_dep >= E_visible
                if !cluster_started
                    z_cluster = z; E_cluster = E_dep; cluster_started = true
                else
                    if abs(z - z_cluster) < Δz_thresh
                        E_cluster += E_dep
                    else
                        return EventRecord(deposits, :MS_rejected, 0.0, 0.0, n_int)
                    end
                end
            end

            if outcome == :in_progress
                if cluster_started
                    σ_E = params.phys_σ_E_over_E * E_cluster
                    E_smeared_keV = (E_cluster + σ_E * randn(rng)) * 1000.0
                    if abs(E_smeared_keV - Q_bb) ≤ ROI_hw
                        outcome = :SS_in_ROI
                    else
                        outcome = :SS_outside_ROI
                    end
                    return EventRecord(deposits, outcome, E_cluster * 1000.0,
                                       E_smeared_keV, n_int)
                else
                    return EventRecord(deposits, :escaped, 0.0, 0.0, n_int)
                end
            end

        else  # Compton
            E_scatt, cos_θ = sample_klein_nishina(rng, E)
            E_dep = E - E_scatt
            push!(deposits, Deposit(x, y, z, E_dep * 1000.0, :compton))

            if E_dep >= E_visible
                if !cluster_started
                    z_cluster = z; E_cluster = E_dep; cluster_started = true
                else
                    if abs(z - z_cluster) < Δz_thresh
                        E_cluster += E_dep
                    else
                        return EventRecord(deposits, :MS_rejected, 0.0, 0.0, n_int)
                    end
                end
            end

            φ_az = 2π * rand(rng)
            dx, dy, dz = rotate_direction(dx, dy, dz, cos_θ, φ_az)
            E = E_scatt

            if E < E_cutoff
                push!(deposits, Deposit(x, y, z, E * 1000.0, :cutoff))
                if E >= E_visible
                    if !cluster_started
                        z_cluster = z; E_cluster += E; cluster_started = true
                    else
                        if abs(z - z_cluster) < Δz_thresh
                            E_cluster += E
                        else
                            return EventRecord(deposits, :MS_rejected, 0.0, 0.0, n_int)
                        end
                    end
                end
                if cluster_started
                    σ_E = params.phys_σ_E_over_E * E_cluster
                    E_smeared_keV = (E_cluster + σ_E * randn(rng)) * 1000.0
                    if abs(E_smeared_keV - Q_bb) ≤ ROI_hw
                        outcome = :SS_in_ROI
                    else
                        outcome = :SS_outside_ROI
                    end
                    return EventRecord(deposits, outcome, E_cluster * 1000.0,
                                       E_smeared_keV, n_int)
                else
                    return EventRecord(deposits, :escaped, 0.0, 0.0, n_int)
                end
            end
        end
    end

    # Should not reach here
    return EventRecord(deposits, :escaped, 0.0, 0.0, n_int)
end

# ─── Parse CLI ────────────────────────────────────────────────────────

function parse_args()
    isotope = "Bi214"
    n_events = 500000
    outdir = "output/histdiag"
    for i in 1:length(ARGS)
        if ARGS[i] == "--isotope" && i < length(ARGS)
            isotope = ARGS[i+1]
        elseif ARGS[i] == "--n-events" && i < length(ARGS)
            n_events = parse(Int, ARGS[i+1])
        elseif ARGS[i] == "--output" && i < length(ARGS)
            outdir = ARGS[i+1]
        end
    end
    return isotope, n_events, outdir
end

# ─── Main ─────────────────────────────────────────────────────────────

function main()
    isotope, n_events, outdir = parse_args()
    mkpath(outdir)

    if isotope == "Bi214"
        source = SourceConfig(
            label="histdiag_Bi214", E_MeV=2.448,
            entry=:barrel, angular=:flat,
            R_entry=72.8, H_entry=145.6, z_min_entry=0.0,
            gammas_per_yr=1.0e5,
        )
    else
        source = SourceConfig(
            label="histdiag_Tl208", E_MeV=2.615,
            entry=:barrel, angular=:flat,
            R_entry=72.8, H_entry=145.6, z_min_entry=0.0,
            gammas_per_yr=1.0e6, is_Tl208=true,
        )
    end

    params = Params(cut_E_visible_threshold_keV=10.0)
    geom = build_geometry(params)
    xcom = load_xcom(params.phys_xcom_data_path)
    rng = MersenneTwister(42)

    @printf("\n=== HISTOGRAM DIAGNOSTIC: %s, %d events ===\n\n", isotope, n_events)

    # Accumulators
    n_ss = 0; n_ms = 0; n_escaped = 0; n_ss_roi = 0

    ss_n_int = Int[]                      # n_interactions for SS events
    ss_dep_E_multi = Float64[]            # deposit energies for SS n>1 (inclusive)
    ss_dep_dz_multi = Float64[]           # |Δz| from cluster center for SS n>1
    ms_dep_E = Float64[]                  # deposit energies for MS (inclusive)
    ms_dep_dz = Float64[]                 # |Δz| from first deposit for MS
    ss_E_true = Float64[]                 # true cluster energy (SS)
    ss_E_smeared = Float64[]             # smeared cluster energy (SS)

    for i in 1:n_events
        evt = track_with_deposits(rng, source, geom, xcom, params)

        if evt.outcome == :escaped
            n_escaped += 1
            continue
        end

        if evt.outcome == :MS_rejected
            n_ms += 1
            # MS: inclusive deposit energies and distances from first deposit
            if length(evt.deposits) >= 1
                z0 = evt.deposits[1].z
                for d in evt.deposits
                    if d.E_keV >= params.cut_E_visible_threshold_keV
                        push!(ms_dep_E, d.E_keV)
                        push!(ms_dep_dz, abs(d.z - z0) * 10.0)  # mm
                    end
                end
            end
        else
            # SS event (in_ROI or outside_ROI)
            n_ss += 1
            push!(ss_n_int, evt.n_interactions)
            push!(ss_E_true, evt.E_cluster_true)
            push!(ss_E_smeared, evt.E_cluster_smeared)

            if evt.outcome == :SS_in_ROI
                n_ss_roi += 1
            end

            # SS n>1: inclusive deposit energies and distances
            if evt.n_interactions > 1
                z_cluster = mean(d.z for d in evt.deposits
                                 if d.E_keV >= params.cut_E_visible_threshold_keV)
                for d in evt.deposits
                    if d.E_keV >= params.cut_E_visible_threshold_keV
                        push!(ss_dep_E_multi, d.E_keV)
                        push!(ss_dep_dz_multi, abs(d.z - z_cluster) * 10.0)  # mm
                    end
                end
            end
        end

        if i % (n_events ÷ 10) == 0
            @printf("  Progress: %.0f%%\r", 100.0 * i / n_events)
        end
    end

    n_total_noescape = n_ss + n_ms
    Q_bb = params.phys_Q_betabeta_keV
    ROI_hw = params.cut_ROI_halfwidth_keV

    # ─── Print summary ────────────────────────────────────────────────
    @printf("\n\n=== SUMMARY (%s, %d events) ===\n", isotope, n_events)
    @printf("  Escaped:        %8d  (%.2f%%)\n", n_escaped, 100.0 * n_escaped / n_events)
    @printf("  SS total:       %8d  (%.2f%% of interacting)\n", n_ss,
            n_total_noescape > 0 ? 100.0 * n_ss / n_total_noescape : 0)
    @printf("  MS total:       %8d  (%.2f%% of interacting)\n", n_ms,
            n_total_noescape > 0 ? 100.0 * n_ms / n_total_noescape : 0)
    @printf("  SS/MS ratio:    %.3f\n", n_ms > 0 ? n_ss / n_ms : Inf)
    @printf("  SS in ROI:      %8d  (%.4f%% of total)\n", n_ss_roi,
            100.0 * n_ss_roi / n_events)
    @printf("  SS/total:       %.4f%%\n", 100.0 * n_ss / n_events)
    @printf("  SS_ROI/total:   %.6f%%\n", 100.0 * n_ss_roi / n_events)

    # Cross section comparison
    xcom_data = load_xcom(params.phys_xcom_data_path)
    E_gamma = source.E_MeV
    pe = σ_photo(xcom_data, E_gamma)
    co = σ_Compton(xcom_data, E_gamma)
    pp = σ_pair(xcom_data, E_gamma)
    tot = pe + co + pp
    @printf("\n  Cross sections at %.3f MeV:\n", E_gamma)
    @printf("    Photo:   %.1f%%\n", 100 * pe / tot)
    @printf("    Compton: %.1f%%\n", 100 * co / tot)
    @printf("    Pair:    %.1f%%\n", 100 * pp / tot)
    @printf("    → Expected SS if only photo: %.1f%%\n", 100 * pe / tot)
    @printf("    → Observed SS:               %.1f%%\n",
            n_total_noescape > 0 ? 100.0 * n_ss / n_total_noescape : 0)
    @printf("    → Excess SS from Compton chains within 3mm: %.1f%%\n",
            n_total_noescape > 0 ? 100.0 * n_ss / n_total_noescape - 100 * pe / tot : 0)

    # ─── Plot 1: SS vs MS bar chart ──────────────────────────────────
    p1 = bar(["SS", "MS"], [n_ss, n_ms];
             ylabel="counts", title="$isotope: SS vs MS",
             color=[:steelblue, :coral], legend=false)

    # ─── Plot 2: SS n_interactions ───────────────────────────────────
    if !isempty(ss_n_int)
        max_n = min(maximum(ss_n_int), 15)
        h_nint = [count(==(n), ss_n_int) for n in 1:max_n]
        p2 = bar(1:max_n, h_nint;
                 xlabel="n_interactions", ylabel="counts",
                 title="$isotope: SS interactions", color=:steelblue, legend=false)
    else
        p2 = plot(title="$isotope: SS interactions (no data)")
    end

    # ─── Plot 3: SS (n>1) deposit energies and |Δz| ─────────────────
    if !isempty(ss_dep_E_multi)
        p3a = histogram(ss_dep_E_multi; bins=100, xlabel="E_dep (keV)",
                        ylabel="counts", title="$isotope: SS(n>1) deposit E",
                        color=:steelblue, legend=false)
        p3b = histogram(ss_dep_dz_multi; bins=0:0.1:5.0, xlabel="|Δz| (mm)",
                        ylabel="counts", title="$isotope: SS(n>1) |Δz|",
                        color=:steelblue, legend=false)
        vline!(p3b, [params.cut_Δz_threshold_mm]; color=:red, linewidth=2,
               label="3mm cut", linestyle=:dash)
    else
        p3a = plot(title="SS(n>1) E (no data)")
        p3b = plot(title="SS(n>1) Δz (no data)")
    end

    # ─── Plot 4: MS deposit energies and |Δz| ────────────────────────
    if !isempty(ms_dep_E)
        p4a = histogram(ms_dep_E; bins=100, xlabel="E_dep (keV)",
                        ylabel="counts", title="$isotope: MS deposit E",
                        color=:coral, legend=false)
        p4b = histogram(ms_dep_dz; bins=0:0.5:50.0, xlabel="|Δz| (mm)",
                        ylabel="counts", title="$isotope: MS |Δz|",
                        color=:coral, legend=false)
        vline!(p4b, [params.cut_Δz_threshold_mm]; color=:red, linewidth=2,
               label="3mm cut", linestyle=:dash)
    else
        p4a = plot(title="MS E (no data)")
        p4b = plot(title="MS Δz (no data)")
    end

    # ─── Plot 5: SS energy spectrum (true + smeared) ─────────────────
    if !isempty(ss_E_true)
        roi_lo = Q_bb - ROI_hw
        roi_hi = Q_bb + ROI_hw
        E_lo = max(0, Q_bb - 500)
        E_hi = Q_bb + 500

        p5 = histogram(ss_E_true; bins=range(E_lo, E_hi, length=200),
                        xlabel="E (keV)", ylabel="counts",
                        title="$isotope: SS energy (true vs smeared)",
                        label="true", color=:steelblue, alpha=0.6)
        histogram!(p5, ss_E_smeared; bins=range(E_lo, E_hi, length=200),
                   label="smeared", color=:orange, alpha=0.5)
        vspan!(p5, [roi_lo, roi_hi]; fillalpha=0.15, fillcolor=:green, label="ROI")
        vline!(p5, [Q_bb]; color=:black, linewidth=1, linestyle=:dash, label="Q_ββ")
    else
        p5 = plot(title="SS energy (no data)")
    end

    # ─── Compose and save ────────────────────────────────────────────
    layout = @layout [
        a b
        c d
        e f
        g{0.5w}
    ]
    p_all = plot(p1, p2, p3a, p3b, p4a, p4b, p5;
                 layout=(4, 2), size=(1200, 1600),
                 plot_title="$isotope Diagnostic Histograms")

    outpath = joinpath(outdir, "histdiag_$(isotope).png")
    savefig(p_all, outpath)
    @printf("\n  Saved: %s\n", outpath)

    # Also save individual panels for clarity
    savefig(p5, joinpath(outdir, "energy_$(isotope).png"))
    @printf("  Saved: %s\n", joinpath(outdir, "energy_$(isotope).png"))
end

main()
