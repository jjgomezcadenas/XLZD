"""
Diagnostic script: trace individual gamma events step by step.

Runs a small number of events with NO BFV cut (all first interactions
accepted) and prints detailed tracking information for every SS event.
Helps understand:
  - Why Bi-214 SS_in_ROI rate seems high
  - Why Tl-208 produces zero SS_in_ROI events

Usage:
  julia --project=. scripts/diagnose.jl --isotope Bi214 --n-events 10000
  julia --project=. scripts/diagnose.jl --isotope Tl208 --n-events 100000
"""

using XLZD, Random, Printf

# ─── Configuration ────────────────────────────────────────────────────

function parse_args()
    isotope = "Bi214"
    n_events = 10000
    for i in 1:length(ARGS)
        if ARGS[i] == "--isotope" && i < length(ARGS)
            isotope = ARGS[i+1]
        elseif ARGS[i] == "--n-events" && i < length(ARGS)
            n_events = parse(Int, ARGS[i+1])
        end
    end
    return isotope, n_events
end

# ─── Verbose single-photon tracker ────────────────────────────────────

"""
Track one gamma with full diagnostic output. No BFV cut.
Returns (outcome, details_string).
"""
function track_verbose(rng, source, geom, xcom, params)
    E_visible = params.cut_E_visible_threshold_keV / 1000.0
    E_cutoff  = params.cut_E_tracking_cutoff_keV / 1000.0
    Δz_thresh = params.cut_Δz_threshold_mm / 10.0
    Q_bb      = params.phys_Q_betabeta_keV
    ROI_hw    = params.cut_ROI_halfwidth_keV

    # Sample entry
    x, y, z = sample_entry_point(rng, source, geom.L_lxe)
    dx, dy, dz = sample_entry_direction(rng, source, x, y, z, geom.L_lxe)
    E = source.E_MeV

    lines = String[]
    push!(lines, @sprintf("  Entry: (%.1f, %.1f, %.1f) cm, r=%.1f cm", x, y, z, sqrt(x^2+y^2)))
    push!(lines, @sprintf("  Dir:   (%.3f, %.3f, %.3f)", dx, dy, dz))
    push!(lines, @sprintf("  E_initial = %.4f MeV", E))

    E_cluster = 0.0
    z_cluster = NaN; x_cluster = NaN; y_cluster = NaN
    cluster_started = false
    n_int = 0
    outcome = :in_progress

    while outcome == :in_progress
        μ_t = μ_total_lin(xcom, E, geom.ρ_LXe)
        d_int = -log(rand(rng)) / μ_t
        d_boundary = path_to_cylinder_exit(x, y, z, dx, dy, dz,
                                           geom.R_lxe, geom.L_lxe)

        if d_int >= d_boundary
            x += d_boundary * dx; y += d_boundary * dy; z += d_boundary * dz
            push!(lines, @sprintf("  → ESCAPED at (%.1f, %.1f, %.1f), r=%.1f, after %.1f cm",
                                  x, y, z, sqrt(x^2+y^2), d_boundary))
            if cluster_started
                E_smeared_keV = (E_cluster + params.phys_σ_E_over_E * E_cluster * randn(rng)) * 1000.0
                push!(lines, @sprintf("  Cluster E = %.1f keV (smeared %.1f), ROI [%.1f, %.1f]",
                                      E_cluster*1000, E_smeared_keV, Q_bb - ROI_hw, Q_bb + ROI_hw))
                if abs(E_smeared_keV - Q_bb) ≤ ROI_hw
                    outcome = :SS_in_ROI
                else
                    outcome = :SS_outside_ROI
                end
            else
                outcome = :escaped
            end
            break
        end

        x += d_int * dx; y += d_int * dy; z += d_int * dz
        n_int += 1
        r = sqrt(x^2 + y^2)

        # Interaction type
        σ_pe = σ_photo(xcom, E)
        σ_co = σ_Compton(xcom, E)
        σ_pp = σ_pair(xcom, E)
        σ_tot = σ_pe + σ_co + σ_pp
        r_branch = rand(rng)

        if r_branch < σ_pp / σ_tot
            push!(lines, @sprintf("  [%d] PAIR at (%.1f,%.1f,%.1f) r=%.1f, E=%.4f MeV",
                                  n_int, x, y, z, r, E))
            outcome = :MS_rejected
            break

        elseif r_branch < (σ_pp + σ_pe) / σ_tot
            E_dep = E
            push!(lines, @sprintf("  [%d] PHOTO at (%.1f,%.1f,%.1f) r=%.1f, E_dep=%.4f MeV",
                                  n_int, x, y, z, r, E_dep))

            if E_dep >= E_visible
                if !cluster_started
                    z_cluster = z; x_cluster = x; y_cluster = y
                    E_cluster = E_dep; cluster_started = true
                    push!(lines, @sprintf("      → cluster started, E=%.1f keV", E_cluster*1000))
                else
                    dz_sep = abs(z - z_cluster) * 10  # mm
                    if dz_sep < Δz_thresh * 10
                        E_cluster += E_dep
                        push!(lines, @sprintf("      → added to cluster (Δz=%.1f mm), E=%.1f keV",
                                              dz_sep, E_cluster*1000))
                    else
                        push!(lines, @sprintf("      → MS rejected (Δz=%.1f mm)", dz_sep))
                        outcome = :MS_rejected; break
                    end
                end
            end

            if outcome == :in_progress
                if cluster_started
                    E_smeared_keV = (E_cluster + params.phys_σ_E_over_E * E_cluster * randn(rng)) * 1000.0
                    push!(lines, @sprintf("  Cluster E = %.1f keV (smeared %.1f), ROI [%.1f, %.1f]",
                                          E_cluster*1000, E_smeared_keV, Q_bb - ROI_hw, Q_bb + ROI_hw))
                    outcome = abs(E_smeared_keV - Q_bb) ≤ ROI_hw ? :SS_in_ROI : :SS_outside_ROI
                else
                    outcome = :escaped
                end
            end
            break

        else  # Compton
            E_scatt, cos_θ = sample_klein_nishina(rng, E)
            E_dep = E - E_scatt
            push!(lines, @sprintf("  [%d] COMPTON at (%.1f,%.1f,%.1f) r=%.1f, E_dep=%.1f keV, E_scatt=%.4f MeV, cosθ=%.3f",
                                  n_int, x, y, z, r, E_dep*1000, E_scatt, cos_θ))

            if E_dep >= E_visible
                if !cluster_started
                    z_cluster = z; x_cluster = x; y_cluster = y
                    E_cluster = E_dep; cluster_started = true
                    push!(lines, @sprintf("      → cluster started, E=%.1f keV", E_cluster*1000))
                else
                    dz_sep = abs(z - z_cluster) * 10
                    if dz_sep < Δz_thresh * 10
                        E_cluster += E_dep
                        push!(lines, @sprintf("      → added to cluster (Δz=%.1f mm), E=%.1f keV",
                                              dz_sep, E_cluster*1000))
                    else
                        push!(lines, @sprintf("      → MS rejected (Δz=%.1f mm)", dz_sep))
                        outcome = :MS_rejected; break
                    end
                end
            else
                push!(lines, @sprintf("      → invisible (%.1f keV < %.1f keV threshold)",
                                      E_dep*1000, E_visible*1000))
            end

            φ_az = 2π * rand(rng)
            dx, dy, dz = rotate_direction(dx, dy, dz, cos_θ, φ_az)
            E = E_scatt

            if E < E_cutoff
                push!(lines, @sprintf("  [%d] CUTOFF E=%.1f keV deposited locally", n_int, E*1000))
                if E >= E_visible
                    if !cluster_started
                        z_cluster = z; x_cluster = x; y_cluster = y
                        E_cluster = E; cluster_started = true
                    else
                        dz_sep = abs(z - z_cluster) * 10
                        if dz_sep < Δz_thresh * 10
                            E_cluster += E
                        else
                            outcome = :MS_rejected; break
                        end
                    end
                end
                if outcome == :in_progress
                    if cluster_started
                        E_smeared_keV = (E_cluster + params.phys_σ_E_over_E * E_cluster * randn(rng)) * 1000.0
                        push!(lines, @sprintf("  Cluster E = %.1f keV (smeared %.1f), ROI [%.1f, %.1f]",
                                              E_cluster*1000, E_smeared_keV, Q_bb - ROI_hw, Q_bb + ROI_hw))
                        outcome = abs(E_smeared_keV - Q_bb) ≤ ROI_hw ? :SS_in_ROI : :SS_outside_ROI
                    else
                        outcome = :escaped
                    end
                end
                break
            end
        end
    end

    push!(lines, @sprintf("  OUTCOME: %s (n_int=%d)", outcome, n_int))
    return outcome, join(lines, "\n")
end

# ─── Main ─────────────────────────────────────────────────────────────

function main()
    isotope, n_events = parse_args()

    if isotope == "Bi214"
        source = SourceConfig(
            label="diag_Bi214", E_MeV=2.448,
            entry=:barrel, angular=:flat,
            R_entry=72.8, H_entry=145.6, z_min_entry=0.0,
            gammas_per_yr=1.0e5,
        )
    else
        source = SourceConfig(
            label="diag_Tl208", E_MeV=2.615,
            entry=:barrel, angular=:flat,
            R_entry=72.8, H_entry=145.6, z_min_entry=0.0,
            gammas_per_yr=1.0e6,
            is_Tl208=true,
        )
    end

    params = Params(
        cut_E_visible_threshold_keV = 10.0,
    )
    geom = build_geometry(params)
    xcom = load_xcom(params.phys_xcom_data_path)
    rng = MersenneTwister(42)

    # Counters
    counts = Dict{Symbol,Int}()
    n_ss_total = 0
    n_ss_roi = 0
    n_ss_outside = 0

    # Energy distribution of SS events
    ss_energies = Float64[]

    @printf("\n=== DIAGNOSTIC: %s, %d events, E_visible=%.0f keV ===\n\n",
            isotope, n_events, params.cut_E_visible_threshold_keV)

    n_printed = 0
    max_print = 50  # print first N SS events in detail

    for i in 1:n_events
        outcome, details = track_verbose(rng, source, geom, xcom, params)
        counts[outcome] = get(counts, outcome, 0) + 1

        if outcome == :SS_in_ROI || outcome == :SS_outside_ROI
            n_ss_total += 1
            # Extract cluster energy from details
            m = match(r"Cluster E = ([\d.]+) keV", details)
            if m !== nothing
                E_keV = parse(Float64, m.captures[1])
                push!(ss_energies, E_keV)
            end

            if outcome == :SS_in_ROI
                n_ss_roi += 1
                if n_printed < max_print
                    @printf("── SS_in_ROI event #%d (event %d) ──\n", n_ss_roi, i)
                    println(details)
                    println()
                    n_printed += 1
                end
            else
                n_ss_outside += 1
            end
        end
    end

    # Summary
    @printf("\n=== SUMMARY (%d events) ===\n", n_events)
    for (k, v) in sort(collect(counts))
        @printf("  %-20s %8d  (%.4f%%)\n", k, v, 100.0 * v / n_events)
    end

    @printf("\n  SS total:       %d\n", n_ss_total)
    @printf("  SS in ROI:      %d\n", n_ss_roi)
    @printf("  SS outside ROI: %d\n", n_ss_outside)

    if !isempty(ss_energies)
        @printf("\n  SS energy distribution:\n")
        @printf("    min:    %.1f keV\n", minimum(ss_energies))
        @printf("    max:    %.1f keV\n", maximum(ss_energies))
        @printf("    mean:   %.1f keV\n", sum(ss_energies)/length(ss_energies))

        # Histogram in 50 keV bins near Q_bb
        Q = params.phys_Q_betabeta_keV
        @printf("\n  Energy bins near Q_ββ (%.0f keV):\n", Q)
        for lo in (Q-200):50:(Q+200)
            hi = lo + 50
            n = count(e -> lo ≤ e < hi, ss_energies)
            bar = repeat("█", min(n, 60))
            @printf("    [%7.0f, %7.0f) keV: %5d %s\n", lo, hi, n, bar)
        end

        # Zoom into ROI
        roi_lo = Q - params.cut_ROI_halfwidth_keV
        roi_hi = Q + params.cut_ROI_halfwidth_keV
        n_in_roi_true = count(e -> roi_lo ≤ e ≤ roi_hi, ss_energies)
        @printf("\n  True energy in ROI [%.1f, %.1f]: %d (before smearing)\n",
                roi_lo, roi_hi, n_in_roi_true)
        @printf("  ROI half-width: %.1f keV (±1σ at σ/E=%.1f%%)\n",
                params.cut_ROI_halfwidth_keV, params.phys_σ_E_over_E*100)
    end
end

main()
