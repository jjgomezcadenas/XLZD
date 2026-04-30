# output.jl — HDF5, JSON, plots, summary.txt

"""
    integrate_in_FV_box(H, bin_edges_z, bin_edges_r2, params::Params) -> Int

Sum H[i,j] over bins whose centers fall inside the FV box.
"""
function integrate_in_FV_box(H::Matrix{Float64},
                             bin_edges_z::Vector{Float64},
                             bin_edges_r2::Vector{Float64},
                             params::Params)::Int
    total = 0.0
    for i in 1:size(H, 1)
        z_center = 0.5 * (bin_edges_z[i] + bin_edges_z[i+1])
        (z_center < params.fv_z_min_cm || z_center > params.fv_z_max_cm) && continue
        for j in 1:size(H, 2)
            r2_center = 0.5 * (bin_edges_r2[j] + bin_edges_r2[j+1])
            if r2_center ≤ params.fv_r2_max_cm2
                total += H[i, j]
            end
        end
    end
    round(Int, total)
end

# ── Derived quantities ──

function _compute_derived(result::Result)
    p = result.params
    s = result.source
    n_emitted = p.mc_N_samples

    n_ss_roi = result.counts[:SS_in_ROI]
    f_ss_roi = n_ss_roi / n_emitted

    # Background rate: f_SS_ROI × gammas_per_yr from the source CSV
    events_per_yr = f_ss_roi * s.gammas_per_yr

    return (;
        n_ss_roi, f_ss_roi, events_per_yr,
        n_emitted
    )
end

# ── HDF5 ──

"""
    write_h5(result::Result, path::String)

Write the HDF5 output file with geometry, source, MC params, physics,
counts, histograms, and derived calibration numbers.
"""
function write_h5(result::Result, path::String)
    p = result.params
    g = result.geometry
    s = result.source
    d = _compute_derived(result)

    h5open(path, "w") do f
        # /geometry/
        grp = create_group(f, "geometry")
        grp["R_lxe"] = g.R_lxe
        grp["L_lxe"] = g.L_lxe
        grp["rho_LXe"] = g.ρ_LXe
        grp["LXe_mass_kg"] = g.LXe_mass_kg
        grp["FV_LXe_mass_kg"] = g.FV_LXe_mass_kg

        # /source/
        grp = create_group(f, "source")
        grp["label"] = s.label
        grp["E_MeV"] = s.E_MeV
        grp["entry"] = string(s.entry)
        grp["angular"] = string(s.angular)
        grp["R_entry"] = s.R_entry
        grp["H_entry"] = s.H_entry
        grp["gammas_per_yr"] = s.gammas_per_yr
        if s.angular == :shaped
            grp["a"] = s.a
            grp["b"] = s.b
            grp["u_min"] = s.u_min
            grp["norm"] = s.norm
        end

        # /mc_params/
        grp = create_group(f, "mc_params")
        grp["n_emitted"] = p.mc_N_samples
        grp["n_threads"] = Threads.nthreads()
        grp["seed"] = p.mc_seed
        grp["runtime_seconds"] = result.runtime_seconds

        # /physics/
        grp = create_group(f, "physics")
        grp["E_gamma_initial_MeV"] = s.E_MeV
        grp["Q_bb_keV"] = p.phys_Q_betabeta_keV
        grp["sigma_E_over_E"] = p.phys_σ_E_over_E
        grp["ROI_halfwidth_keV"] = p.cut_ROI_halfwidth_keV
        grp["dz_threshold_mm"] = p.cut_Δz_threshold_mm
        grp["E_visible_threshold_keV"] = p.cut_E_visible_threshold_keV
        grp["E_tracking_cutoff_keV"] = p.cut_E_tracking_cutoff_keV

        # /fv/
        grp = create_group(f, "fv")
        grp["z_min_cm"] = p.fv_z_min_cm
        grp["z_max_cm"] = p.fv_z_max_cm
        grp["r2_max_cm2"] = p.fv_r2_max_cm2
        grp["bfv_z_min_cm"] = g.bfv.z_min
        grp["bfv_z_max_cm"] = g.bfv.z_max
        grp["bfv_r2_max_cm2"] = g.bfv.r2_max

        # /counts/
        grp = create_group(f, "counts")
        for (k, v) in result.counts
            grp["n_$(k)"] = v
        end

        # /histograms_2D/
        grp = create_group(f, "histograms_2D")
        grp["bin_edges_z"] = result.bin_edges_z
        grp["bin_edges_r2"] = result.bin_edges_r2
        grp["H_first_interaction"] = result.H_first_interaction
        grp["H_cluster_position"] = result.H_cluster_position
        grp["H_signal"] = result.H_signal
        grp["H_total_energy"] = result.H_total_energy

        # /histograms_1D/
        grp = create_group(f, "histograms_1D")
        grp["bin_edges_E_keV"] = result.bin_edges_E_keV
        grp["E_total_cluster_all_SS"] = result.E_total_cluster_all_SS
        grp["n_interactions_per_photon"] = result.n_interactions_per_photon
        grp["path_length_in_LXe_cm"] = result.path_length_in_LXe

        # /derived/
        grp = create_group(f, "derived")
        grp["f_SS_in_ROI"] = d.f_ss_roi
        grp["events_per_yr"] = d.events_per_yr
    end
end

# ── Summary text ──

"""
    write_summary(result::Result, path::String)

Write human-readable summary.txt.
"""
function write_summary(result::Result, path::String)
    p = result.params
    g = result.geometry
    s = result.source
    d = _compute_derived(result)
    n = p.mc_N_samples

    open(path, "w") do io
        @printf(io, "═══ XLZD Background MC — Summary ═══\n\n")

        @printf(io, "── Source ──\n")
        @printf(io, "  Label         = %s\n", s.label)
        @printf(io, "  E_γ           = %.3f MeV\n", s.E_MeV)
        @printf(io, "  Entry         = %s (R=%.1f cm)\n", s.entry, s.R_entry)
        @printf(io, "  Angular       = %s\n", s.angular)
        @printf(io, "  γ/yr (input)  = %.3e\n", s.gammas_per_yr)
        @printf(io, "\n")

        @printf(io, "── Geometry ──\n")
        @printf(io, "  R_lxe         = %.1f cm\n", g.R_lxe)
        @printf(io, "  L_lxe         = %.1f cm\n", g.L_lxe)
        @printf(io, "  LXe mass      = %.1f t\n", g.LXe_mass_kg / 1000)
        @printf(io, "  FV:  z ∈ [%.1f, %.1f] cm, r ≤ %.1f cm\n",
                p.fv_z_min_cm, p.fv_z_max_cm, sqrt(p.fv_r2_max_cm2))
        @printf(io, "  FV LXe mass   = %.1f t\n", g.FV_LXe_mass_kg / 1000)
        @printf(io, "  BFV: z ∈ [%.1f, %.1f] cm, r ≤ %.1f cm\n",
                g.bfv.z_min, g.bfv.z_max, sqrt(g.bfv.r2_max))
        @printf(io, "\n")

        @printf(io, "── Cuts ──\n")
        @printf(io, "  Q_ββ          = %.1f keV\n", p.phys_Q_betabeta_keV)
        @printf(io, "  ROI           = Q_ββ ± %.1f keV\n", p.cut_ROI_halfwidth_keV)
        @printf(io, "  Δz threshold  = %.1f mm\n", p.cut_Δz_threshold_mm)
        @printf(io, "  E_visible     = %.1f keV\n", p.cut_E_visible_threshold_keV)
        @printf(io, "  E_cutoff      = %.1f keV\n", p.cut_E_tracking_cutoff_keV)
        @printf(io, "  σ_E/E         = %.4f\n", p.phys_σ_E_over_E)
        @printf(io, "\n")

        @printf(io, "── MC run ──\n")
        @printf(io, "  N_emitted     = %d\n", n)
        @printf(io, "  N_threads     = %d\n", Threads.nthreads())
        @printf(io, "  Seed          = %d\n", p.mc_seed)
        @printf(io, "  Runtime       = %.2f s\n", result.runtime_seconds)
        @printf(io, "\n")

        @printf(io, "── Outcome counts ──\n")
        for k in [:outside_bfv, :escaped, :MS_rejected, :SS_outside_FV,
                  :SS_outside_ROI, :companion_vetoed, :skin_vetoed, :SS_in_ROI]
            v = get(result.counts, k, 0)
            @printf(io, "  %-20s %12d  (%7.4f%%)\n", k, v, 100.0 * v / n)
        end
        total = sum(values(result.counts))
        @printf(io, "  %-20s %12d\n", "TOTAL", total)
        @printf(io, "\n")

        @printf(io, "── Result ──\n")
        @printf(io, "  f(SS-in-ROI)       = %.6e  (per emitted gamma)\n", d.f_ss_roi)
        @printf(io, "  γ/yr from source   = %.3e\n", s.gammas_per_yr)
        @printf(io, "  Background rate    = %.4e events/yr\n", d.events_per_yr)
        @printf(io, "═══════════════════════════════════════════════════════\n")
    end
end

# ── Trajectories JSON ──

function write_trajectories_json(result::Result, path::String)
    json_str = to_json(result.trajectories, result.trajectories_fv,
                       result.geometry, result.params)
    open(path, "w") do io
        write(io, json_str)
    end
end

# ── Viewer copy ──

function copy_viewer(template_path::String, output_path::String)
    cp(template_path, output_path; force=true)
end

# ── Plots ──

function _heatmap_with_fv(H, bin_edges_z, bin_edges_r2, params, title_str)
    z_centers  = 0.5 .* (bin_edges_z[1:end-1] .+ bin_edges_z[2:end])
    r2_centers = 0.5 .* (bin_edges_r2[1:end-1] .+ bin_edges_r2[2:end])

    gr()
    p = heatmap(z_centers, r2_centers, H';
                xlabel="z (cm)", ylabel="r² (cm²)",
                title=title_str, color=:viridis,
                colorbar_title="counts")

    fv_z = [params.fv_z_min_cm, params.fv_z_max_cm, params.fv_z_max_cm,
            params.fv_z_min_cm, params.fv_z_min_cm]
    fv_r2 = [0.0, 0.0, params.fv_r2_max_cm2, params.fv_r2_max_cm2, 0.0]
    plot!(p, fv_z, fv_r2; linecolor=:green, linewidth=2, label="FV", linestyle=:dash)

    return p
end

function plot_heatmap_signal(result::Result, path::String)
    d = _compute_derived(result)
    title = @sprintf("SS-in-ROI (N=%d) — %s", d.n_ss_roi, result.source.label)
    p = _heatmap_with_fv(result.H_signal, result.bin_edges_z,
                         result.bin_edges_r2, result.params, title)
    savefig(p, path)
end

function plot_heatmap_first_interaction(result::Result, path::String)
    title = @sprintf("First interaction — %s", result.source.label)
    p = _heatmap_with_fv(result.H_first_interaction, result.bin_edges_z,
                         result.bin_edges_r2, result.params, title)
    savefig(p, path)
end

function plot_spectrum_SS(result::Result, path::String)
    pr = result.params
    edges = result.bin_edges_E_keV
    centers = 0.5 .* (edges[1:end-1] .+ edges[2:end])

    gr()
    p = plot(centers, result.E_total_cluster_all_SS;
             xlabel="E (keV)", ylabel="counts",
             title="SS energy spectrum — $(result.source.label)",
             label="all SS", linecolor=:blue, legend=:topleft)

    roi_lo = pr.phys_Q_betabeta_keV - pr.cut_ROI_halfwidth_keV
    roi_hi = pr.phys_Q_betabeta_keV + pr.cut_ROI_halfwidth_keV
    vspan!(p, [roi_lo, roi_hi]; fillalpha=0.2, fillcolor=:green, label="ROI")

    savefig(p, path)
end

function plot_diagnostic(result::Result, path::String)
    pr = result.params
    gr()

    labels = ["out_bfv", "escaped", "MS_rej", "SS_noFV", "SS_noROI", "comp_veto", "skin_veto", "SS_ROI"]
    syms   = [:outside_bfv, :escaped, :MS_rejected, :SS_outside_FV,
              :SS_outside_ROI, :companion_vetoed, :skin_vetoed, :SS_in_ROI]
    vals   = [get(result.counts, s, 0) for s in syms]
    p1 = bar(labels, vals; ylabel="counts", title="Outcomes",
             legend=false, color=:steelblue, xrotation=30)

    n_bins_int = length(result.n_interactions_per_photon)
    p2 = bar(1:n_bins_int, result.n_interactions_per_photon;
             xlabel="n_interactions", ylabel="counts",
             title="Interactions per photon", legend=false, color=:coral)

    n_bins_pl = length(result.path_length_in_LXe)
    pl_edges = range(0.0, 500.0, length=n_bins_pl + 1)
    pl_centers = 0.5 .* (pl_edges[1:end-1] .+ pl_edges[2:end])
    p3 = bar(pl_centers, result.path_length_in_LXe;
             xlabel="path length (cm)", ylabel="counts",
             title="Path length in LXe", legend=false, color=:mediumpurple)

    z_c  = 0.5 .* (result.bin_edges_z[1:end-1] .+ result.bin_edges_z[2:end])
    r2_c = 0.5 .* (result.bin_edges_r2[1:end-1] .+ result.bin_edges_r2[2:end])
    p4 = heatmap(z_c, r2_c, result.H_total_energy';
                 xlabel="z (cm)", ylabel="r² (cm²)",
                 title="Energy deposits", color=:inferno)

    p_all = plot(p1, p2, p3, p4; layout=(2, 2), size=(1000, 800))
    savefig(p_all, path)
end

# ── Top-level driver ──

"""
    write_outputs(result::Result)

Top-level output driver: creates output dir, then calls all writers and plotters.
"""
function write_outputs(result::Result)
    dir = result.params.out_dir
    mkpath(dir)

    println("Writing HDF5...")
    write_h5(result, joinpath(dir, "results.h5"))

    println("Writing summary...")
    write_summary(result, joinpath(dir, "summary.txt"))

    println("Writing trajectories JSON...")
    write_trajectories_json(result, joinpath(dir, "trajectories.json"))

    println("Copying viewer...")
    try
        copy_viewer(result.params.out_viewer_template, joinpath(dir, "viewer.html"))
    catch e
        @warn "Could not copy viewer: $e"
    end

    println("Plotting...")
    plot_heatmap_signal(result, joinpath(dir, "heatmap_signal.png"))
    plot_heatmap_first_interaction(result, joinpath(dir, "heatmap_first_interaction.png"))
    plot_spectrum_SS(result, joinpath(dir, "spectrum_SS.png"))
    plot_diagnostic(result, joinpath(dir, "diagnostic.png"))

    # Print final result
    d = _compute_derived(result)
    @printf("\n═══ RESULT: %s ═══\n", result.source.label)
    @printf("  f(SS-in-ROI) = %.6e\n", d.f_ss_roi)
    @printf("  γ/yr         = %.3e\n", result.source.gammas_per_yr)
    @printf("  Background   = %.4e events/yr\n", d.events_per_yr)
    @printf("═══════════════════════════════════════\n")

    println("All outputs written to $(dir)")
end
