# mc.jl — Main MC loop, threading, single-photon tracker

"""
    ThreadResult

Per-thread accumulator: outcome counts, 2D/1D histograms, trajectory buffer.
"""
mutable struct ThreadResult
    counts::Dict{Symbol, Int}
    H_first_interaction::Matrix{Float64}
    H_cluster_position::Matrix{Float64}
    H_signal::Matrix{Float64}
    H_total_energy::Matrix{Float64}
    E_total_cluster_all_SS::Vector{Float64}
    n_interactions_per_photon::Vector{Float64}
    path_length_in_LXe::Vector{Float64}
    traj_buffer::TrajectoryBuffer
    traj_buffer_fv::TrajectoryBuffer
end

"""
    ThreadResult(params::Params) -> ThreadResult

Create a zeroed-out per-thread accumulator with histograms sized from `params`.
"""
function ThreadResult(params::Params)::ThreadResult
    nz  = params.hist_n_z_bins
    nr2 = params.hist_n_r2_bins
    nE  = params.hist_n_E_bins

    counts = Dict{Symbol, Int}(
        :outside_bfv => 0, :escaped => 0, :MS_rejected => 0,
        :SS_outside_ROI => 0, :SS_in_ROI => 0, :SS_outside_FV => 0,
        :companion_vetoed => 0, :skin_vetoed => 0
    )

    ThreadResult(
        counts,
        zeros(Float64, nz, nr2),   # H_first_interaction
        zeros(Float64, nz, nr2),   # H_cluster_position
        zeros(Float64, nz, nr2),   # H_signal
        zeros(Float64, nz, nr2),   # H_total_energy
        zeros(Float64, nE),        # E_total_cluster_all_SS
        zeros(Float64, 50),        # n_interactions_per_photon
        zeros(Float64, 100),       # path_length_in_LXe
        TrajectoryBuffer(params.mc_n_traj_per_outcome),
        TrajectoryBuffer(params.mc_n_traj_per_outcome)
    )
end

# ── Histogram helpers ──

"""Compute bin index (1-based) for value in [lo, hi) with n bins. Returns 0 if out of range."""
@inline function _bin_index(val::Float64, lo::Float64, hi::Float64, n::Int)::Int
    val < lo && return 0
    val >= hi && return 0
    return clamp(floor(Int, (val - lo) / (hi - lo) * n) + 1, 1, n)
end

@inline function _fill_2d!(H::Matrix{Float64}, z::Float64, r2::Float64, params::Params)
    iz  = _bin_index(z,  params.hist_z_min_cm,   params.hist_z_max_cm,   params.hist_n_z_bins)
    ir2 = _bin_index(r2, params.hist_r2_min_cm2,  params.hist_r2_max_cm2, params.hist_n_r2_bins)
    if iz > 0 && ir2 > 0
        H[iz, ir2] += 1.0
    end
end

@inline function _fill_2d_weighted!(H::Matrix{Float64}, z::Float64, r2::Float64,
                                     w::Float64, params::Params)
    iz  = _bin_index(z,  params.hist_z_min_cm,   params.hist_z_max_cm,   params.hist_n_z_bins)
    ir2 = _bin_index(r2, params.hist_r2_min_cm2,  params.hist_r2_max_cm2, params.hist_n_r2_bins)
    if iz > 0 && ir2 > 0
        H[iz, ir2] += w
    end
end

# ── SS energy classification ──

"""
Return :SS_in_ROI or :SS_outside_ROI based on cluster energy vs ROI window.
The true deposited energy is smeared by the detector resolution (Gaussian,
σ_E/E = phys_σ_E_over_E) before the ROI check.
"""
@inline function classify_ss_energy(E_cluster::Float64, rng::AbstractRNG,
                                    params::Params)::Symbol
    σ_E = params.phys_σ_E_over_E * E_cluster
    E_smeared_keV = (E_cluster + σ_E * randn(rng)) * 1000.0
    if abs(E_smeared_keV - params.phys_Q_betabeta_keV) ≤ params.cut_ROI_halfwidth_keV
        return :SS_in_ROI
    else
        return :SS_outside_ROI
    end
end

"""
    track_companion_gamma(rng, x0, y0, z0, E_companion, geom, xcom, params) -> Bool

Track the ²⁰⁸Tl companion gamma (583/860/763 keV) emitted from the same
source point as the 2.6 MeV gamma, in a random direction (full 4π).
Returns true if the companion is VISIBLE (event should be vetoed).

The companion is visible if:
  - It deposits ≥ companion_veto_keV in the active LXe (inside FC), OR
  - It deposits ≥ skin_veto_keV in the LXe skin region

The companion is invisible (returns false) if:
  - It exits the LXe cylinder without depositing above either threshold
  - It is absorbed in the Ti walls (never enters LXe)
  - E_companion is 0.0 (no companion emitted, 1% of decays)
"""
function track_companion_gamma(rng::AbstractRNG,
                                x0::Float64, y0::Float64, z0::Float64,
                                E_companion::Float64,
                                geom::Geometry, xcom::XCOMTable,
                                params::Params)::Bool
    # No companion emitted
    E_companion <= 0.0 && return false

    E_comp_veto = params.cut_companion_veto_keV / 1000.0  # MeV
    E_skin_veto = params.cut_skin_veto_keV / 1000.0       # MeV
    E_cutoff    = params.cut_E_tracking_cutoff_keV / 1000.0

    # Emit in random direction (full 4π)
    u = 2.0 * rand(rng) - 1.0
    φ = 2π * rand(rng)
    s = sqrt(1.0 - u * u)
    dx = s * cos(φ);  dy = s * sin(φ);  dz = u

    x = x0; y = y0; z = z0
    E = E_companion

    # Accumulate deposits in active LXe and skin separately
    E_dep_active = 0.0
    E_dep_skin = 0.0

    while true
        μ_t = μ_total_lin(xcom, E, geom.ρ_LXe)
        d_int = -log(rand(rng)) / μ_t
        # Companion tracked in full LXe volume (TPC + skin)
        d_boundary = path_to_cylinder_exit(x, y, z, dx, dy, dz,
                                           geom.R_lxe_full, 0.0, geom.L_lxe_full)

        if d_int >= d_boundary
            # Companion escapes full LXe — check accumulated deposits
            break
        end

        # Move to interaction point
        x += d_int * dx;  y += d_int * dy;  z += d_int * dz

        # Determine interaction type
        σ_pe = σ_photo(xcom, E)
        σ_co = σ_Compton(xcom, E)
        σ_pp = σ_pair(xcom, E)
        σ_tot = σ_pe + σ_co + σ_pp
        r_branch = rand(rng)

        if r_branch < σ_pp / σ_tot
            # Pair production at companion energy is below threshold (< 1.022 MeV)
            # but handle gracefully: deposit full energy
            E_dep = E
            if in_skin(x, y, z, params, geom.L_lxe)
                E_dep_skin += E_dep
            elseif in_active_lxe(x, y, z, params, geom.L_lxe)
                E_dep_active += E_dep
            end
            break

        elseif r_branch < (σ_pp + σ_pe) / σ_tot
            # Photoelectric: full absorption
            E_dep = E
            if in_skin(x, y, z, params, geom.L_lxe)
                E_dep_skin += E_dep
            elseif in_active_lxe(x, y, z, params, geom.L_lxe)
                E_dep_active += E_dep
            end
            break

        else
            # Compton scatter
            E_scatt, cos_θ = sample_klein_nishina(rng, E)
            E_dep = E - E_scatt

            if in_skin(x, y, z, params, geom.L_lxe)
                E_dep_skin += E_dep
            elseif in_active_lxe(x, y, z, params, geom.L_lxe)
                E_dep_active += E_dep
            end

            # Early exit: if either threshold already exceeded, companion is visible
            if E_dep_active >= E_comp_veto || E_dep_skin >= E_skin_veto
                return true
            end

            # Update direction and energy
            φ_az = 2π * rand(rng)
            dx, dy, dz = rotate_direction(dx, dy, dz, cos_θ, φ_az)
            E = E_scatt

            # Energy cutoff: deposit remaining locally
            if E < E_cutoff
                if in_skin(x, y, z, params, geom.L_lxe)
                    E_dep_skin += E
                elseif in_active_lxe(x, y, z, params, geom.L_lxe)
                    E_dep_active += E
                end
                break
            end
        end
    end

    # Check accumulated deposits against thresholds
    return E_dep_active >= E_comp_veto || E_dep_skin >= E_skin_veto
end

"""
    Result

Merged MC result across all threads. Contains params, geometry, source config,
summed counts, summed histograms, merged trajectories, and total runtime.
"""
struct Result
    params::Params
    geometry::Geometry
    source::SourceConfig
    xcom::XCOMTable
    counts::Dict{Symbol, Int}
    H_first_interaction::Matrix{Float64}
    H_cluster_position::Matrix{Float64}
    H_signal::Matrix{Float64}
    H_total_energy::Matrix{Float64}
    E_total_cluster_all_SS::Vector{Float64}
    n_interactions_per_photon::Vector{Float64}
    path_length_in_LXe::Vector{Float64}
    bin_edges_z::Vector{Float64}
    bin_edges_r2::Vector{Float64}
    bin_edges_E_keV::Vector{Float64}
    trajectories::Vector{Trajectory}
    trajectories_fv::Vector{Trajectory}
    runtime_seconds::Float64
end

"""
    track_one_photon!(state::ThreadResult, rng, geom::Geometry,
                      source::SourceConfig, xcom::XCOMTable,
                      params::Params) -> Symbol

Simulate one gamma entering the LXe from a configured source and track it
through the liquid xenon until it escapes, is absorbed, or falls below
the energy cutoff. Returns a Symbol classifying the outcome.

## Algorithm

**1. Entry point.**
Sampled uniformly on the source's entry surface (barrel cylinder or
endcap disc). Unlike the old ring code, the gamma starts at the LXe
boundary, not inside a copper ring.

**2. Inward direction.**
u = cos θ is sampled from the source's angular distribution:
  - Flat: u ~ Uniform(0, 1) — for field cage and PMT sources.
  - Shaped: inverse CDF of linear fit — for cryostat sources that
    traversed Ti walls and LXe skin.
All gammas point inward by construction. No backward rejection needed.

**3. BFV pre-check on first interaction.**
After the first interaction, check if the point is inside the BFV
(FV + 1 cm buffer). If outside → :outside_bfv, stop immediately.
This avoids tracking gammas that interact far from the fiducial volume.

**4. Tracking loop.**
Identical to the original code: exponential step length, ray-cylinder
intersection, photoelectric / Compton / pair branching, Klein-Nishina
scattering, energy cutoff, SS/MS Δz clustering.

**5. FV check on final cluster.**
After tracking completes with a single-site event, check if the cluster
position is inside the FV (strict). If outside → :SS_outside_FV.
If inside → apply energy smearing and ROI check.

## Outcomes

| Symbol            | Meaning                                             |
|:------------------|:----------------------------------------------------|
| `:outside_bfv`    | First interaction outside Buffer Fiducial Volume     |
| `:escaped`        | Entered LXe but no visible energy deposited          |
| `:MS_rejected`    | Multiple deposits separated by |Δz| ≥ threshold     |
| `:SS_outside_FV`  | Single-site, but cluster outside Fiducial Volume     |
| `:SS_outside_ROI` | Single-site in FV, but energy outside ROI window     |
| `:SS_in_ROI`      | Single-site in FV, energy in ROI — the signal        |
"""
function track_one_photon!(state::ThreadResult, rng::AbstractRNG,
                           geom::Geometry, source::SourceConfig,
                           xcom::XCOMTable, params::Params)::Symbol
    # Convert cut parameters to MeV/cm for inner loop
    E_visible  = params.cut_E_visible_threshold_keV / 1000.0
    E_cutoff   = params.cut_E_tracking_cutoff_keV / 1000.0
    Δz_thresh  = params.cut_Δz_threshold_mm / 10.0

    # 1. Sample entry point on the source surface
    x, y, z = sample_entry_point(rng, source, geom.L_lxe)

    # 2. Sample inward direction from angular distribution
    dx, dy, dz = sample_entry_direction(rng, source, x, y, z, geom.L_lxe)

    # 3. Initialize photon state
    E = source.E_MeV
    E_cluster = 0.0
    z_cluster = NaN;  x_cluster = NaN;  y_cluster = NaN
    cluster_started = false
    E_skin_deposit = 0.0   # accumulated energy in skin region
    outcome = :in_progress
    n_int = 0
    total_path = 0.0
    first_int_x = NaN;  first_int_y = NaN;  first_int_z = NaN
    bfv_checked = false

    # Trajectory recording
    has_buffer_space = buffer_space_available(state.traj_buffer, :SS_in_ROI) ||
                buffer_space_available(state.traj_buffer, :SS_outside_ROI) ||
                buffer_space_available(state.traj_buffer, :MS_rejected) ||
                buffer_space_available(state.traj_buffer, :escaped) ||
                buffer_space_available(state.traj_buffer, :outside_bfv)
    traj_points = InteractionPoint[]
    if has_buffer_space
        push!(traj_points, InteractionPoint(x, y, z, E, :entry, 0.0))
    end

    # Tracking volume: for extended sources (PMTs, cryostat), expand
    # the cylinder to include the source position.
    #   - Radially: R_lxe_full (includes skin) for cryostat barrel sources
    #   - Axially: extend below cathode for bottom PMTs, above gate for top PMTs
    # Visibility rules:
    #   - Deposits in skin (r > R_skin_inner): apply skin_veto_threshold → :skin_vetoed
    #   - Deposits in TPC drift (r ≤ R_skin_inner, 0 ≤ z ≤ L_lxe): normal SS/MS logic
    #   - Deposits below cathode (z < 0) or above gate (z > L_lxe): invisible
    if source.use_extended_volume
        track_R = geom.R_lxe_full                       # include skin radially
        track_z_min = min(geom.z_bottom, source.z_entry) # bottom head or PMT
        track_z_max = max(geom.z_top, source.z_entry)    # top head or PMT or gate
    else
        track_R = geom.R_lxe
        track_z_min = 0.0
        track_z_max = geom.L_lxe
    end

    # Skin veto threshold in MeV
    E_skin_veto = params.cut_skin_veto_keV / 1000.0
    R_skin_inner2 = params.geom_R_skin_inner^2  # = 74.3² cm²

    # 4a. If gamma starts in gas (above liquid surface), propagate to
    # LXe surface in a straight line (no interactions in gas).
    if z > geom.L_lxe && dz < 0.0
        d_to_surface = (z - geom.L_lxe) / (-dz)
        x += d_to_surface * dx
        y += d_to_surface * dy
        z = geom.L_lxe  # exactly at liquid surface
        total_path += d_to_surface
        # Check if gamma went outside radially during gas propagation
        if x*x + y*y > track_R * track_R
            outcome = :escaped
        end
    elseif z > geom.L_lxe && dz >= 0.0
        # Gamma above liquid surface going up — escapes
        outcome = :escaped
    end

    # 4b. Tracking loop in LXe
    while outcome == :in_progress
        μ_t = μ_total_lin(xcom, E, geom.ρ_LXe)
        d_int = -log(rand(rng)) / μ_t
        d_boundary = path_to_cylinder_exit(x, y, z, dx, dy, dz,
                                           track_R, track_z_min, track_z_max)

        if d_int ≥ d_boundary
            # Photon escapes
            total_path += d_boundary
            if cluster_started
                # Check FV before classifying
                if in_fv(x_cluster, y_cluster, z_cluster, params)
                    outcome = classify_ss_energy(E_cluster, rng, params)
                else
                    outcome = :SS_outside_FV
                end
            else
                outcome = :escaped
            end
            break
        end

        # Move to interaction point
        x += d_int * dx;  y += d_int * dy;  z += d_int * dz
        total_path += d_int
        n_int += 1
        if n_int == 1
            first_int_x = x;  first_int_y = y;  first_int_z = z
        end

        # BFV check on first interaction
        if !bfv_checked
            bfv_checked = true
            if !in_bfv(x, y, z, geom.bfv)
                outcome = :outside_bfv
                break
            end
        end

        # Sample interaction type
        σ_pe = σ_photo(xcom, E)
        σ_co = σ_Compton(xcom, E)
        σ_pp = σ_pair(xcom, E)
        σ_tot = σ_pe + σ_co + σ_pp
        r_branch = rand(rng)

        if r_branch < σ_pp / σ_tot
            # Pair production → immediate MS rejection
            if has_buffer_space
                push!(traj_points, InteractionPoint(x, y, z, 0.0, :pair, E))
            end
            outcome = :MS_rejected
            break

        elseif r_branch < (σ_pp + σ_pe) / σ_tot
            # Photoelectric absorption
            E_dep = E
            if has_buffer_space
                push!(traj_points, InteractionPoint(x, y, z, 0.0, :photo, E_dep))
            end

            # Check deposit region
            r2_dep = x*x + y*y
            if r2_dep > R_skin_inner2
                # In skin: accumulate for skin veto
                E_skin_deposit += E_dep
                if E_skin_deposit ≥ E_skin_veto
                    outcome = :skin_vetoed;  break
                end
            elseif E_dep ≥ E_visible && z ≥ geom.z_cathode && z ≤ geom.L_lxe
                # In drift region: normal SS/MS logic
                if !cluster_started
                    z_cluster = z;  x_cluster = x;  y_cluster = y
                    E_cluster = E_dep;  cluster_started = true
                else
                    if abs(z - z_cluster) < Δz_thresh
                        E_cluster += E_dep
                    else
                        outcome = :MS_rejected;  break
                    end
                end
            end
            # Photon fully absorbed
            if outcome == :in_progress
                if cluster_started
                    if in_fv(x_cluster, y_cluster, z_cluster, params)
                        outcome = classify_ss_energy(E_cluster, rng, params)
                    else
                        outcome = :SS_outside_FV
                    end
                else
                    outcome = :escaped
                end
            end
            break

        else  # Compton
            E_scatt, cos_θ = sample_klein_nishina(rng, E)
            E_dep = E - E_scatt
            if has_buffer_space
                push!(traj_points, InteractionPoint(x, y, z, E_scatt, :compton, E_dep))
            end

            # Check deposit region
            r2_dep = x*x + y*y
            if r2_dep > R_skin_inner2
                E_skin_deposit += E_dep
                if E_skin_deposit ≥ E_skin_veto
                    outcome = :skin_vetoed;  break
                end
            elseif E_dep ≥ E_visible && z ≥ geom.z_cathode && z ≤ geom.L_lxe
                if !cluster_started
                    z_cluster = z;  x_cluster = x;  y_cluster = y
                    E_cluster = E_dep;  cluster_started = true
                else
                    if abs(z - z_cluster) < Δz_thresh
                        E_cluster += E_dep
                    else
                        outcome = :MS_rejected;  break
                    end
                end
            end

            # Update direction
            φ_az = 2π * rand(rng)
            dx, dy, dz = rotate_direction(dx, dy, dz, cos_θ, φ_az)
            E = E_scatt

            # Energy cutoff
            if E < E_cutoff
                E_dep_final = E
                r2_dep_final = x*x + y*y
                if r2_dep_final > R_skin_inner2
                    E_skin_deposit += E_dep_final
                    if E_skin_deposit ≥ E_skin_veto
                        outcome = :skin_vetoed;  break
                    end
                elseif E_dep_final ≥ E_visible && z ≥ geom.z_cathode && z ≤ geom.L_lxe
                    if !cluster_started
                        z_cluster = z;  x_cluster = x;  y_cluster = y
                        E_cluster = E_dep_final;  cluster_started = true
                    else
                        if abs(z - z_cluster) < Δz_thresh
                            E_cluster += E_dep_final
                        else
                            outcome = :MS_rejected;  break
                        end
                    end
                end
                if has_buffer_space
                    push!(traj_points, InteractionPoint(x, y, z, 0.0, :cutoff, E_dep_final))
                end
                if outcome == :in_progress
                    if cluster_started
                        if in_fv(x_cluster, y_cluster, z_cluster, params)
                            outcome = classify_ss_energy(E_cluster, rng, params)
                        else
                            outcome = :SS_outside_FV
                        end
                    else
                        outcome = :escaped
                    end
                end
                break
            end
        end
    end

    # 5. Companion gamma veto for ²⁰⁸Tl events
    # If the 2.6 MeV gamma passed all cuts (SS_in_ROI), check whether the
    # companion gamma (583/860/763 keV) is visible anywhere in the LXe.
    # If visible → event is vetoed (:companion_vetoed).
    if outcome == :SS_in_ROI && source.is_Tl208
        E_comp = sample_companion_energy(rng)
        if E_comp > 0.0
            # Companion emitted from same source point, random 4π direction
            x0 = traj_points[1].x
            y0 = traj_points[1].y
            z0 = traj_points[1].z
            if track_companion_gamma(rng, x0, y0, z0, E_comp, geom, xcom, params)
                outcome = :companion_vetoed
            end
        end
    end

    # 6. Record outcome
    state.counts[outcome] += 1

    # Fill histograms
    if n_int > 0
        _fill_2d!(state.H_first_interaction, first_int_z,
                  first_int_x^2 + first_int_y^2, params)
    end

    if cluster_started
        r2_cl = x_cluster^2 + y_cluster^2
        _fill_2d!(state.H_cluster_position, z_cluster, r2_cl, params)

        if outcome == :SS_in_ROI
            _fill_2d!(state.H_signal, z_cluster, r2_cl, params)
        end

        _fill_2d_weighted!(state.H_total_energy, z_cluster, r2_cl,
                           E_cluster * 1000.0, params)

        if outcome == :SS_in_ROI || outcome == :SS_outside_ROI
            iE = _bin_index(E_cluster * 1000.0, 0.0, params.hist_E_max_keV,
                            params.hist_n_E_bins)
            if iE > 0
                state.E_total_cluster_all_SS[iE] += 1.0
            end
        end
    end

    # n_interactions histogram
    if n_int > 0
        ib = min(n_int, length(state.n_interactions_per_photon))
        state.n_interactions_per_photon[ib] += 1.0
    end

    # Path length histogram (5 cm bins)
    if total_path > 0.0
        ib = _bin_index(total_path, 0.0, 500.0, length(state.path_length_in_LXe))
        if ib > 0
            state.path_length_in_LXe[ib] += 1.0
        end
    end

    # Trajectory commit
    if has_buffer_space && outcome != :in_progress
        traj = Trajectory(
            traj_points[1].x, traj_points[1].y, traj_points[1].z,
            source.label,
            dx, dy, dz,
            traj_points, outcome,
            cluster_started ? z_cluster : NaN,
            cluster_started ? x_cluster : NaN,
            cluster_started ? y_cluster : NaN,
            cluster_started ? E_cluster : 0.0,
            cluster_started, n_int
        )
        commit!(state.traj_buffer, traj)

        # Also commit to FV buffer if cluster is inside the fiducial volume
        if cluster_started && in_fv(x_cluster, y_cluster, z_cluster, params)
            commit!(state.traj_buffer_fv, traj)
        end
    end

    return outcome
end

"""
    run_mc(params::Params, source::SourceConfig) -> Result

Top-level MC driver. Builds geometry and XCOM table, partitions samples
across threads, runs the MC loop, merges results.
"""
function run_mc(params::Params, source::SourceConfig)::Result
    geom = build_geometry(params)
    xcom = load_xcom(params.phys_xcom_data_path)

    @printf("\n── Source: %s ──\n", source.label)
    @printf("  E_γ     = %.3f MeV\n", source.E_MeV)
    @printf("  Entry   = %s\n", source.entry)
    @printf("  Angular = %s\n", source.angular)
    if source.angular == :shaped
        @printf("  a=%.4f, b=%.4f, u_min=%.2f, norm=%.4f\n",
                source.a, source.b, source.u_min, source.norm)
    end
    @printf("  γ/yr    = %.3e\n", source.gammas_per_yr)
    @printf("──────────────────────\n\n")

    n_t = Threads.nthreads()
    samples_per_thread = div(params.mc_N_samples, n_t)
    remainder = params.mc_N_samples - samples_per_thread * n_t

    thread_results = [ThreadResult(params) for _ in 1:n_t]

    N = params.mc_N_samples
    report_every = max(1, N ÷ (n_t * 10))

    t_start = time()

    Threads.@threads for tid in 1:n_t
        n_local = samples_per_thread + (tid ≤ remainder ? 1 : 0)
        rng = MersenneTwister(params.mc_seed + tid)
        tr = thread_results[tid]
        for i in 1:n_local
            track_one_photon!(tr, rng, geom, source, xcom, params)
            if tid == 1 && i % report_every == 0
                est_done = i * n_t
                elapsed = time() - t_start
                @printf("  Progress: ~%.1f%%  (%d / %d)  elapsed: %.1f s\r",
                        100.0 * est_done / N, est_done, N, elapsed)
            end
        end
    end
    @printf("  Progress: 100.0%%  (%d / %d)  elapsed: %.1f s\n",
            N, N, time() - t_start)

    runtime = time() - t_start
    merge_results(thread_results, geom, source, xcom, params, runtime)
end

"""
    merge_results(...) -> Result

Sum histograms element-wise, sum counts, merge trajectory buffers.
"""
function merge_results(thread_results::Vector{ThreadResult},
                       geom::Geometry, source::SourceConfig,
                       xcom::XCOMTable, params::Params,
                       runtime_seconds::Float64)::Result
    ref = thread_results[1]

    counts = Dict{Symbol, Int}()
    for k in keys(ref.counts)
        counts[k] = sum(tr.counts[k] for tr in thread_results)
    end

    H_fi  = sum(tr.H_first_interaction for tr in thread_results)
    H_cp  = sum(tr.H_cluster_position  for tr in thread_results)
    H_sig = sum(tr.H_signal            for tr in thread_results)
    H_te  = sum(tr.H_total_energy      for tr in thread_results)
    E_ss  = sum(tr.E_total_cluster_all_SS for tr in thread_results)
    n_ip  = sum(tr.n_interactions_per_photon for tr in thread_results)
    pl    = sum(tr.path_length_in_LXe for tr in thread_results)

    bin_edges_z   = range(params.hist_z_min_cm,   params.hist_z_max_cm,
                          length=params.hist_n_z_bins + 1) |> collect
    bin_edges_r2  = range(params.hist_r2_min_cm2,  params.hist_r2_max_cm2,
                          length=params.hist_n_r2_bins + 1) |> collect
    bin_edges_E   = range(0.0, params.hist_E_max_keV,
                          length=params.hist_n_E_bins + 1) |> collect

    trajs = merge_buffers([tr.traj_buffer for tr in thread_results],
                          params.mc_n_traj_per_outcome)
    trajs_fv = merge_buffers([tr.traj_buffer_fv for tr in thread_results],
                             params.mc_n_traj_per_outcome)

    Result(params, geom, source, xcom, counts,
           H_fi, H_cp, H_sig, H_te, E_ss, n_ip, pl,
           bin_edges_z, bin_edges_r2, bin_edges_E,
           trajs, trajs_fv, runtime_seconds)
end
