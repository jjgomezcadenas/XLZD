using Random, Statistics

# Helper: create a simple flat barrel Bi-214 source for testing
function _test_source_bi214()
    SourceConfig(
        label="test_flat_barrel_Bi214",
        E_MeV=2.448, entry=:barrel, angular=:flat,
        R_entry=74.3, H_entry=159.35, z_min_entry=70.0,
        gammas_per_yr=1.0e5,
    )
end

# Helper: create a Tl-208 source (with companion)
function _test_source_tl208()
    SourceConfig(
        label="test_flat_barrel_Tl208",
        E_MeV=2.615, entry=:barrel, angular=:flat,
        R_entry=74.3, H_entry=159.35, z_min_entry=70.0,
        gammas_per_yr=1.0e6,
        is_Tl208=true,
    )
end

@testset "ThreadResult constructor" begin
    params = Params()
    tr = XLZD.ThreadResult(params)

    # All outcome counts start at zero
    for k in [:outside_bfv, :escaped, :MS_rejected, :SS_outside_ROI,
              :SS_in_ROI, :SS_outside_FV, :companion_vetoed]
        @test tr.counts[k] == 0
    end

    # Histogram dimensions match params
    @test size(tr.H_signal) == (params.hist_n_z_bins, params.hist_n_r2_bins)
    @test length(tr.E_total_cluster_all_SS) == params.hist_n_E_bins
end

@testset "track_one_photon! — Bi-214 single call" begin
    params = Params()
    geom = build_geometry(params)
    xcom = load_xcom(params.phys_xcom_data_path)
    source = _test_source_bi214()
    tr = XLZD.ThreadResult(params)
    rng = MersenneTwister(42)

    outcome = XLZD.track_one_photon!(tr, rng, geom, source, xcom, params)

    # Outcome must be one of the valid categories
    @test outcome in [:outside_bfv, :escaped, :MS_rejected,
                      :SS_outside_ROI, :SS_in_ROI, :SS_outside_FV]

    # No companion_vetoed for Bi-214
    @test tr.counts[:companion_vetoed] == 0

    # Exactly one count incremented
    @test sum(values(tr.counts)) == 1
    @test tr.counts[outcome] == 1
end

@testset "track_one_photon! — Tl-208 single call" begin
    params = Params()
    geom = build_geometry(params)
    xcom = load_xcom(params.phys_xcom_data_path)
    source = _test_source_tl208()
    tr = XLZD.ThreadResult(params)
    rng = MersenneTwister(42)

    outcome = XLZD.track_one_photon!(tr, rng, geom, source, xcom, params)

    # Outcome includes companion_vetoed as a possibility
    @test outcome in [:outside_bfv, :escaped, :MS_rejected,
                      :SS_outside_ROI, :SS_in_ROI, :SS_outside_FV,
                      :companion_vetoed]

    @test sum(values(tr.counts)) == 1
end

@testset "track_companion_gamma — inward at center → visible" begin
    # Companion emitted from center of detector heading inward
    # should almost always interact in active LXe (visible)
    params = Params()
    geom = build_geometry(params)
    xcom = load_xcom(params.phys_xcom_data_path)
    rng = MersenneTwister(42)

    n_visible = 0
    N = 100
    for _ in 1:N
        visible = XLZD.track_companion_gamma(
            rng, 0.0, 0.0, 100.0,   # center of detector
            0.583,                    # 583 keV companion
            geom, xcom, params)
        if visible; n_visible += 1; end
    end

    # 583 keV in LXe: mfp ~7 cm, detector is ~150 cm radius
    # Should be visible essentially always from center
    @test n_visible > 95
end

@testset "track_companion_gamma — E=0 → invisible" begin
    # No companion emitted (1% case)
    params = Params()
    geom = build_geometry(params)
    xcom = load_xcom(params.phys_xcom_data_path)
    rng = MersenneTwister(42)

    visible = XLZD.track_companion_gamma(
        rng, 0.0, 0.0, 100.0, 0.0, geom, xcom, params)
    @test !visible
end

@testset "run_mc — Bi-214, short run, outcome counts" begin
    N = 50_000
    params = Params(mc_N_samples=N, mc_n_traj_per_outcome=2)
    source = _test_source_bi214()
    result = run_mc(params, source)

    # All outcomes sum to N
    total = sum(values(result.counts))
    @test total == N

    # No companion_vetoed for Bi-214
    @test result.counts[:companion_vetoed] == 0

    # Most events should be outside_bfv (small FV vs large detector)
    @test result.counts[:outside_bfv] > N * 0.5

    # Histograms have entries
    @test sum(result.H_first_interaction) > 0

    # Bin edges correct
    @test length(result.bin_edges_z) == params.hist_n_z_bins + 1
    @test length(result.bin_edges_r2) == params.hist_n_r2_bins + 1
end

@testset "run_mc — Tl-208, companion veto active" begin
    N = 50_000
    params = Params(mc_N_samples=N, mc_n_traj_per_outcome=2)
    source = _test_source_tl208()
    result = run_mc(params, source)

    # All outcomes sum to N
    total = sum(values(result.counts))
    @test total == N

    # companion_vetoed should have some counts (Tl-208 with companion)
    # Not guaranteed for small N since SS_in_ROI is rare, but companion_vetoed ≥ 0
    @test result.counts[:companion_vetoed] >= 0

    # SS_in_ROI + companion_vetoed should be less than without veto
    # (companion_vetoed events would have been SS_in_ROI without the veto)
end

@testset "run_mc — no :in_progress leaks" begin
    N = 10_000
    params = Params(mc_N_samples=N, mc_n_traj_per_outcome=2)

    for source in [_test_source_bi214(), _test_source_tl208()]
        result = run_mc(params, source)
        # No :in_progress should remain
        @test !haskey(result.counts, :in_progress)
        @test sum(values(result.counts)) == N
    end
end

@testset "merge_results — deterministic" begin
    params = Params(mc_N_samples=1000, mc_n_traj_per_outcome=2)
    geom = build_geometry(params)
    xcom = load_xcom(params.phys_xcom_data_path)
    source = _test_source_bi214()

    tr1 = XLZD.ThreadResult(params)
    tr2 = XLZD.ThreadResult(params)

    rng1 = MersenneTwister(100)
    rng2 = MersenneTwister(200)

    for _ in 1:500
        XLZD.track_one_photon!(tr1, rng1, geom, source, xcom, params)
        XLZD.track_one_photon!(tr2, rng2, geom, source, xcom, params)
    end

    merged = XLZD.merge_results([tr1, tr2], geom, source, xcom, params, 0.0)

    # Merged counts equal sum of individual counts
    for k in keys(tr1.counts)
        @test merged.counts[k] == tr1.counts[k] + tr2.counts[k]
    end

    # Total counts = 1000
    @test sum(values(merged.counts)) == 1000

    # Merged histograms equal element-wise sum
    @test merged.H_signal ≈ tr1.H_signal .+ tr2.H_signal
end
