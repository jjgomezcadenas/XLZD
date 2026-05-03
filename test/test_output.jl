@testset "integrate_in_FV_box" begin
    # Use the new bbonu FV defaults: z ∈ [26, 96], r² ≤ 1521 (r ≤ 39)
    # 4×4 histogram with bin edges chosen for easy reasoning.
    # z edges: [0, 30, 60, 120, 300] → centers: 15, 45, 90, 210
    # r² edges: [0, 500, 1000, 2000, 5000] → centers: 250, 750, 1500, 3500
    # FV: z ∈ [26, 96], r² ≤ 1521
    # Bins inside FV:
    #   z_center=45 (yes), z_center=90 (yes), z_center=15 (no), z_center=210 (no)
    #   r²_center=250 (yes), r²_center=750 (yes), r²_center=1500 (yes), r²_center=3500 (no)
    # So 2 × 3 = 6 bins inside FV

    params = Params(
        hist_n_z_bins=4, hist_n_r2_bins=4,
        hist_z_min_cm=0.0, hist_z_max_cm=300.0,
        hist_r2_min_cm2=0.0, hist_r2_max_cm2=5000.0
    )

    bin_edges_z  = [0.0, 30.0, 60.0, 120.0, 300.0]
    bin_edges_r2 = [0.0, 500.0, 1000.0, 2000.0, 5000.0]

    # Fill all bins with 1.0
    H = ones(Float64, 4, 4)
    n = XLZD.integrate_in_FV_box(H, bin_edges_z, bin_edges_r2, params)
    @test n == 6

    # Specific bins
    H2 = zeros(Float64, 4, 4)
    H2[2, 1] = 10.0  # z=45,  r²=250  → inside
    H2[3, 2] = 10.0  # z=90,  r²=750  → inside
    H2[1, 1] = 99.0  # z=15   → outside (z < 26)
    H2[4, 4] = 99.0  # z=210, r²=3500 → outside
    n2 = XLZD.integrate_in_FV_box(H2, bin_edges_z, bin_edges_r2, params)
    @test n2 == 20  # 10 + 10
end

@testset "to_json — structure" begin
    params = Params(mc_N_samples=1000, mc_n_traj_per_outcome=2)
    source = SourceConfig(
        label="test_json", E_MeV=2.448, entry=:barrel, angular=:flat,
        R_entry=74.3, H_entry=159.35, z_min_entry=70.0,
        gammas_per_yr=1.0e5,
    )
    result = run_mc(params, source)

    json_str = XLZD.to_json(result.trajectories, result.trajectories_fv,
                            result.geometry, result.params)
    data = JSON3.read(json_str)

    # Top-level keys
    @test haskey(data, :geometry)
    @test haskey(data, :trajectories)
    @test haskey(data, :trajectories_fv)

    # Geometry metadata matches
    @test data.geometry.R_lxe == result.geometry.R_lxe

    # Trajectories have correct structure
    if length(data.trajectories) > 0
        t = data.trajectories[1]
        @test haskey(t, :outcome)
        @test haskey(t, :entry)
        @test haskey(t, :interactions)
        @test haskey(t, :n_interactions)
    end
end

@testset "write_outputs — smoke test" begin
    test_dir = mktempdir()
    params = Params(
        mc_N_samples=1000,
        mc_n_traj_per_outcome=2,
        out_dir=test_dir,
    )
    source = SourceConfig(
        label="test_smoke", E_MeV=2.448, entry=:barrel, angular=:flat,
        R_entry=74.3, H_entry=159.35, z_min_entry=70.0,
        gammas_per_yr=1.0e5,
    )

    result = run_mc(params, source)
    write_outputs(result)

    # Check expected files exist
    @test isfile(joinpath(test_dir, "results.h5"))
    @test isfile(joinpath(test_dir, "summary.txt"))
    @test isfile(joinpath(test_dir, "trajectories.json"))
    @test isfile(joinpath(test_dir, "cut1_dNdu.png"))
    @test isfile(joinpath(test_dir, "cut2_cluster_pre_fv.png"))
    @test isfile(joinpath(test_dir, "cut3_ssms.png"))
    @test isfile(joinpath(test_dir, "cut4_ss_fv.png"))
    @test isfile(joinpath(test_dir, "cut4_energy.png"))
    @test isfile(joinpath(test_dir, "cut4_signal.png"))
    @test isfile(joinpath(test_dir, "diagnostics.png"))

    # summary.txt should contain key info
    txt = read(joinpath(test_dir, "summary.txt"), String)
    @test occursin("Outcome counts", txt)
    @test occursin("SS_in_ROI", txt)
    @test occursin("companion_vetoed", txt)
    @test occursin("test_smoke", txt)
end
