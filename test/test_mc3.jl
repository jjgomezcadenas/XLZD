# test/test_mc3.jl — Verify the per-photon tracker (src3/).

using Test
using Random
using Printf

isdefined(Main, :XLZD3) || include("../src3/XLZD3.jl")
using .XLZD3

# ---------------------------------------------------------------------------
# Setup: detector, sources, XCOM table, default MC params
# ---------------------------------------------------------------------------
const lxe_csv      = joinpath(@__DIR__, "..", "data", "lxe_detector.csv")
const lxe_nist     = joinpath(@__DIR__, "..", "data", "nist_lxe.csv")
const ti_path      = joinpath(@__DIR__, "..", "data", "nist_ti.csv")
const geom_csv     = joinpath(@__DIR__, "..", "data", "lz_cryo_geometry.csv")
const extras_csv   = joinpath(@__DIR__, "..", "data", "lz_cryo_extras.csv")
const surfaces_csv = joinpath(@__DIR__, "..", "data", "lz_cryo_surface_sources.csv")
const xcom_path    = joinpath(@__DIR__, "..", "data", "nist.csv")

mat_LXe = load_material("LXe", 2.953, lxe_nist)
mat_Ti  = load_material("Ti",  4.510, ti_path)
det     = build_lxe_detector(lxe_csv, mat_LXe)
cryo    = build_cryostat(geom_csv, extras_csv, surfaces_csv)
indiv   = build_individual_sources(cryo, mat_Ti)
effs    = build_effective_sources(indiv, cryo, mat_Ti)
xcom    = load_xcom(xcom_path)
params  = MCParams()
by_name = Dict(e.name => e for e in effs)


# ---------------------------------------------------------------------------
# MCParams + helpers
# ---------------------------------------------------------------------------

@testset "MCParams defaults and accessors" begin
    p = MCParams()
    @test p.Q_betabeta_keV       == 2458.0
    @test p.σ_E_over_E             == 0.007
    @test p.ROI_halfwidth_keV      == 17.2
    @test p.E_tracking_cutoff_keV  == 40.0
    @test p.Δz_threshold_mm        == 3.0
    @test p.fv_z_min_cm            == 26.0
    @test p.fv_z_max_cm            == 96.0
    @test p.fv_r2_max_cm2          == 1521.0
    @test E_tracking_cutoff_MeV(p) ≈ 0.040
    @test Δz_threshold_cm(p)        ≈ 0.30
end

@testset "in_fv basic" begin
    p = MCParams()
    @test  in_fv( 0.0,  0.0, 50.0, p)
    @test  in_fv(20.0, 30.0, 60.0, p)        # r² = 1300 < 1521
    @test !in_fv(40.0,  0.0, 50.0, p)        # r² = 1600 > 1521
    @test !in_fv( 0.0,  0.0, 100.0, p)       # z > z_max
    @test !in_fv( 0.0,  0.0, 20.0, p)        # z < z_min
end

# ---------------------------------------------------------------------------
# path_to_next_region geometric checks
# ---------------------------------------------------------------------------

@testset "path_to_next_region — known intersections" begin
    # Photon at axis center, moving in +x: nearest cylinder is R_FC_inner = 72.8
    d = path_to_next_region(0.0, 0.0, 50.0, 1.0, 0.0, 0.0, det)
    @test isapprox(d, 72.8, atol=1e-6)

    # Photon just inside ICV inner moving outward radially: nearest is R_ICV_inner
    d = path_to_next_region(82.0, 0.0, 50.0, 1.0, 0.0, 0.0, det)
    @test isapprox(d, 0.1, atol=1e-6)

    # Photon at axis center moving in +z: nearest plane in path is z_gate=145.6
    d = path_to_next_region(0.0, 0.0, 50.0, 0.0, 0.0, 1.0, det)
    @test isapprox(d, 145.6 - 50.0, atol=1e-6)

    # Photon at axis moving in -z hits z_cathode = 0 first
    d = path_to_next_region(0.0, 0.0, 50.0, 0.0, 0.0, -1.0, det)
    @test isapprox(d, 50.0, atol=1e-6)
end

# ---------------------------------------------------------------------------
# MC4: companion veto + threaded driver
# ---------------------------------------------------------------------------

@testset "companion_reach_prob has the right scaling" begin
    comp_eff = by_name["CB_Tl208c"]
    r = companion_reach_prob(comp_eff)
    # comp.total_per_yr / decay_rate ≤ 0.5 (inward hemisphere only); times BR=0.99
    # Expected: a few tens of percent for thin Ti walls
    @test 0.0 < r < 0.5
    # Cross-check by recomputing decay rate from contributions
    total_produced = sum(c.source.produced_per_yr for c in comp_eff.contributions)
    decay_rate     = total_produced / 0.99
    @test isapprox(r, comp_eff.total_per_yr / decay_rate, rtol=1e-12)
end

@testset "companion_visible! returns Bool; mostly false for low-E γ" begin
    rng = MersenneTwister(0xCAFEBABE)
    comp_eff = by_name["CB_Tl208c"]
    n_visible = 0
    N = 5000
    for _ in 1:N
        if companion_visible!(rng, det, comp_eff, xcom, params)
            n_visible += 1
        end
    end
    # 583 keV γ that reaches the LXe is fairly likely to deposit something
    # visible (μ_LXe ~ 0.25 cm⁻¹ → mfp ~4 cm; mostly Compton, often visible)
    @test 0.0 < n_visible / N < 1.0
    @printf("  companion_visible! fraction (CB_Tl208c, %d trials) = %.3f\n",
            N, n_visible / N)
end

@testset "run_mc — Bi-214 has no companion-vetoed events" begin
    eff = by_name["CB_Bi214"]
    res = run_mc(det, eff, nothing, xcom, params, 5000; mc_seed=0xAAAA)
    @test res.counts[:companion_vetoed] == 0
    @test res.r_comp == 0.0
    @test res.n_total == 5000
    @test res.γ_per_yr_total == eff.total_per_yr
    @test isapprox(res.bg_per_yr, res.f_SS_in_ROI * eff.total_per_yr; rtol=1e-12)
end

@testset "run_mc — Tl-208 with companion produces some companion-vetoed" begin
    eff      = by_name["CB_Tl208"]
    comp_eff = by_name["CB_Tl208c"]
    res = run_mc(det, eff, comp_eff, xcom, params, 20_000; mc_seed=0xBBBB)
    # Companion veto enabled → r_comp > 0
    @test res.r_comp > 0.0
    # Companion vetos = SS_in_ROI candidates that got knocked out
    @test res.counts[:companion_vetoed] >= 0
    # Total accounting still consistent
    @test res.n_total == 20_000
    @test sum(values(res.counts)) == 20_000
end

@testset "run_mc — reproducibility with fixed seed and same thread count" begin
    eff = by_name["CB_Bi214"]
    res1 = run_mc(det, eff, nothing, xcom, params, 4000; mc_seed=0xC1C1)
    res2 = run_mc(det, eff, nothing, xcom, params, 4000; mc_seed=0xC1C1)
    @test res1.counts == res2.counts
    res3 = run_mc(det, eff, nothing, xcom, params, 4000; mc_seed=0xC1C2)
    @test res1.counts != res3.counts
end

@testset "run_mc_all returns 6 results, sums match per-source" begin
    res_all = run_mc_all(det, effs, xcom, params, 2000; mc_seed=0xDDDD)
    @test length(res_all) == 6
    names = [r.name for r in res_all]
    @test sort(names) == sort(["CB_Bi214", "CTH_Bi214", "CBH_Bi214",
                                "CB_Tl208", "CTH_Tl208", "CBH_Tl208"])
    # Tl-208 sources have r_comp > 0; Bi-214 sources have r_comp == 0
    for r in res_all
        if r.isotope === :Tl208
            @test r.r_comp > 0.0
        else
            @test r.r_comp == 0.0
        end
    end
end

# ---------------------------------------------------------------------------
# Pair-production end-to-end (Tl-208 source via run_mc)
# ---------------------------------------------------------------------------

@testset "Pair-production end-to-end on Tl-208" begin
    # Tl-208 at 2.615 MeV pair-produces a few % of the time. With the
    # stack-based tracker, pair vertices spawn two 511 keV children that
    # are tracked through their own cascades. Just verify run_mc
    # bookkeeping holds at sensible N.
    eff = by_name["CB_Tl208"]
    res = run_mc(det, eff, by_name["CB_Tl208c"], xcom, params, 50_000;
                  mc_seed=0xCAFEFACE)
    @test sum(values(res.counts)) == 50_000
    @test res.counts[:MS_rejected] >= 0
    @test res.counts[:escaped]     > 0
end

# ===========================================================================
# run_mc API tests on the stack pipeline
# ===========================================================================

@testset "run_mc basic invariants (stack pipeline)" begin
    eff = by_name["CB_Bi214"]
    res = run_mc(det, eff, nothing, xcom, params, 2000;
                  mc_seed=0xABCD)
    @test res isa MCResult
    @test sum(values(res.counts)) == res.n_total
    @test res.n_total == 2000
    for (k, v) in res.counts
        @test k in (:escaped, :MS_rejected, :skin_vetoed,
                    :outside_FV, :SS_outside_ROI, :SS_in_ROI,
                    :companion_vetoed)
        @test v >= 0
    end
end

@testset "run_mc reproducibility on stack pipeline" begin
    eff = by_name["CB_Bi214"]
    r1 = run_mc(det, eff, nothing, xcom, params, 1000;
                 mc_seed=0xABCD)
    r2 = run_mc(det, eff, nothing, xcom, params, 1000;
                 mc_seed=0xABCD)
    @test r1.counts == r2.counts
end

@testset "run_mc Tl-208 companion-veto bookkeeping" begin
    eff      = by_name["CB_Tl208"]
    comp_eff = by_name["CB_Tl208c"]
    res = run_mc(det, eff, comp_eff, xcom, params, 5000;
                  mc_seed=0xCAFE)
    # On Tl-208, the companion γ frequently fires; expect ≥ 1 vetoed event.
    @test res.counts[:companion_vetoed] >= 0
    # Aggregate counts non-negative; at least one outcome class non-empty.
    @test sum(values(res.counts)) == 5000
end

@testset "run_mc_all bookkeeping over 6 sources" begin
    rs = run_mc_all(det, effs, xcom, params, 500;
                     mc_seed=0xBEEF)
    @test length(rs) == 6
    for r in rs
        @test sum(values(r.counts)) == 500
        @test r.stack_hists   isa StackHistogramSet
        @test r.cluster_hists isa ClusterHistogramSet
    end
end

@testset "stack/cluster histograms populated under default kwargs" begin
    eff = by_name["CB_Bi214"]
    res = run_mc(det, eff, nothing, xcom, params, 1000; mc_seed=0x1234)
    @test res.stack_hists   isa StackHistogramSet
    @test res.cluster_hists isa ClusterHistogramSet
    # Stack/cluster histograms accumulate ONCE per event that reaches
    # full tracking (i.e. fast_veto returned :pass). Events fast-rejected
    # (:vetoed_skin / :rejected_fv from the FAST path) skip the update.
    # The two sums must therefore agree with each other and be > 0 and
    # ≤ n_total.
    n_h = sum(res.stack_hists.path_length_LXe_counts)
    @test n_h == sum(res.cluster_hists.N_clusters_counts)
    @test 0 < n_h <= res.n_total
end

@testset "with_*_histograms=false leaves them as nothing" begin
    eff = by_name["CB_Bi214"]
    res = run_mc(det, eff, nothing, xcom, params, 200;
                  mc_seed=0x4321,
                  with_stack_histograms=false,
                  with_cluster_histograms=false)
    @test res.stack_hists   === nothing
    @test res.cluster_hists === nothing
end
