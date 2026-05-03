# test/test_cut_histograms3.jl — Unit tests for CutHistograms.

using Test

isdefined(Main, :XLZD3) || include("../src3/XLZD3.jl")
using .XLZD3

const params = MCParams()
const E_VIS_MEV = params.E_visible_keV / 1000.0

mk(; x=0.0, y=0.0, z=50.0, ec=1.0, es=1.0) = Cluster(x, y, z, ec, es)

# ---------------------------------------------------------------------------
# Construction defaults
# ---------------------------------------------------------------------------

@testset "1. CutHistograms() defaults are all zero" begin
    ch = CutHistograms()
    @test sum(ch.h_u_sampled)              == 0
    @test sum(ch.first_interaction_r_z)    == 0
    @test sum(ch.dz_inclusive)             == 0
    @test sum(ch.n_visible)                == 0
    @test sum(ch.E_total)                  == 0
    @test sum(ch.ss_ec_pre_roi)            == 0
    @test sum(ch.ss_es_pre_roi)            == 0
    @test sum(ch.ss_r_z)                   == 0
end

# ---------------------------------------------------------------------------
# Cut 1 — fill_cut1_u!
# ---------------------------------------------------------------------------

@testset "2. fill_cut1_u! bins u into [0,1] correctly" begin
    ch = CutHistograms(u_n_bins=10)
    fill_cut1_u!(ch, 0.05)   # bin 1
    fill_cut1_u!(ch, 0.15)   # bin 2
    fill_cut1_u!(ch, 0.95)   # bin 10
    @test ch.h_u_sampled[1]  == 1
    @test ch.h_u_sampled[2]  == 1
    @test ch.h_u_sampled[10] == 1
    @test sum(ch.h_u_sampled) == 3
end

@testset "3. fill_cut1_u! out-of-range u skipped silently" begin
    ch = CutHistograms(u_n_bins=10)
    fill_cut1_u!(ch, -0.1)
    fill_cut1_u!(ch,  1.1)
    @test sum(ch.h_u_sampled) == 0
end

# ---------------------------------------------------------------------------
# Cut 2 — fill_cut2_first_interaction!
# ---------------------------------------------------------------------------

@testset "4. fill_cut2_first_interaction! bins (r, z) correctly" begin
    ch = CutHistograms(r_n_bins=10, r_max_cm=100.0,
                        z_n_bins=10, z_min_cm=0.0, z_max_cm=100.0)
    fill_cut2_first_interaction!(ch, 30.0, 40.0, 50.0)   # r = 50, z = 50 → bin (6, 6)
    @test ch.first_interaction_r_z[6, 6] == 1
    @test sum(ch.first_interaction_r_z)  == 1
end

@testset "5. fill_cut2_first_interaction! out-of-range silently dropped" begin
    ch = CutHistograms(r_n_bins=10, r_max_cm=10.0,
                        z_n_bins=10, z_min_cm=0.0, z_max_cm=10.0)
    fill_cut2_first_interaction!(ch, 100.0, 0.0, 5.0)    # r = 100 ≫ r_max
    @test sum(ch.first_interaction_r_z) == 0
end

# ---------------------------------------------------------------------------
# Cut 3 — fill_cut3!
# ---------------------------------------------------------------------------

@testset "6. fill_cut3!: 1 visible cluster — n_visible=1, no Δz, E_total=ec" begin
    ch = CutHistograms()
    fill_cut3!(ch, [mk(z=50.0, ec=1.0)], params)
    @test ch.n_visible[2]            == 1                 # bin 2 = 1 visible
    @test sum(ch.dz_inclusive)        == 0                 # no pairs
    bin = E -> floor(Int, E / 2.7 * 270) + 1
    @test ch.E_total[bin(1.0)]        == 1
end

@testset "7. fill_cut3!: 3 visible clusters z-spaced — 2 dz entries" begin
    ch = CutHistograms()
    cs = [mk(z=10.0, ec=0.5), mk(z=10.5, ec=0.4), mk(z=12.0, ec=0.6)]
    fill_cut3!(ch, cs, params)
    @test ch.n_visible[4] == 1                             # bin 4 = 3 visible
    # |Δz| values: 0.5 cm and 1.5 cm. Bins: dz_max = 5 cm, 200 bins → bin width = 0.025 cm
    bin_dz = d -> floor(Int, d / 5.0 * 200) + 1
    @test ch.dz_inclusive[bin_dz(0.5)] == 1
    @test ch.dz_inclusive[bin_dz(1.5)] == 1
    @test sum(ch.dz_inclusive) == 2
    bin_E = E -> floor(Int, E / 2.7 * 270) + 1
    @test ch.E_total[bin_E(0.5 + 0.4 + 0.6)] == 1
end

@testset "8. fill_cut3!: sub-visible clusters skipped" begin
    ch = CutHistograms()
    sub = E_VIS_MEV * 0.5
    cs = [mk(z=10.0, ec=sub), mk(z=20.0, ec=1.0), mk(z=30.0, ec=sub)]
    fill_cut3!(ch, cs, params)
    @test ch.n_visible[2] == 1                             # only the middle cluster counts
    @test sum(ch.dz_inclusive) == 0                        # 1 visible → no pairs
    bin_E = E -> floor(Int, E / 2.7 * 270) + 1
    @test ch.E_total[bin_E(1.0)] == 1
end

@testset "9. fill_cut3!: empty cluster vector — n_visible[0+1]=1, no other entries" begin
    ch = CutHistograms()
    fill_cut3!(ch, Cluster[], params)
    @test ch.n_visible[1] == 1                             # bin 1 = 0 visible
    @test sum(ch.dz_inclusive) == 0
    @test sum(ch.E_total) == 0
end

# ---------------------------------------------------------------------------
# Cut 4 — fill_cut4!
# ---------------------------------------------------------------------------

@testset "10. fill_cut4!: visible cluster ec/es/r-z all binned" begin
    ch = CutHistograms()
    cs = [mk(x=10.0, y=0.0, z=50.0, ec=2.448, es=2.460)]
    fill_cut4!(ch, cs, params)
    bin_E = E -> floor(Int, E / 2.7 * 270) + 1
    @test ch.ss_ec_pre_roi[bin_E(2.448)] == 1
    @test ch.ss_es_pre_roi[bin_E(2.460)] == 1
    @test sum(ch.ss_r_z) == 1
end

@testset "11. fill_cut4!: picks visible cluster, not clusters[1]" begin
    ch = CutHistograms()
    sub = E_VIS_MEV * 0.5
    cs = [mk(z=10.0, ec=sub, es=sub),                       # sub-vis (would mis-bin if used)
          mk(x=5.0, y=0.0, z=50.0, ec=2.448, es=2.460)]    # visible
    fill_cut4!(ch, cs, params)
    bin_E = E -> floor(Int, E / 2.7 * 270) + 1
    @test ch.ss_ec_pre_roi[bin_E(2.448)] == 1
    @test ch.ss_es_pre_roi[bin_E(2.460)] == 1
end

@testset "12. fill_cut4!: no visible cluster → no fill" begin
    ch = CutHistograms()
    sub = E_VIS_MEV * 0.5
    cs = [mk(ec=sub, es=sub)]
    fill_cut4!(ch, cs, params)
    @test sum(ch.ss_ec_pre_roi) == 0
    @test sum(ch.ss_es_pre_roi) == 0
    @test sum(ch.ss_r_z)        == 0
end

# ---------------------------------------------------------------------------
# Merge
# ---------------------------------------------------------------------------

@testset "13. merge_cut_histograms! sums every field element-wise" begin
    ch1 = CutHistograms()
    ch2 = CutHistograms()
    fill_cut1_u!(ch1, 0.5)
    fill_cut1_u!(ch2, 0.5)
    fill_cut2_first_interaction!(ch1, 0.0, 0.0, 50.0)
    fill_cut2_first_interaction!(ch2, 0.0, 0.0, 50.0)
    fill_cut3!(ch1, [mk(z=10.0, ec=1.0)], params)
    fill_cut3!(ch2, [mk(z=10.0, ec=1.0)], params)
    fill_cut4!(ch1, [mk(x=5.0, z=50.0, ec=2.448, es=2.448)], params)
    fill_cut4!(ch2, [mk(x=5.0, z=50.0, ec=2.448, es=2.448)], params)

    merge_cut_histograms!(ch1, ch2)
    @test sum(ch1.h_u_sampled)              == 2
    @test sum(ch1.first_interaction_r_z)    == 2
    @test sum(ch1.n_visible)                == 2
    @test sum(ch1.E_total)                  == 2
    @test sum(ch1.ss_ec_pre_roi)            == 2
    @test sum(ch1.ss_es_pre_roi)            == 2
    @test sum(ch1.ss_r_z)                   == 2
end

@testset "14. merge_cut_histograms! asserts compatible binning" begin
    ch1 = CutHistograms(u_n_bins=100)
    ch2 = CutHistograms(u_n_bins=50)            # mismatched
    @test_throws AssertionError merge_cut_histograms!(ch1, ch2)
end

println("\n  ── test_cut_histograms3.jl: CutHistograms OK ──\n")
