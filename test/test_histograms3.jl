# test/test_histograms3.jl — Verify the control-histogram pipeline.

using Test
using Random
using Printf

isdefined(Main, :XLZD3) || include("../src3/XLZD3.jl")
using .XLZD3

# ---------------------------------------------------------------------------
# HistogramSet basic shape
# ---------------------------------------------------------------------------

@testset "HistogramSet defaults" begin
    h = HistogramSet()
    @test length(h.ssms_counts)        == 3
    @test length(h.Δz_counts)          == 100
    @test length(h.E_first_counts)     == 270
    @test length(h.E_cluster_counts)   == 270
    @test length(h.N_clusters_counts)  == 21
    @test length(h.N_extra_counts)     == 21
    @test all(h.ssms_counts        .== 0)
    @test all(h.Δz_counts          .== 0)
    @test all(h.E_first_counts     .== 0)
    @test all(h.E_cluster_counts   .== 0)
    @test all(h.N_clusters_counts  .== 0)
    @test all(h.N_extra_counts     .== 0)
end

# ---------------------------------------------------------------------------
# update_histograms! on hand-built scratch buffers
# ---------------------------------------------------------------------------

@testset "no deposits → no_cluster bucket only" begin
    h = HistogramSet()
    s = PhotonScratch()
    p = MCParams()
    update_histograms!(h, s, p)
    @test h.ssms_counts == [0, 0, 1]              # no_cluster
    @test h.N_clusters_counts[1] == 1             # bin 1 = 0 clusters
    @test h.N_extra_counts[1]    == 1
    @test sum(h.E_first_counts)   == 0
    @test sum(h.E_cluster_counts) == 0
    @test sum(h.Δz_counts)        == 0
end

@testset "single deposit → SS, 1 cluster, no Δz fill" begin
    h = HistogramSet()
    s = PhotonScratch()
    p = MCParams()
    push!(s.deposits, LXeDeposit(0.0, 0.0, 5.0, 0.5, :TPC))
    update_histograms!(h, s, p)

    @test h.ssms_counts == [1, 0, 0]              # SS
    @test h.N_clusters_counts[2] == 1             # bin 2 = 1 cluster
    @test h.N_extra_counts[1]    == 1             # 0 extra
    # E_first = 0.5 MeV → bin (0.5/2.7*270)+1 = 51
    @test h.E_first_counts[51]   == 1
    @test sum(h.E_first_counts)  == 1
    # one cluster of energy 0.5 MeV → same bin
    @test h.E_cluster_counts[51] == 1
    @test sum(h.Δz_counts)       == 0             # only 1 deposit, no Δz
end

@testset "two deposits within Δz < 3mm → SS, 1 cluster summed" begin
    h = HistogramSet()
    s = PhotonScratch()
    p = MCParams()
    push!(s.deposits, LXeDeposit(0.0, 0.0, 5.00, 0.5, :TPC))
    push!(s.deposits, LXeDeposit(0.0, 0.0, 5.10, 0.3, :TPC))   # 1 mm away
    update_histograms!(h, s, p)

    @test h.ssms_counts == [1, 0, 0]              # SS
    @test h.N_clusters_counts[2] == 1             # 1 cluster
    @test h.N_extra_counts[1]    == 1             # 0 extra
    # E_first = 0.5 MeV → bin 51
    @test h.E_first_counts[51]   == 1
    # cluster sum = 0.8 MeV → bin (0.8/2.7*270)+1 = 81
    @test h.E_cluster_counts[81] == 1
    # Δz = |5.1 − 5.0| = 0.1 cm → bin 1
    @test h.Δz_counts[1]         == 1
end

@testset "two deposits with Δz > 3mm → MS, 2 clusters" begin
    h = HistogramSet()
    s = PhotonScratch()
    p = MCParams()
    push!(s.deposits, LXeDeposit(0.0, 0.0,  5.0, 0.5, :TPC))
    push!(s.deposits, LXeDeposit(0.0, 0.0, 10.0, 0.7, :TPC))   # 5 cm away
    update_histograms!(h, s, p)

    @test h.ssms_counts == [0, 1, 0]              # MS
    @test h.N_clusters_counts[3] == 1             # 2 clusters → bin 3
    @test h.N_extra_counts[2]    == 1             # 1 extra
    # Δz = 5 cm → bin (5/50*100)+1 = 11
    @test h.Δz_counts[11]        == 1
    # cluster 1: 0.5 MeV → bin 51
    # cluster 2: 0.7 MeV → bin (0.7/2.7*270)+1 = 71
    @test h.E_cluster_counts[51] == 1
    @test h.E_cluster_counts[71] == 1
end

@testset "merge_histograms! sums counts elementwise" begin
    h1 = HistogramSet()
    h2 = HistogramSet()
    s  = PhotonScratch()
    p  = MCParams()

    push!(s.deposits, LXeDeposit(0.0, 0.0,  5.0, 0.5, :TPC));  update_histograms!(h1, s, p)
    empty!(s.deposits)
    push!(s.deposits, LXeDeposit(0.0, 0.0, 10.0, 0.7, :TPC)); update_histograms!(h2, s, p)

    merge_histograms!(h1, h2)
    @test h1.ssms_counts == [2, 0, 0]
    @test sum(h1.N_clusters_counts) == 2
end

# Note: integration tests of run_mc + control HistogramSet were removed
# in step 11i. The legacy tracker — the only consumer that populated
# PhotonScratch.deposits — was deleted; control histograms will be
# re-sourced from the PhotonStack in a later step. The unit tests above
# (on hand-built PhotonScratch buffers) still exercise the histogram
# accumulation logic itself.
