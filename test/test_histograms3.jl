# test/test_histograms3.jl — Tests for the rejection histograms.
#
# The legacy 6-bin HistogramSet (along with LXeDeposit, PhotonScratch,
# update_histograms!, merge_histograms!) was deleted in the histogram
# refactor. Stack/Cluster diagnostic histograms now live in their own
# test files (test_stack_histograms3.jl, test_cluster_histograms3.jl).
# This file covers RejectionHistograms only.

using Test

isdefined(Main, :XLZD3) || include("../src3/XLZD3.jl")
using .XLZD3

# ---------------------------------------------------------------------------
# RejectionHistograms shape and accumulation
# ---------------------------------------------------------------------------

@testset "RejectionHistograms defaults" begin
    rh = RejectionHistograms()
    @test size(rh.skin_r2z_counts) == (60, 100)
    @test size(rh.fv_r2z_counts)   == (60, 100)
    @test length(rh.skin_E_counts) == 270
    @test length(rh.fv_E_counts)   == 270
    @test all(rh.skin_r2z_counts  .== 0)
    @test all(rh.skin_E_counts    .== 0)
    @test all(rh.fv_r2z_counts    .== 0)
    @test all(rh.fv_E_counts      .== 0)
end

@testset "fill_rejected_skin! increments r²×z and E bins" begin
    rh = RejectionHistograms()
    fill_rejected_skin!(rh, 50.0, 0.0, 50.0, 0.5)
    @test sum(rh.skin_r2z_counts) == 1
    @test sum(rh.skin_E_counts)   == 1
    @test sum(rh.fv_r2z_counts)   == 0
    @test sum(rh.fv_E_counts)     == 0
end

@testset "fill_rejected_fv! increments fv-side bins only" begin
    rh = RejectionHistograms()
    fill_rejected_fv!(rh, 0.0, 0.0, 50.0, 0.020)
    @test sum(rh.fv_r2z_counts) == 1
    @test sum(rh.fv_E_counts)   == 1
    @test sum(rh.skin_r2z_counts) == 0
    @test sum(rh.skin_E_counts)   == 0
end

@testset "merge_rejection_histograms! sums elementwise" begin
    rh1 = RejectionHistograms()
    rh2 = RejectionHistograms()
    fill_rejected_skin!(rh1, 50.0, 0.0,  50.0, 0.5)
    fill_rejected_fv!(  rh2,  0.0, 0.0,  50.0, 0.020)
    merge_rejection_histograms!(rh1, rh2)
    @test sum(rh1.skin_r2z_counts) == 1
    @test sum(rh1.skin_E_counts)   == 1
    @test sum(rh1.fv_r2z_counts)   == 1
    @test sum(rh1.fv_E_counts)     == 1
end

println("\n  ── test_histograms3.jl: RejectionHistograms OK ──\n")
