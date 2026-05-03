# test/test_cluster_histograms3.jl — Unit tests for ClusterHistogramSet.
#
# Stripped to the single field that survives the cut-flow refactor: the
# per-cluster Ec spectrum. SS pre-ROI fills moved to CutHistograms
# (tested in test_cut_histograms3.jl).

using Test

isdefined(Main, :XLZD3) || include("../src3/XLZD3.jl")
using .XLZD3

const params = MCParams()

mk(; x=0.0, y=0.0, z=50.0, ec=1.0, es=1.0) = Cluster(x, y, z, ec, es)

# ---------------------------------------------------------------------------
# 1. Empty cluster vector: no Ec entries.
# ---------------------------------------------------------------------------

@testset "1. empty clusters → no Ec fill" begin
    ch = ClusterHistogramSet()
    update_cluster_histograms!(ch, Cluster[], params)
    @test sum(ch.Ec_counts) == 0
end

# ---------------------------------------------------------------------------
# 2. Single cluster: Ec[bin] = 1
# ---------------------------------------------------------------------------

@testset "2. single cluster bins Ec" begin
    ch = ClusterHistogramSet()
    update_cluster_histograms!(ch, [mk(ec=2.448, es=2.448, z=50.0)], params)
    bin_E = floor(Int, 2.448 / 2.7 * 270) + 1
    @test ch.Ec_counts[bin_E] == 1
    @test sum(ch.Ec_counts)   == 1
end

# ---------------------------------------------------------------------------
# 3. Multiple clusters: one Ec entry per cluster.
# ---------------------------------------------------------------------------

@testset "3. one Ec entry per cluster" begin
    ch = ClusterHistogramSet()
    cs = [mk(ec=0.3), mk(ec=0.5), mk(ec=1.0)]
    update_cluster_histograms!(ch, cs, params)
    bin = E -> floor(Int, E / 2.7 * 270) + 1
    @test ch.Ec_counts[bin(0.3)] == 1
    @test ch.Ec_counts[bin(0.5)] == 1
    @test ch.Ec_counts[bin(1.0)] == 1
    @test sum(ch.Ec_counts) == 3
end

# ---------------------------------------------------------------------------
# 4. merge_cluster_histograms!
# ---------------------------------------------------------------------------

@testset "4. merge_cluster_histograms! sums Ec counts" begin
    ch1 = ClusterHistogramSet()
    ch2 = ClusterHistogramSet()
    update_cluster_histograms!(ch1, [mk(ec=1.0)], params)
    update_cluster_histograms!(ch2, [mk(ec=2.0)], params)
    merge_cluster_histograms!(ch1, ch2)
    bin = E -> floor(Int, E / 2.7 * 270) + 1
    @test ch1.Ec_counts[bin(1.0)] == 1
    @test ch1.Ec_counts[bin(2.0)] == 1
    @test sum(ch1.Ec_counts)      == 2
end

@testset "5. merge_cluster_histograms! asserts compatible binning" begin
    ch1 = ClusterHistogramSet(E_n_bins=270)
    ch2 = ClusterHistogramSet(E_n_bins=100)
    @test_throws AssertionError merge_cluster_histograms!(ch1, ch2)
end

println("\n  ── test_cluster_histograms3.jl: ClusterHistogramSet OK ──\n")
