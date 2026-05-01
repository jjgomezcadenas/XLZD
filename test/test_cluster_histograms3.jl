# test/test_cluster_histograms3.jl — Unit tests for ClusterHistogramSet.

using Test

isdefined(Main, :XLZD3) || include("../src3/XLZD3.jl")
using .XLZD3

const params = MCParams()

mk(; x=0.0, y=0.0, z=50.0, ec=1.0, es=1.0) = Cluster(x, y, z, ec, es)

# ---------------------------------------------------------------------------
# 1. Empty cluster vector: only N_clusters_counts[1] increments
# ---------------------------------------------------------------------------

@testset "1. empty clusters" begin
    ch = ClusterHistogramSet()
    update_cluster_histograms!(ch, Cluster[], params)
    @test ch.N_clusters_counts[1] == 1
    @test sum(ch.Ec_counts)        == 0
    @test sum(ch.Emax_counts)      == 0
    @test sum(ch.Einc_counts)      == 0
    @test sum(ch.closest_D3_counts) == 0
    @test sum(ch.r2_vs_z_2d_counts) == 0
end

# ---------------------------------------------------------------------------
# 2. Single cluster: Ec[bin] = 1; Emax/Emin/Einc all in same bin; no pairs
# ---------------------------------------------------------------------------

@testset "2. single cluster" begin
    ch = ClusterHistogramSet()
    cs = [mk(ec=2.448, es=2.448, z=50.0)]
    update_cluster_histograms!(ch, cs, params)
    @test ch.N_clusters_counts[2] == 1            # 1 cluster → bin 2
    bin_E = floor(Int, 2.448 / 2.7 * 270) + 1
    @test ch.Ec_counts[bin_E]   == 1
    @test ch.Emax_counts[bin_E] == 1
    @test ch.Emin_counts[bin_E] == 1
    @test ch.Einc_counts[bin_E] == 1
    @test sum(ch.closest_D3_counts) == 0   # no pairs from 1 cluster
end

# ---------------------------------------------------------------------------
# 3. Two clusters: pair distance fills closest_D3 = furthest_D3 (= the pair)
# ---------------------------------------------------------------------------

@testset "3. two clusters: closest_D3 == furthest_D3 (single pair)" begin
    ch = ClusterHistogramSet()
    cs = [mk(x=0.0, y=0.0, z=30.0, ec=1.0, es=1.0),
          mk(x=0.0, y=0.0, z=80.0, ec=0.5, es=0.5)]
    update_cluster_histograms!(ch, cs, params)
    # 3D distance = |Δz| = 50 cm; D3 bin: 50 / 200 * 100 + 1 = 26
    @test ch.closest_D3_counts[26]   == 1
    @test ch.furthest_D3_counts[26]  == 1
    # Δz = 50 cm: |Δz| equals dz_max_cm (50.0), so _bin_idx returns 0
    # (x >= hi → out of range). The dz histogram should remain all zeros.
    @test sum(ch.closest_dz_counts)  == 0
end

@testset "3b. two clusters with smaller Δz, fits in dz bins" begin
    ch = ClusterHistogramSet()
    cs = [mk(x=0.0, y=0.0, z=30.0, ec=1.0, es=1.0),
          mk(x=0.0, y=0.0, z=60.0, ec=0.5, es=0.5)]
    update_cluster_histograms!(ch, cs, params)
    # |Δz| = 30 → bin 61 (1-indexed: floor(30/50*100)+1)
    @test ch.closest_dz_counts[61]   == 1
    @test ch.furthest_dz_counts[61]  == 1
end

# ---------------------------------------------------------------------------
# 4. Emax / Emin per event
# ---------------------------------------------------------------------------

@testset "4. Emax/Emin/Einc per event" begin
    ch = ClusterHistogramSet()
    # Keep the inclusive sum < E_max_MeV (= 2.7) so it lands in-range.
    cs = [mk(ec=0.3, es=0.3),
          mk(ec=0.5, es=0.5),
          mk(ec=1.0, es=1.0)]
    update_cluster_histograms!(ch, cs, params)
    bin = E -> floor(Int, E / 2.7 * 270) + 1
    @test ch.Emax_counts[bin(1.0)] == 1
    @test ch.Emin_counts[bin(0.3)] == 1
    @test ch.Einc_counts[bin(0.3 + 0.5 + 1.0)] == 1
    @test ch.Ec_counts[bin(0.3)] == 1
    @test ch.Ec_counts[bin(0.5)] == 1
    @test ch.Ec_counts[bin(1.0)] == 1
    @test sum(ch.Ec_counts) == 3
end

# ---------------------------------------------------------------------------
# 5. r²-vs-z and D-vs-z 2D heatmaps
# ---------------------------------------------------------------------------

@testset "5. 2D heatmaps r²-vs-z and D-vs-z" begin
    ch = ClusterHistogramSet()
    # r = 30 cm → r² = 900 → bin = floor(900 / 6740 * 100) + 1
    cs = [mk(x=30.0, y=0.0, z=50.0, ec=1.0, es=1.0)]
    update_cluster_histograms!(ch, cs, params)
    iz  = floor(Int, (50.0 - ch.z_min_cm) / (ch.z_max_cm - ch.z_min_cm) * 100) + 1
    ir2 = floor(Int, 900.0 / (82.1^2) * 100) + 1
    iD  = floor(Int, 30.0 / 100.0 * 100) + 1
    @test ch.r2_vs_z_2d_counts[ir2, iz] == 1
    @test ch.D_vs_z_2d_counts[iD, iz]   == 1
end

# ---------------------------------------------------------------------------
# 6. merge_cluster_histograms!
# ---------------------------------------------------------------------------

@testset "6. merge_cluster_histograms!" begin
    ch1 = ClusterHistogramSet()
    ch2 = ClusterHistogramSet()
    update_cluster_histograms!(ch1, [mk(ec=1.0, es=1.0)], params)
    update_cluster_histograms!(ch2, [mk(ec=2.0, es=2.0)], params)
    merge_cluster_histograms!(ch1, ch2)
    @test sum(ch1.Ec_counts) == 2
    bin = E -> floor(Int, E / 2.7 * 270) + 1
    @test ch1.Ec_counts[bin(1.0)] == 1
    @test ch1.Ec_counts[bin(2.0)] == 1
end

println("\n  ── test_cluster_histograms3.jl: ClusterHistogramSet OK ──\n")
