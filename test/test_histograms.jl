# test/test_histograms.jl — Verify the control-histogram pipeline.

using Test
using Random
using Printf

include("../src2/XLZD2.jl")
using .XLZD2

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
    push!(s.deposits, LXeDeposit(0.0, 0.0, 5.0, 0.5, :active))
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
    push!(s.deposits, LXeDeposit(0.0, 0.0, 5.00, 0.5, :active))
    push!(s.deposits, LXeDeposit(0.0, 0.0, 5.10, 0.3, :active))   # 1 mm away
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
    push!(s.deposits, LXeDeposit(0.0, 0.0,  5.0, 0.5, :active))
    push!(s.deposits, LXeDeposit(0.0, 0.0, 10.0, 0.7, :active))   # 5 cm away
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

    push!(s.deposits, LXeDeposit(0.0, 0.0,  5.0, 0.5, :active));  update_histograms!(h1, s, p)
    empty!(s.deposits)
    push!(s.deposits, LXeDeposit(0.0, 0.0, 10.0, 0.7, :active)); update_histograms!(h2, s, p)

    merge_histograms!(h1, h2)
    @test h1.ssms_counts == [2, 0, 0]
    @test sum(h1.N_clusters_counts) == 2
end

# ---------------------------------------------------------------------------
# Integration: run a small MC and verify histogram bookkeeping
# ---------------------------------------------------------------------------

@testset "run_mc with_histograms=true populates the histograms" begin
    lxe_csv      = joinpath(@__DIR__, "..", "data", "lxe_detector.csv")
    lxe_nist     = joinpath(@__DIR__, "..", "data", "nist_lxe.csv")
    ti_path      = joinpath(@__DIR__, "..", "data", "nist_ti.csv")
    geom_csv     = joinpath(@__DIR__, "..", "data", "lz_cryo_geometry.csv")
    extras_csv   = joinpath(@__DIR__, "..", "data", "lz_cryo_extras.csv")
    surfaces_csv = joinpath(@__DIR__, "..", "data", "lz_cryo_surface_sources.csv")
    xcom_path    = joinpath(@__DIR__, "..", "data", "nist.csv")

    mat_LXe = load_material("LXe", 2.953, lxe_nist)
    mat_Ti  = load_material("Ti",  4.510, ti_path)
    det     = build_lxe_detector(lxe_csv, mat_LXe)
    cryo    = build_cryostat(geom_csv, extras_csv, surfaces_csv)
    indiv   = build_individual_sources(cryo, mat_Ti)
    effs    = build_effective_sources(indiv, cryo, mat_Ti)
    xcom    = load_xcom(xcom_path)
    params  = MCParams()
    by_name = Dict(e.name => e for e in effs)

    eff = by_name["CB_Bi214"]
    res = run_mc(det, eff, nothing, xcom, params, 5000;
                  mc_seed=0xABCD, with_histograms=true)
    @test res.histograms !== nothing
    h = res.histograms
    # SS + MS + no_cluster sums to n_total
    @test sum(h.ssms_counts) == res.n_total
    # N_clusters histogram sums to n_total
    @test sum(h.N_clusters_counts) == res.n_total
    @test sum(h.N_extra_counts)    == res.n_total
    # Some events have visible deposits → E_first non-empty
    @test sum(h.E_first_counts) > 0
    # E_cluster has at least as many entries as SS+MS events
    @test sum(h.E_cluster_counts) >= h.ssms_counts[1] + h.ssms_counts[2]

    println()
    @printf("  CB_Bi214, 5000 photons:\n")
    @printf("    SS=%d  MS=%d  no_cluster=%d\n",
            h.ssms_counts[1], h.ssms_counts[2], h.ssms_counts[3])
    @printf("    Σ E_first = %d   Σ E_cluster = %d\n",
            sum(h.E_first_counts), sum(h.E_cluster_counts))
    println()
end

@testset "run_mc with_histograms=false keeps histograms = nothing" begin
    lxe_csv      = joinpath(@__DIR__, "..", "data", "lxe_detector.csv")
    lxe_nist     = joinpath(@__DIR__, "..", "data", "nist_lxe.csv")
    ti_path      = joinpath(@__DIR__, "..", "data", "nist_ti.csv")
    geom_csv     = joinpath(@__DIR__, "..", "data", "lz_cryo_geometry.csv")
    extras_csv   = joinpath(@__DIR__, "..", "data", "lz_cryo_extras.csv")
    surfaces_csv = joinpath(@__DIR__, "..", "data", "lz_cryo_surface_sources.csv")
    xcom_path    = joinpath(@__DIR__, "..", "data", "nist.csv")

    mat_LXe = load_material("LXe", 2.953, lxe_nist)
    mat_Ti  = load_material("Ti",  4.510, ti_path)
    det     = build_lxe_detector(lxe_csv, mat_LXe)
    cryo    = build_cryostat(geom_csv, extras_csv, surfaces_csv)
    indiv   = build_individual_sources(cryo, mat_Ti)
    effs    = build_effective_sources(indiv, cryo, mat_Ti)
    xcom    = load_xcom(xcom_path)
    params  = MCParams()
    by_name = Dict(e.name => e for e in effs)

    res = run_mc(det, by_name["CB_Bi214"], nothing, xcom, params, 1000;
                  mc_seed=0xBEEF, with_histograms=false)
    @test res.histograms === nothing
end
