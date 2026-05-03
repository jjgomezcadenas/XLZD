# test/test_fc_effective_sources3.jl — FC EffectiveSources + sampling
# primitives. Verifies the 6 FC components × 3 isotopes are wired
# correctly and that the new :fc_barrel / :fc_annular sampling branches
# emit at the right geometry with inward-going directions.

using Test
using Printf
using Random

isdefined(Main, :XLZD3) || include("../src3/XLZD3.jl")
using .XLZD3

ti_path     = joinpath(@__DIR__, "..", "data", "nist_ti.csv")
fe_path     = joinpath(@__DIR__, "..", "data", "nist_fe.csv")
teflon_path = joinpath(@__DIR__, "..", "data", "nist_teflon.csv")
lxe_path    = joinpath(@__DIR__, "..", "data", "nist_lxe.csv")
lxe_csv     = joinpath(@__DIR__, "..", "data", "lxe_detector.csv")

mat_Ti   = load_material("Ti",   4.510, ti_path)
mat_SS   = load_material("SS",   7.930, fe_path)
mat_PTFE = load_material("PTFE", 2.200, teflon_path)
mat_LXe  = load_material("LXe",  2.953, lxe_path)
det      = build_lxe_detector(lxe_csv, mat_LXe)

fc       = build_field_cage(mat_Ti, mat_SS, mat_PTFE)
fc_effs  = build_field_cage_effective_sources(fc)

@testset "FC effective sources: count and names" begin
    @test length(fc_effs) == 6 * 3      # 6 components × 3 isotopes
    expected = String[]
    for cmp in ("FCRN", "FCRS", "FCSE", "FCPT", "FCTG", "FCBG")
        for iso in ("Bi214", "Tl208", "Tl208c")
            push!(expected, "$(cmp)_$(iso)")
        end
    end
    actual = [e.name for e in fc_effs]
    @test sort(actual) == sort(expected)
end

@testset "FC effective sources: regions and isotopes" begin
    by_name = Dict(e.name => e for e in fc_effs)
    # Barrel components → :fc_barrel
    for cmp in ("FCRN", "FCRS", "FCSE", "FCPT")
        for iso in ("Bi214", "Tl208", "Tl208c")
            e = by_name["$(cmp)_$(iso)"]
            @test e.region == :fc_barrel
        end
    end
    # Grid holders → :fc_annular
    for cmp in ("FCTG", "FCBG")
        for iso in ("Bi214", "Tl208", "Tl208c")
            e = by_name["$(cmp)_$(iso)"]
            @test e.region == :fc_annular
        end
    end
    # Isotope tag matches the name suffix
    for e in fc_effs
        if endswith(e.name, "_Bi214")
            @test e.isotope === :Bi214
            @test e.E_MeV == E_BI214_MEV
        elseif endswith(e.name, "_Tl208c")
            @test e.isotope === :Tl208c
            @test e.E_MeV == E_TL208_COMPANION_MEV
        elseif endswith(e.name, "_Tl208")
            @test e.isotope === :Tl208
            @test e.E_MeV == E_TL208_MEV
        end
    end
end

@testset "FC effective sources: single contribution, empty downstream" begin
    for e in fc_effs
        @test length(e.contributions) == 1
        @test isempty(e.contributions[1].downstream_slabs)
    end
end

@testset "FC dN/du integrates to total_per_yr" begin
    for e in fc_effs
        # Trapezoidal integration of dN/du over u_bins
        s = 0.0
        for i in 1:length(e.u_bins)-1
            s += 0.5 * (e.dNdu[i] + e.dNdu[i+1]) *
                       (e.u_bins[i+1] - e.u_bins[i])
        end
        @test isapprox(s, e.total_per_yr; rtol=1e-12)
    end
end

@testset "FC barrel sampling: position on shell, direction inward" begin
    rng = MersenneTwister(20251104)
    by_name = Dict(e.name => e for e in fc_effs)
    e = by_name["FCRN_Bi214"]
    cdf = build_cdf(e.u_bins, e.dNdu)

    # PCyl geometry as recorded
    p = e.contributions[1].source.producer::PCyl
    R = p.geom.R_inner
    zmin, zmax = p.geom.z_min, p.geom.z_max

    for _ in 1:300
        x, y, z, dx, dy, dz, u = sample_entry(rng, det, e, cdf)
        # Position on the inner cylindrical surface
        r = sqrt(x*x + y*y)
        @test isapprox(r, R; rtol=1e-10)
        @test zmin <= z <= zmax
        # Direction unit length
        @test isapprox(dx*dx + dy*dy + dz*dz, 1.0; rtol=1e-10)
        # Inward-going (dot with inward normal −r̂ should be positive)
        @test (-dx*x - dy*y) / r > 0
        # Sampled u in [0, 1]
        @test 0 ≤ u ≤ 1
    end
end

@testset "FC annular sampling: position on annulus, direction inward" begin
    rng = MersenneTwister(98765)
    by_name = Dict(e.name => e for e in fc_effs)

    for src_name in ("FCTG_Bi214", "FCBG_Bi214")
        e = by_name[src_name]
        cdf = build_cdf(e.u_bins, e.dNdu)
        p = e.contributions[1].source.producer::PAnnularDisk

        for _ in 1:300
            x, y, z, dx, dy, dz, u = sample_entry(rng, det, e, cdf)
            r2 = x*x + y*y
            # Inside the annulus (allowing for floating roundoff at edges)
            @test p.R_in^2 - 1e-9 ≤ r2 ≤ p.R_out^2 + 1e-9
            # On the LXe-facing face (z is exact)
            @test z == p.z_face
            # Unit-length direction
            @test isapprox(dx*dx + dy*dy + dz*dz, 1.0; rtol=1e-10)
            # Inward-going: dz · normal_sign > 0
            @test dz * p.normal_sign > 0
            @test 0 ≤ u ≤ 1
        end
    end
end

@testset "FC annular: r distribution uniform in r²" begin
    # Verify that annular sampling is uniform in r² (the way uniform-by-area
    # samples a flat annulus). A coarse two-bin chi-square is enough.
    rng = MersenneTwister(123)
    by_name = Dict(e.name => e for e in fc_effs)
    e = by_name["FCBG_Bi214"]
    cdf = build_cdf(e.u_bins, e.dNdu)
    p = e.contributions[1].source.producer::PAnnularDisk

    N = 5000
    r2_mid = 0.5 * (p.R_in^2 + p.R_out^2)
    n_inner = 0
    for _ in 1:N
        x, y, _z, _dx, _dy, _dz, _u = sample_entry(rng, det, e, cdf)
        if x*x + y*y < r2_mid
            n_inner += 1
        end
    end
    # Inner half (in r²) should hold ~½ of the points; allow ±5σ
    σ = sqrt(N * 0.25)
    @test abs(n_inner - N/2) < 5σ
end

@testset "Sanity: γ/yr per FC source vs. py/lz_fieldcage.py" begin
    # Inward γ/yr (after self-shielding) at each source's exit surface.
    # Compare against the half-rate (mass × spec_act × BR × ½) for the
    # thin-source limit; FCTG/FCBG should fall short due to self-shielding.
    by_name = Dict(e.name => e for e in fc_effs)
    println()
    println("── FC Bi-214 inward γ/yr (post self-shielding) ──")
    @printf("  %-6s  %12s  %12s  %5s\n",
            "src", "inward γ/yr", "no-shield ½R", "ratio")
    println("  ", "─"^50)
    for cmp in ("FCRN", "FCRS", "FCSE", "FCPT", "FCTG", "FCBG")
        e = by_name["$(cmp)_Bi214"]
        # Total γ/yr produced (4π) for this source
        p = e.contributions[1].source.producer
        R_total = gamma_rate_Bi214(p)
        no_shield_inward = 0.5 * R_total
        ratio = e.total_per_yr / no_shield_inward
        @printf("  %-6s  %12.3e  %12.3e  %5.2f\n",
                cmp, e.total_per_yr, no_shield_inward, ratio)
        # Sanity: inward fraction should be in (0, 0.5] (strictly < 0.5
        # for shielded sources; ≈ 0.5 for negligible shielding).
        @test 0 < ratio ≤ 1.0
    end
    println()
end
