# test/test_field_cage3.jl — Field-cage component masses, activities,
# γ-rates against bb0nu Table I; geometry sanity for FCTG/FCBG annular slabs.

using Test
using Printf

isdefined(Main, :XLZD3) || include("../src3/XLZD3.jl")
using .XLZD3

ti_path        = joinpath(@__DIR__, "..", "data", "nist_ti.csv")
fe_path        = joinpath(@__DIR__, "..", "data", "nist_fe.csv")
teflon_path    = joinpath(@__DIR__, "..", "data", "nist_teflon.csv")
fc_barrels_csv = joinpath(@__DIR__, "..", "data", "lz_fc_barrels.csv")
fc_grids_csv   = joinpath(@__DIR__, "..", "data", "lz_fc_grids.csv")

mat_Ti   = load_material("Ti",   4.510, ti_path)
mat_SS   = load_material("SS",   7.930, fe_path)
mat_PTFE = load_material("PTFE", 2.200, teflon_path)
fc       = build_field_cage(fc_barrels_csv, fc_grids_csv,
                             Dict("Ti"   => mat_Ti,
                                  "SS"   => mat_SS,
                                  "PTFE" => mat_PTFE))

named = Dict(p.name => p for p in fc_components(fc))

@testset "Component identity and order" begin
    comps = fc_components(fc)
    @test length(comps) == 6
    names = [p.name for p in comps]
    @test names == ["FCRN", "FCRS", "FCSE", "FCPT", "FCTG", "FCBG"]

    # Barrel components are PCyl, grids are PAnnularDisk
    for n in ("FCRN", "FCRS", "FCSE", "FCPT")
        @test named[n] isa PCyl
    end
    for n in ("FCTG", "FCBG")
        @test named[n] isa PAnnularDisk
    end

    # FieldCage struct exposes barrels + grids vectors
    @test length(fc.barrels) == 4
    @test length(fc.grids)   == 2
end

@testset "Masses match bb0nu Table I (CSV-driven)" begin
    # Each component's back-derived geometry should yield exactly the
    # mass listed in the source CSV row.
    @test isapprox(mass(named["FCRN"]),  93.0;   rtol=1e-10)
    @test isapprox(mass(named["FCRS"]),   0.06;  rtol=1e-10)
    @test isapprox(mass(named["FCSE"]),   5.02;  rtol=1e-10)
    @test isapprox(mass(named["FCPT"]), 184.0;   rtol=1e-10)
    @test isapprox(mass(named["FCTG"]),  44.55;  rtol=1e-10)
    @test isapprox(mass(named["FCBG"]),  22.275; rtol=1e-10)

    # Total: 93 + 0.06 + 5.02 + 184 + 44.55 + 22.275 = 348.905 kg
    @test isapprox(fc_total_mass(fc), 348.905; rtol=1e-10)
end

@testset "Activities match bb0nu Table I" begin
    # Bi-214 chain
    @test isapprox(activity_U238_late(named["FCRN"]),  93.0  * 0.35;   rtol=1e-10)
    @test isapprox(activity_U238_late(named["FCRS"]),   0.06 * 1350.0; rtol=1e-10)
    @test isapprox(activity_U238_late(named["FCSE"]),   5.02 * 5.82;   rtol=1e-10)
    @test isapprox(activity_U238_late(named["FCPT"]), 184.0  * 0.04;   rtol=1e-10)
    @test isapprox(activity_U238_late(named["FCTG"]),  44.55 * 2.63;   rtol=1e-10)
    @test isapprox(activity_U238_late(named["FCBG"]),  22.275* 2.63;   rtol=1e-10)

    # Tl-208 chain (Th-late)
    @test isapprox(activity_Th232_late(named["FCRN"]),  93.0  * 0.24;  rtol=1e-10)
    @test isapprox(activity_Th232_late(named["FCRS"]),   0.06 * 2010.0; rtol=1e-10)
    @test isapprox(activity_Th232_late(named["FCSE"]),   5.02 * 1.88;  rtol=1e-10)
    @test isapprox(activity_Th232_late(named["FCPT"]), 184.0  * 0.01;  rtol=1e-10)
    @test isapprox(activity_Th232_late(named["FCTG"]),  44.55 * 1.46;  rtol=1e-10)
    @test isapprox(activity_Th232_late(named["FCBG"]),  22.275* 1.46;  rtol=1e-10)
end

@testset "Barrel geometry: radii and z extent (CSV-driven)" begin
    # All barrel components share the FC barrel z extent (per CSV)
    for n in ("FCRN", "FCRS", "FCSE", "FCPT")
        p = named[n]
        @test p.geom.z_min ≈ -13.75
        @test p.geom.z_max ≈ 145.6
        @test height(p.geom) ≈ 159.35
    end
    # Rings/resistors/sensors at R=74.3 (FC outer); PTFE at R=72.8 (TPC inner)
    @test named["FCRN"].geom.R_inner ≈ 74.3
    @test named["FCRS"].geom.R_inner ≈ 74.3
    @test named["FCSE"].geom.R_inner ≈ 74.3
    @test named["FCPT"].geom.R_inner ≈ 72.8
end

@testset "Back-derived wall thickness sanity" begin
    # FCRN (Ti, ρ=4.51): t ≈ 0.277 cm
    @test isapprox(named["FCRN"].geom.wall_thickness, 0.277; rtol=5e-2)
    # FCRS (SS, 60 g): essentially 0
    @test named["FCRS"].geom.wall_thickness < 1e-3
    # FCSE (SS, 5 kg): t ≈ 0.0087 cm
    @test isapprox(named["FCSE"].geom.wall_thickness, 0.0087; rtol=5e-2)
    # FCPT (PTFE, 184 kg): t ≈ 1.15 cm
    @test isapprox(named["FCPT"].geom.wall_thickness, 1.15; rtol=5e-2)
end

@testset "FCTG/FCBG annular slab geometry (CSV-driven)" begin
    # Both share inner=72.8 (TPC inner), outer=80.3 (skin outer)
    @test named["FCTG"].R_in  ≈ 72.8
    @test named["FCTG"].R_out ≈ 80.3
    @test named["FCBG"].R_in  ≈ 72.8
    @test named["FCBG"].R_out ≈ 80.3

    # FCTG: face at gate plane (z=145.6), slab opens UP (normal -z)
    @test named["FCTG"].z_face ≈ 145.6
    @test named["FCTG"].normal_sign == -1
    # FCBG: face at cathode plane (z=0), slab opens DOWN (normal +z)
    @test named["FCBG"].z_face ≈ 0.0
    @test named["FCBG"].normal_sign == +1

    # Slab heights from mass / (ρ_SS × annular footprint)
    expected_H_FCTG = 44.55  * 1000 / (mat_SS.density * π * (80.3^2 - 72.8^2))
    expected_H_FCBG = 22.275 * 1000 / (mat_SS.density * π * (80.3^2 - 72.8^2))
    @test isapprox(named["FCTG"].H_axial, expected_H_FCTG; rtol=1e-10)
    @test isapprox(named["FCBG"].H_axial, expected_H_FCBG; rtol=1e-10)
    # Sanity: design values 1.56 cm and 0.78 cm
    @test isapprox(named["FCTG"].H_axial, 1.558; rtol=2e-3)
    @test isapprox(named["FCBG"].H_axial, 0.779; rtol=2e-3)
end

@testset "Slab regions stay clear of active LXe" begin
    # FCTG slab spans z ∈ [z_gate, z_gate + H], i.e. fully ABOVE the gate
    z_top_FCTG = named["FCTG"].z_face + named["FCTG"].H_axial
    @test z_top_FCTG > 145.6                      # opens upward
    # FCBG slab spans z ∈ [-H, 0], i.e. fully BELOW the cathode
    z_bot_FCBG = named["FCBG"].z_face - named["FCBG"].H_axial
    @test z_bot_FCBG < 0.0                        # opens downward
end

@testset "γ-rates per source" begin
    println()
    println("── Field-cage γ rates (γ/yr) at full 4π ──")
    @printf("  %-6s %10s %10s %10s\n", "src", "mass(kg)", "Bi-214", "Tl-208")
    println("  ", "─"^46)
    for p in fc_components(fc)
        # 4π rates (no inward halving, no shielding) — sanity inputs
        γ_Bi = gamma_rate_Bi214(p)
        γ_Tl = gamma_rate_Tl208(p)
        @printf("  %-6s %10.4f %10.3e %10.3e\n", p.name, mass(p), γ_Bi, γ_Tl)
        # Both rates positive
        @test γ_Bi > 0
        @test γ_Tl > 0
    end
    println()

    # Total Bi-214 = (93×0.35 + 0.06×1350 + 5.02×5.82 + 184×0.04 + 44.55×2.63 + 22.275×2.63)
    #              = 32.55 + 81.0 + 29.22 + 7.36 + 117.16 + 58.58 = 325.87 mBq
    expected_total_Bi_mBq = 93.0*0.35 + 0.06*1350 + 5.02*5.82 + 184*0.04 +
                            44.55*2.63 + 22.275*2.63
    actual_total_Bi_mBq = sum(activity_U238_late(p) for p in fc_components(fc))
    @test isapprox(actual_total_Bi_mBq, expected_total_Bi_mBq; rtol=1e-10)
end

@testset "Self-shielding non-negligible only for FCRN/FCPT/FCTG/FCBG" begin
    # μ·t at 2.448 MeV; FCRS and FCSE should have negligible shielding
    # (μt < 0.05 per the py/lz_fieldcage.py analysis).
    μ_Ti   = mat_Ti.μ_lin(E_BI214_MEV)
    μ_SS   = mat_SS.μ_lin(E_BI214_MEV)
    μ_PTFE = mat_PTFE.μ_lin(E_BI214_MEV)
    μt_FCRN = μ_Ti   * named["FCRN"].geom.wall_thickness
    μt_FCRS = μ_SS   * named["FCRS"].geom.wall_thickness
    μt_FCSE = μ_SS   * named["FCSE"].geom.wall_thickness
    μt_FCPT = μ_PTFE * named["FCPT"].geom.wall_thickness
    μt_FCTG = μ_SS   * named["FCTG"].H_axial
    μt_FCBG = μ_SS   * named["FCBG"].H_axial

    @test μt_FCRS < 0.05
    @test μt_FCSE < 0.05
    @test μt_FCRN > 0.04          # ≈ 0.05
    @test μt_FCPT > 0.05          # ≈ 0.10
    @test μt_FCTG > 0.10          # ≈ 0.49
    @test μt_FCBG > 0.05          # ≈ 0.24

    println()
    println("── μ·t at 2.448 MeV (self-shielding optical depth) ──")
    @printf("  FCRN  %.4f   FCRS  %.5f   FCSE  %.4f\n", μt_FCRN, μt_FCRS, μt_FCSE)
    @printf("  FCPT  %.4f   FCTG  %.4f   FCBG  %.4f\n", μt_FCPT, μt_FCTG, μt_FCBG)
    println()
end
