# test/test_pobjects3.jl — Verify PCyl/PDisk masses, activities, γ-rates.

using Test
using Printf

isdefined(Main, :XLZD3) || include("../src3/XLZD3.jl")
using .XLZD3

ti_path     = joinpath(@__DIR__, "..", "data", "nist_ti.csv")
geom_csv    = joinpath(@__DIR__, "..", "data", "lz_cryo_geometry.csv")
extras_csv  = joinpath(@__DIR__, "..", "data", "lz_cryo_extras.csv")

mat_Ti = load_material("Ti", 4.510, ti_path)
cryo   = build_cryostat(geom_csv, extras_csv)

@testset "PCyl/PDisk basic" begin
    g = GCyl(10.0, 0.5, 0.0, 100.0)
    p = PCyl(g, mat_Ti, 0.08, 0.22; count=2, name="test")
    # mass = volume × density × count
    @test mass(p) ≈ mass(g, mat_Ti.density) * 2
    # activities
    @test activity_U238_late(p) ≈ mass(p) * 0.08
    @test activity_Th232_late(p) ≈ mass(p) * 0.22
    # γ rates
    @test gamma_rate_Bi214(p) ≈
        activity_U238_late(p) * 1e-3 * BR_BI214_GAMMA * SEC_PER_YEAR
    @test gamma_rate_Tl208(p) ≈
        activity_Th232_late(p) * 1e-3 * BR_TL208_FROM_CHAIN * SEC_PER_YEAR

    d = GDisk(50.0, 0.8, 0.0, 2.0, :up)
    pd = PDisk(d, mat_Ti, 0.08, 0.22; count=1, name="test_disc")
    @test mass(pd) ≈ mass(d, mat_Ti.density)
    @test activity_U238_late(pd) ≈ mass(pd) * 0.08
end

@testset "Total Ti activity matches bb0nu Table I" begin
    # All Ti: barrels + heads + ALL extras (mc_only=false) → 2572 kg
    ps_all = pobjects_from_cryostat(cryo, mat_Ti; mc_only=false)
    m_total = sum(mass(p) for p in ps_all)
    @test isapprox(m_total, 2590.0; rtol=0.05)

    # bb0nu Ti specific activities → total activity per chain
    aU_total  = sum(activity_U238_late(p)  for p in ps_all)
    aTh_total = sum(activity_Th232_late(p) for p in ps_all)
    @test isapprox(aU_total,  m_total * TI_BB0NU_U238_LATE_MBQKG;  rtol=1e-12)
    @test isapprox(aTh_total, m_total * TI_BB0NU_TH232_LATE_MBQKG; rtol=1e-12)

    # γ production rates
    γ_Bi214 = sum(gamma_rate_Bi214(p) for p in ps_all)
    γ_Tl208 = sum(gamma_rate_Tl208(p) for p in ps_all)
    # Expected: 2590 kg × 0.08 mBq/kg × 1e-3 × 0.0155 × 3.1557e7 ≈ 1.013e5 γ/yr
    @test isapprox(γ_Bi214, 1.013e5; rtol=0.05)
    # 2590 × 0.22 × 1e-3 × 0.359 × 3.1557e7 ≈ 6.456e6 γ/yr
    @test isapprox(γ_Tl208, 6.456e6; rtol=0.05)

    # Documentation block
    println()
    println("── Total Ti activity / γ rate (full cryostat) ──")
    @printf("  Ti mass                     %8.1f kg\n", m_total)
    @printf("  ²³⁸U-late activity          %8.2f mBq\n", aU_total)
    @printf("  ²³²Th-late activity         %8.2f mBq\n", aTh_total)
    @printf("  Bi-214 γ rate (2.448 MeV)   %.3e γ/yr\n", γ_Bi214)
    @printf("  Tl-208 γ rate (2.615 MeV)   %.3e γ/yr\n", γ_Tl208)
    println()
end

@testset "MC-active subset" begin
    ps_mc = pobjects_from_cryostat(cryo, mat_Ti; mc_only=true)
    # 2 barrels + 4 heads + 11 active extras = 17
    @test length(ps_mc) == 17
    # MC-active mass ~84% of total Ti
    m_mc = sum(mass(p) for p in ps_mc)
    @test 0.80 < m_mc / 2572.1 < 0.86
end

@testset "PSurface basic + MLI activity" begin
    s = PSurface("MLI_test", :ICV_body, "MLI", 13.8, 11.1, 7.79)
    @test mass(s) == 13.8
    @test activity_U238_late(s) ≈ 13.8 * 11.1
    @test activity_Th232_late(s) ≈ 13.8 * 7.79
    # Bi-214 γ produced = 13.8 × 11.1 × 1e-3 × 0.0155 × 3.1557e7 ≈ 7.495e4 γ/yr
    @test isapprox(gamma_rate_Bi214(s), 7.495e4; rtol=1e-3)
    # Tl-208: 13.8 × 7.79 × 1e-3 × 0.359 × 3.1557e7 ≈ 1.218e6 γ/yr
    @test isapprox(gamma_rate_Tl208(s), 1.218e6; rtol=1e-3)

    # Surface spectrum is constant in u (= R/2 for u > 0)
    u_bins = collect(range(0.005, 0.995, length=10))
    R = gamma_rate_Bi214(s)
    dNdu = self_shielded_spectrum(s, R, E_BI214_MEV, u_bins)
    @test all(d ≈ R / 2 for d in dNdu)
end

@testset "Count multiplicity" begin
    # Tie-bar ports (count=3) and conduit ports (count=2) are in the
    # MC-active set. Verify mass scales by count.
    ps_mc = pobjects_from_cryostat(cryo, mat_Ti; mc_only=true)
    tb = only(p for p in ps_mc if p.name == "OCV_top_tie_bar_ports")
    cd = only(p for p in ps_mc if p.name == "OCV_top_conduit_ports")
    @test tb.count == 3
    @test cd.count == 2
    @test mass(tb) ≈ mass(tb.geom, tb.material.density) * 3
    @test mass(cd) ≈ mass(cd.geom, cd.material.density) * 2
end
