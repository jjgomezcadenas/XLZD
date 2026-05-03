# test/test_field_cage3.jl — Field-cage component masses, activities,
# γ-rates against bb0nu Table I; geometry sanity for FCTG/FCBG annular slabs.

using Test
using Printf

isdefined(Main, :XLZD3) || include("../src3/XLZD3.jl")
using .XLZD3

ti_path     = joinpath(@__DIR__, "..", "data", "nist_ti.csv")
fe_path     = joinpath(@__DIR__, "..", "data", "nist_fe.csv")
teflon_path = joinpath(@__DIR__, "..", "data", "nist_teflon.csv")

mat_Ti   = load_material("Ti",   4.510, ti_path)
mat_SS   = load_material("SS",   7.930, fe_path)
mat_PTFE = load_material("PTFE", 2.200, teflon_path)
fc       = build_field_cage(mat_Ti, mat_SS, mat_PTFE)

@testset "Component identity and order" begin
    comps = fc_components(fc)
    @test length(comps) == 6
    names = [p.name for p in comps]
    @test names == ["FCRN", "FCRS", "FCSE", "FCPT", "FCTG", "FCBG"]

    @test fc.FCRN isa PCyl
    @test fc.FCRS isa PCyl
    @test fc.FCSE isa PCyl
    @test fc.FCPT isa PCyl
    @test fc.FCTG isa PAnnularDisk
    @test fc.FCBG isa PAnnularDisk
end

@testset "Masses match bb0nu Table I" begin
    # Each component's back-derived geometry should yield exactly the
    # target mass to within float roundoff.
    @test isapprox(mass(fc.FCRN), 93.0;        rtol=1e-10)
    @test isapprox(mass(fc.FCRS),  0.06;       rtol=1e-10)
    @test isapprox(mass(fc.FCSE),  5.02;       rtol=1e-10)
    @test isapprox(mass(fc.FCPT), 184.0;       rtol=1e-10)
    @test isapprox(mass(fc.FCTG),  44.55;      rtol=1e-10)
    @test isapprox(mass(fc.FCBG),  22.275;     rtol=1e-10)

    # Total: 93 + 0.06 + 5.02 + 184 + 44.55 + 22.275 = 348.905 kg
    @test isapprox(fc_total_mass(fc), 348.905; rtol=1e-10)
end

@testset "Activities match bb0nu Table I" begin
    # Bi-214 chain
    @test isapprox(activity_U238_late(fc.FCRN), 93.0   * 0.35;    rtol=1e-10)
    @test isapprox(activity_U238_late(fc.FCRS),  0.06  * 1350.0;  rtol=1e-10)
    @test isapprox(activity_U238_late(fc.FCSE),  5.02  * 5.82;    rtol=1e-10)
    @test isapprox(activity_U238_late(fc.FCPT), 184.0  * 0.04;    rtol=1e-10)
    @test isapprox(activity_U238_late(fc.FCTG),  44.55 * 2.63;    rtol=1e-10)
    @test isapprox(activity_U238_late(fc.FCBG),  22.275* 2.63;    rtol=1e-10)

    # Tl-208 chain (Th-late)
    @test isapprox(activity_Th232_late(fc.FCRN), 93.0   * 0.24;   rtol=1e-10)
    @test isapprox(activity_Th232_late(fc.FCRS),  0.06  * 2010.0; rtol=1e-10)
    @test isapprox(activity_Th232_late(fc.FCSE),  5.02  * 1.88;   rtol=1e-10)
    @test isapprox(activity_Th232_late(fc.FCPT), 184.0  * 0.01;   rtol=1e-10)
    @test isapprox(activity_Th232_late(fc.FCTG),  44.55 * 1.46;   rtol=1e-10)
    @test isapprox(activity_Th232_late(fc.FCBG),  22.275* 1.46;   rtol=1e-10)
end

@testset "Barrel geometry: radii and z extent" begin
    # All barrel components share the FC barrel z extent
    for p in (fc.FCRN, fc.FCRS, fc.FCSE, fc.FCPT)
        @test p.geom.z_min ≈ FC_BARREL_Z_BOT_CM
        @test p.geom.z_max ≈ FC_BARREL_Z_TOP_CM
        @test height(p.geom) ≈ FC_BARREL_HEIGHT_CM
    end
    # Rings/resistors/sensors at R=74.3 (FC outer); PTFE at R=72.8 (TPC inner)
    @test fc.FCRN.geom.R_inner ≈ FC_R_RING_CM
    @test fc.FCRS.geom.R_inner ≈ FC_R_RING_CM
    @test fc.FCSE.geom.R_inner ≈ FC_R_RING_CM
    @test fc.FCPT.geom.R_inner ≈ FC_R_TPC_INNER_CM
end

@testset "Back-derived wall thickness sanity" begin
    # FCRN (Ti, ρ=4.51): t ≈ 0.277 cm
    @test isapprox(fc.FCRN.geom.wall_thickness, 0.277; rtol=5e-2)
    # FCRS (SS, 60 g): essentially 0
    @test fc.FCRS.geom.wall_thickness < 1e-3
    # FCSE (SS, 5 kg): t ≈ 0.0087 cm
    @test isapprox(fc.FCSE.geom.wall_thickness, 0.0087; rtol=5e-2)
    # FCPT (PTFE, 184 kg): t ≈ 1.15 cm
    @test isapprox(fc.FCPT.geom.wall_thickness, 1.15; rtol=5e-2)
end

@testset "FCTG/FCBG annular slab geometry" begin
    # Both share inner=72.8 (TPC inner), outer=80.3 (skin outer)
    @test fc.FCTG.R_in  ≈ FC_R_TPC_INNER_CM
    @test fc.FCTG.R_out ≈ FC_R_HOLDER_OUT_CM
    @test fc.FCBG.R_in  ≈ FC_R_TPC_INNER_CM
    @test fc.FCBG.R_out ≈ FC_R_HOLDER_OUT_CM

    # FCTG: face at gate plane, slab opens UP (into gas/anode region) —
    # inward normal = −ẑ (emit into LXe drift below)
    @test fc.FCTG.z_face ≈ FC_Z_GATE_CM
    @test fc.FCTG.normal_sign == -1
    # FCBG: face at cathode plane, slab opens DOWN (into RFR) —
    # inward normal = +ẑ (emit into LXe drift above)
    @test fc.FCBG.z_face ≈ FC_Z_CATHODE_CM
    @test fc.FCBG.normal_sign == +1

    # Slab heights from mass / (ρ_SS × annular footprint)
    expected_H_FCTG = 44.55 * 1000 / (mat_SS.density *
                       π * (FC_R_HOLDER_OUT_CM^2 - FC_R_TPC_INNER_CM^2))
    expected_H_FCBG = 22.275 * 1000 / (mat_SS.density *
                       π * (FC_R_HOLDER_OUT_CM^2 - FC_R_TPC_INNER_CM^2))
    @test isapprox(fc.FCTG.H_axial, expected_H_FCTG; rtol=1e-10)
    @test isapprox(fc.FCBG.H_axial, expected_H_FCBG; rtol=1e-10)
    # Sanity: design values 1.56 cm and 0.78 cm
    @test isapprox(fc.FCTG.H_axial, 1.558; rtol=2e-3)
    @test isapprox(fc.FCBG.H_axial, 0.779; rtol=2e-3)
end

@testset "Slab regions stay clear of active LXe" begin
    # FCTG slab spans z ∈ [z_gate, z_gate + H], i.e. fully ABOVE the gate
    z_top_FCTG = fc.FCTG.z_face + fc.FCTG.H_axial
    @test z_top_FCTG > FC_Z_GATE_CM           # opens upward
    # FCBG slab spans z ∈ [-H, 0], i.e. fully BELOW the cathode
    z_bot_FCBG = fc.FCBG.z_face - fc.FCBG.H_axial
    @test z_bot_FCBG < FC_Z_CATHODE_CM        # opens downward
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
    μt_FCRN = μ_Ti   * fc.FCRN.geom.wall_thickness
    μt_FCRS = μ_SS   * fc.FCRS.geom.wall_thickness
    μt_FCSE = μ_SS   * fc.FCSE.geom.wall_thickness
    μt_FCPT = μ_PTFE * fc.FCPT.geom.wall_thickness
    μt_FCTG = μ_SS   * fc.FCTG.H_axial
    μt_FCBG = μ_SS   * fc.FCBG.H_axial

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
