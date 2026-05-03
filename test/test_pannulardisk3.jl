# test/test_pannulardisk3.jl — PAnnularDisk geometry, mass, slab integral.
#
# PAnnularDisk is the field-cage grid-holder primitive (FCTG, FCBG):
# a thick-walled annular slab, single LXe-facing emission face, slab
# self-shielding handled by the same integral as PDisk.

using Test
using Printf

isdefined(Main, :XLZD3) || include("../src3/XLZD3.jl")
using .XLZD3

# Material setup. Use the Fe XCOM table as proxy for stainless steel
# (Fe is 70% of SS-304 by mass; Cr/Ni Z values are nearly identical so
# μ/ρ matches to <2%). Density set to ρ_SS = 7.93 g/cm³.
fe_path     = joinpath(@__DIR__, "..", "data", "nist_fe.csv")
teflon_path = joinpath(@__DIR__, "..", "data", "nist_teflon.csv")
mat_SS      = load_material("SS",   7.93, fe_path)
mat_PTFE    = load_material("PTFE", 2.20, teflon_path)

# FCBG cathode-only holder dimensions (per design discussion):
const R_IN_CM     = 72.8
const R_OUT_CM    = 80.3
const FCBG_MASS_KG = 22.28        # ¼ of bb0nu's 89.1 kg "field grids and holders"
const FCBG_BI_MBQKG = 2.63
const FCBG_TH_MBQKG = 1.46

@testset "PAnnularDisk constructor + invariants" begin
    p = PAnnularDisk(R_IN_CM, R_OUT_CM, 0.0, 0.78, +1, mat_SS,
                     FCBG_BI_MBQKG, FCBG_TH_MBQKG; name="FCBG_test")
    @test p.R_in   == R_IN_CM
    @test p.R_out  == R_OUT_CM
    @test p.z_face == 0.0
    @test p.H_axial == 0.78
    @test p.normal_sign == +1
    @test p.count == 1
    @test p.name  == "FCBG_test"

    # Bad inputs
    @test_throws ErrorException PAnnularDisk(-1.0, 80.3, 0.0, 1.0, +1, mat_SS, 1.0, 1.0)
    @test_throws ErrorException PAnnularDisk(80.0, 72.8, 0.0, 1.0, +1, mat_SS, 1.0, 1.0)
    @test_throws ErrorException PAnnularDisk(72.8, 80.3, 0.0, -0.5, +1, mat_SS, 1.0, 1.0)
    @test_throws ErrorException PAnnularDisk(72.8, 80.3, 0.0, 1.0, 0,  mat_SS, 1.0, 1.0)
end

@testset "Volume and mass" begin
    # Volume = π·(R_out² − R_in²)·H_axial = π·(80.3² − 72.8²)·H
    p = PAnnularDisk(R_IN_CM, R_OUT_CM, 0.0, 0.78, +1, mat_SS,
                     FCBG_BI_MBQKG, FCBG_TH_MBQKG)
    expected_V = π * (R_OUT_CM^2 - R_IN_CM^2) * 0.78
    @test isapprox(volume_shell(p), expected_V; rtol=1e-12)

    # Mass = V × ρ × count / 1000  (cm³ × g/cm³ → g → kg)
    expected_mass_kg = expected_V * mat_SS.density / 1000.0
    @test isapprox(mass(p), expected_mass_kg; rtol=1e-12)

    # Pick H_axial so mass equals the FCBG target — verifies the
    # mass-from-density derivation we did in the design.
    H_for_22kg = FCBG_MASS_KG * 1000.0 /
                 (mat_SS.density * π * (R_OUT_CM^2 - R_IN_CM^2))
    p2 = PAnnularDisk(R_IN_CM, R_OUT_CM, 0.0, H_for_22kg, +1, mat_SS,
                      FCBG_BI_MBQKG, FCBG_TH_MBQKG)
    @test isapprox(mass(p2), FCBG_MASS_KG; rtol=1e-10)
    @test isapprox(H_for_22kg, 0.779; rtol=1e-3)        # design value 0.78 cm
end

@testset "count multiplicity" begin
    p1 = PAnnularDisk(R_IN_CM, R_OUT_CM, 0.0, 0.78, +1, mat_SS,
                     FCBG_BI_MBQKG, FCBG_TH_MBQKG; count=1)
    p3 = PAnnularDisk(R_IN_CM, R_OUT_CM, 0.0, 0.78, +1, mat_SS,
                     FCBG_BI_MBQKG, FCBG_TH_MBQKG; count=3)
    @test mass(p3) ≈ 3 * mass(p1)
    @test activity_U238_late(p3) ≈ 3 * activity_U238_late(p1)
end

@testset "Activities and γ rates" begin
    H_for_22kg = FCBG_MASS_KG * 1000.0 /
                 (mat_SS.density * π * (R_OUT_CM^2 - R_IN_CM^2))
    p = PAnnularDisk(R_IN_CM, R_OUT_CM, 0.0, H_for_22kg, +1, mat_SS,
                     FCBG_BI_MBQKG, FCBG_TH_MBQKG; name="FCBG")

    # 22.28 kg × 2.63 mBq/kg = 58.6 mBq (Bi-214 chain)
    @test isapprox(activity_U238_late(p),  FCBG_MASS_KG * FCBG_BI_MBQKG; rtol=1e-10)
    # 22.28 kg × 1.46 mBq/kg = 32.5 mBq (Tl-208 chain)
    @test isapprox(activity_Th232_late(p), FCBG_MASS_KG * FCBG_TH_MBQKG; rtol=1e-10)

    # γ-production rates use the standard chain BR formulas
    γ_Bi  = activity_U238_late(p)  * 1e-3 * BR_BI214_GAMMA      * SEC_PER_YEAR
    γ_Tl  = activity_Th232_late(p) * 1e-3 * BR_TL208_FROM_CHAIN * SEC_PER_YEAR
    γ_Tlc = activity_Th232_late(p) * 1e-3 * BR_TL208_FROM_CHAIN *
            BR_TL208_COMPANION * SEC_PER_YEAR
    @test isapprox(gamma_rate_Bi214(p),            γ_Bi;  rtol=1e-12)
    @test isapprox(gamma_rate_Tl208(p),            γ_Tl;  rtol=1e-12)
    @test isapprox(gamma_rate_Tl208_companion(p),  γ_Tlc; rtol=1e-12)
end

@testset "Slab thickness dispatch" begin
    p_thick = PAnnularDisk(R_IN_CM, R_OUT_CM, 0.0, 1.56, -1, mat_SS, 2.63, 1.46)
    p_thin  = PAnnularDisk(R_IN_CM, R_OUT_CM, 0.0, 0.78, +1, mat_SS, 2.63, 1.46)
    @test source_slab_thickness(p_thick) ≈ 1.56
    @test source_slab_thickness(p_thin)  ≈ 0.78
end

@testset "Self-shielded spectrum: thin-slab limit" begin
    # In the limit μ·t → 0 the slab integral reduces to dN/du = R/2
    # (constant over u), matching the surface source case.
    H_tiny = 1e-6
    p = PAnnularDisk(R_IN_CM, R_OUT_CM, 0.0, H_tiny, +1, mat_SS, 2.63, 1.46)
    R = gamma_rate_Bi214(p)
    u_bins = collect(range(0.005, 0.995, length=50))
    dNdu = self_shielded_spectrum(p, R, E_BI214_MEV, u_bins)
    # All bins should be ≈ R/2 at this tiny thickness
    @test all(isapprox(d, R/2; rtol=1e-3) for d in dNdu)
end

@testset "Self-shielded spectrum: thick-slab regime" begin
    # μ·t ~ 0.26 (FCBG axial) → significant absorption, monotone in u
    H_FCBG = 0.78
    p = PAnnularDisk(R_IN_CM, R_OUT_CM, 0.0, H_FCBG, +1, mat_SS,
                     FCBG_BI_MBQKG, FCBG_TH_MBQKG)
    R = gamma_rate_Bi214(p)
    u_bins = collect(range(0.005, 0.995, length=100))
    dNdu = self_shielded_spectrum(p, R, E_BI214_MEV, u_bins)

    # Monotone increasing in u (less grazing → less self-absorption)
    @test all(diff(dNdu) .≥ -1e-12)

    # Total inward flux: trapezoidal integral. Must be < R/2 (the
    # no-shielding limit) and > 0.
    total_inward = sum(0.5 * (dNdu[i] + dNdu[i+1]) * (u_bins[i+1] - u_bins[i])
                       for i in 1:length(u_bins)-1)
    @test 0 < total_inward < R/2

    # At μ_SS · H ≈ 0.24 (Fe proxy at 2.448 MeV) the inward fraction
    # sits around 0.35 × R (vs. the no-shielding limit of 0.5 × R) —
    # bound checked against an offline numerical integral.
    inward_fraction = total_inward / R
    @test 0.32 < inward_fraction < 0.40

    println()
    println("── PAnnularDisk slab integral sanity (FCBG-like) ──")
    @printf("  H_axial               %.3f cm\n", H_FCBG)
    @printf("  μ_SS at 2.448 MeV     %.4f cm⁻¹\n", mat_SS.μ_lin(E_BI214_MEV))
    @printf("  μ·H                   %.4f\n", mat_SS.μ_lin(E_BI214_MEV) * H_FCBG)
    @printf("  inward fraction       %.4f  (no-shielding limit = 0.5)\n",
            inward_fraction)
    println()
end

@testset "Self-shielded spectrum agrees with PDisk in matching limit" begin
    # PDisk and PAnnularDisk share the same slab integral. With matching
    # H_axial / wall_thickness and the same material, dN/du(u) curves
    # should be identical (the formula does not depend on lateral footprint).
    H = 0.50
    p_anndisk = PAnnularDisk(R_IN_CM, R_OUT_CM, 0.0, H, +1, mat_SS, 2.63, 1.46)
    g_disc    = GDisk(R_OUT_CM, H, 0.0, Inf, :up)            # flat disc
    p_disc    = PDisk(g_disc, mat_SS, 2.63, 1.46; name="ref_disc")

    u_bins = collect(range(0.01, 0.99, length=20))
    R = 1.0    # γ/yr — only the angular shape matters here
    dNdu_a = self_shielded_spectrum(p_anndisk, R, E_BI214_MEV, u_bins)
    dNdu_d = self_shielded_spectrum(p_disc,    R, E_BI214_MEV, u_bins)
    # Same integrand (same μ, same t, same R) → identical arrays
    @test all(isapprox.(dNdu_a, dNdu_d; rtol=1e-12))
end

@testset "PTFE material loads (sanity for FCPT in commit B)" begin
    # μ_lin should give a plausible value at our energies; this is a
    # smoke test for the table loader, not a precision check.
    @test 0.05 < mat_PTFE.μ_lin(E_BI214_MEV) < 0.15
    @test 0.05 < mat_PTFE.μ_lin(E_TL208_MEV) < 0.15
end
