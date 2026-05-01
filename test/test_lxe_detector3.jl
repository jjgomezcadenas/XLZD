# test/test_lxe_detector3.jl — Build the LXeDetector and verify region_at.

using Test
using Printf

isdefined(Main, :XLZD3) || include("../src3/XLZD3.jl")
using .XLZD3

const lxe_csv  = joinpath(@__DIR__, "..", "data", "lxe_detector.csv")
const lxe_nist = joinpath(@__DIR__, "..", "data", "nist_lxe.csv")
@assert isfile(lxe_csv)
@assert isfile(lxe_nist)

mat_LXe = load_material("LXe", 2.953, lxe_nist)
det     = build_lxe_detector(lxe_csv, mat_LXe)

@testset "Build + parameters" begin
    @test det.R_FC_inner   ≈ 72.8
    @test det.R_FC_outer   ≈ 74.3
    @test det.R_ICV_inner  ≈ 82.1
    @test det.z_cathode    ≈ 0.0
    @test det.z_gate       ≈ 145.6
    @test det.z_RFR_bottom ≈ -13.75
    @test det.z_LXe_bottom ≈ -69.0
    @test det.z_ICV_top    ≈ 190.0
    @test det.E_visible_keV   ≈ 10.0
    @test det.E_skin_veto_keV ≈ 100.0
    @test det.material.name == "LXe"
end

@testset "region_at — interior representatives" begin
    # On-axis active drift region
    @test region_at(det, 0.0,  0.0, 50.0)  === :active
    @test region_at(det, 0.0,  0.0, 100.0) === :active
    # Off-axis but inside FC inner radius
    @test region_at(det, 50.0, 0.0, 50.0)  === :active
    @test region_at(det, 0.0, 50.0, 50.0)  === :active

    # FC annulus (transparent)
    @test region_at(det, 73.5, 0.0,  50.0) === :fc_region
    @test region_at(det, 0.0,  73.5, 100.0) === :fc_region
    # FC annulus also exists in the RFR z range
    @test region_at(det, 73.5, 0.0, -10.0) === :fc_region

    # Skin (between FC outer and ICV inner, in the FC-bearing z range)
    @test region_at(det, 78.0, 0.0,  50.0) === :skin
    @test region_at(det, 80.0, 0.0,  100.0) === :skin
    @test region_at(det, 78.0, 0.0, -10.0) === :skin   # skin extends through RFR z

    # RFR + dome inert LXe
    @test region_at(det, 0.0, 0.0, -5.0)  === :inert   # RFR (below cathode, above z_RFR_bottom)
    @test region_at(det, 0.0, 0.0, -30.0) === :inert   # dome (below z_RFR_bottom)
    @test region_at(det, 50.0, 0.0, -50.0) === :inert  # dome
    @test region_at(det, 78.0, 0.0, -30.0) === :inert  # below z_RFR_bottom — no skin here

    # Gas above the liquid
    @test region_at(det, 0.0,  0.0, 150.0) === :gas
    @test region_at(det, 0.0,  0.0, 180.0) === :gas
    @test region_at(det, 50.0, 0.0, 160.0) === :gas

    # Outside the LXe (radially in the ICV wall or beyond)
    @test region_at(det, 85.0, 0.0,  50.0) === :outside_lxe
    @test region_at(det, 0.0, 90.0, -20.0) === :outside_lxe
end

@testset "region_at — small jitters around boundaries" begin
    eps_r = 1e-6
    eps_z = 1e-6

    # r = R_FC_inner ± ε at active z
    @test region_at(det, det.R_FC_inner - eps_r, 0.0, 50.0) === :active
    @test region_at(det, det.R_FC_inner + eps_r, 0.0, 50.0) === :fc_region

    # r = R_FC_outer ± ε
    @test region_at(det, det.R_FC_outer - eps_r, 0.0, 50.0) === :fc_region
    @test region_at(det, det.R_FC_outer + eps_r, 0.0, 50.0) === :skin

    # r = R_ICV_inner ± ε
    @test region_at(det, det.R_ICV_inner - eps_r, 0.0, 50.0) === :skin
    @test region_at(det, det.R_ICV_inner + eps_r, 0.0, 50.0) === :outside_lxe

    # z just above gate → gas
    @test region_at(det, 0.0, 0.0, det.z_gate + eps_z) === :gas
    # z just below z_RFR_bottom (and above the cathode is impossible — z_RFR_bottom < z_cathode):
    @test region_at(det, 0.0, 0.0, det.z_RFR_bottom - eps_z) === :inert
end

@testset "Sanity: μ_LXe and masses" begin
    # μ_LXe at the relevant γ energies
    μ_2448 = μ_LXe(det, 2.448)
    μ_2615 = μ_LXe(det, 2.615)
    @test 0.07 < μ_2448 < 0.12
    @test 0.07 < μ_2615 < 0.12
    @test μ_2615 < μ_2448   # higher E → lower μ in Compton-dominated regime

    # Active LXe mass: π × 72.8² × 145.6 × 2.953 / 1000 ≈ 7.16 t
    m_active = active_mass_kg(det)
    @test isapprox(m_active, 7160.0; rtol=0.02)

    m_skin = skin_mass_kg(det)
    @test m_skin > 0
end

println()
@printf("  Active LXe mass = %7.1f kg  (target ≈ 7000 kg)\n", active_mass_kg(det))
@printf("  Skin LXe mass   = %7.1f kg\n",                   skin_mass_kg(det))
@printf("  μ_LXe(2.448 MeV) = %.4f cm⁻¹\n",                 μ_LXe(det, 2.448))
@printf("  μ_LXe(2.615 MeV) = %.4f cm⁻¹\n",                 μ_LXe(det, 2.615))
@printf("  μ_LXe(0.583 MeV) = %.4f cm⁻¹  (Tl-208 companion)\n", μ_LXe(det, 0.583))
println()
