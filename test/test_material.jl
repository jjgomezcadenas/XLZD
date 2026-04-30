# test/test_material.jl — Verify the Ti material loader and μ_lin(E) interp.

using Test
using Printf

include("../src2/XLZD2.jl")
using .XLZD2

ti_path = joinpath(@__DIR__, "..", "data", "nist_ti.csv")
@assert isfile(ti_path)

mat_Ti = load_material("Ti", 4.510, ti_path)

@testset "Material(Ti) basics" begin
    @test mat_Ti.name == "Ti"
    @test mat_Ti.density ≈ 4.510
    @test mat_Ti.μ_lin isa Function
    @test mat_Ti.μ_lin(2.0) > 0.0
    @test mat_Ti.μ_lin(3.0) > 0.0
end

@testset "μ_lin log-log interpolation" begin
    # Tabulated values in data/nist_ti.csv:
    #   2.000 MeV → 4.180e-2 cm²/g
    #   2.044 MeV → 4.139e-2 cm²/g
    #   3.000 MeV → 3.512e-2 cm²/g
    # μ_lin = (mass attn) × density.
    @test mat_Ti.μ_lin(2.000) ≈ 4.180e-2 * 4.510 rtol=1e-6
    @test mat_Ti.μ_lin(3.000) ≈ 3.512e-2 * 4.510 rtol=1e-6

    # Bi-214 line at 2.448 MeV: log-log between 2.044 and 3.000
    μ_Ti_Bi214 = mat_Ti.μ_lin(2.448)
    @test 0.15 < μ_Ti_Bi214 < 0.20            # cm⁻¹, sanity
    @test μ_Ti_Bi214 ≈ 0.173 atol=0.01

    # Tl-208 line at 2.615 MeV
    μ_Ti_Tl208 = mat_Ti.μ_lin(2.615)
    @test 0.15 < μ_Ti_Tl208 < 0.20
    @test μ_Ti_Tl208 ≈ 0.168 atol=0.01

    # Monotonic decrease from 2 to 3 MeV (Compton-dominated regime)
    @test mat_Ti.μ_lin(2.0) > mat_Ti.μ_lin(2.5) > mat_Ti.μ_lin(3.0)
end

@testset "Out-of-range queries error" begin
    @test_throws ErrorException mat_Ti.μ_lin(0.05)
    @test_throws ErrorException mat_Ti.μ_lin(50.0)
end

println()
@printf("  μ_Ti(0.583 MeV) = %.4f cm⁻¹  (Tl-208 companion)\n", mat_Ti.μ_lin(0.583))
@printf("  μ_Ti(2.448 MeV) = %.4f cm⁻¹  (Bi-214)\n",            mat_Ti.μ_lin(2.448))
@printf("  μ_Ti(2.615 MeV) = %.4f cm⁻¹  (Tl-208)\n",            mat_Ti.μ_lin(2.615))
println()
