# test/test_sources.jl — Build 17 individual GammaSource × 2 isotopes = 34
# self-shielded sources at the cryostat's MC-active Ti volumes.

using Test
using Printf

include("../src2/XLZD2.jl")
using .XLZD2

ti_path     = joinpath(@__DIR__, "..", "data", "nist_ti.csv")
geom_csv    = joinpath(@__DIR__, "..", "data", "lz_cryo_geometry.csv")
extras_csv  = joinpath(@__DIR__, "..", "data", "lz_cryo_extras.csv")

mat_Ti = load_material("Ti", 4.510, ti_path)
cryo   = build_cryostat(geom_csv, extras_csv)

# Trapezoidal helper duplicated for tests
_trapz_t(y, x) = sum(0.5 * (y[i]+y[i+1]) * (x[i+1]-x[i]) for i in 1:length(x)-1)

@testset "build_individual_sources count and shape" begin
    sources = build_individual_sources(cryo, mat_Ti)
    @test length(sources) == 51   # 17 sources × 3 isotopes (Bi214, Tl208, Tl208c)

    # One third each
    @test count(s -> s.isotope === :Bi214,  sources) == 17
    @test count(s -> s.isotope === :Tl208,  sources) == 17
    @test count(s -> s.isotope === :Tl208c, sources) == 17

    # All have positive produced rate
    @test all(s -> s.produced_per_yr > 0, sources)
    # Inward exit ≤ produced/2 (factor 1/2 = inward hemisphere only)
    @test all(s -> s.exit_inward_per_yr ≤ s.produced_per_yr / 2 + 1e-9, sources)

    # Spectrum integrates to exit_inward (within trapezoidal error)
    for s in sources
        I = _trapz_t(s.dNdu, s.u_bins)
        @test isapprox(I, s.exit_inward_per_yr; rtol=1e-10)
    end
end

@testset "Tl-208 companion BR and self-shielding" begin
    sources = build_individual_sources(cryo, mat_Ti)
    by_name = Dict(s.name => s for s in sources)

    # Companion produced rate = Tl-208 produced × BR_TL208_COMPANION
    for base in ("OCV_barrel", "ICV_barrel", "ICV_top_head", "ICV_bottom_head")
        main = by_name["$(base)_Tl208"]
        comp = by_name["$(base)_Tl208c"]
        @test isapprox(comp.produced_per_yr,
                       main.produced_per_yr * BR_TL208_COMPANION;
                       rtol=1e-12)
        @test comp.E_MeV ≈ E_TL208_COMPANION_MEV
        # Companion is more strongly self-shielded than the main γ
        # (μ_Ti(0.583) > μ_Ti(2.615)), so the exit/produced fraction is lower.
        # Surface sources (PSurface) have no self-shielding so frac is constant.
        if !(comp.producer isa PSurface)
            @test comp.exit_inward_per_yr / comp.produced_per_yr <
                  main.exit_inward_per_yr / main.produced_per_yr
        end
    end
end

@testset "Thin-source / thick-source limits" begin
    # Thin source: t → 0, exit_inward/produced → 1/2 over full u ∈ [0,1].
    # Our DEFAULT_U_BINS is [0.005, 0.995] (avoids u=0 divergence) so the
    # truncated integral hits 0.5 × (0.995 − 0.005) = 0.495 in the limit.
    g_thin = GCyl(10.0, 1e-4, 0.0, 100.0)
    p_thin = PCyl(g_thin, mat_Ti, TI_BB0NU_U238_LATE_MBQKG,
                  TI_BB0NU_TH232_LATE_MBQKG; name="thin")
    s_thin = make_gamma_source(p_thin, :Bi214, DEFAULT_U_BINS)
    @test s_thin.exit_inward_per_yr / s_thin.produced_per_yr ≈ 0.495 rtol=1e-3

    # Thick source: t = 10 cm of Ti at 2.448 MeV (μ ≈ 0.17 cm⁻¹ → μt ≈ 1.7),
    # heavy attenuation
    g_thick = GCyl(10.0, 10.0, 0.0, 100.0)
    p_thick = PCyl(g_thick, mat_Ti, TI_BB0NU_U238_LATE_MBQKG,
                   TI_BB0NU_TH232_LATE_MBQKG; name="thick")
    s_thick = make_gamma_source(p_thick, :Bi214, DEFAULT_U_BINS)
    # Analytic depth-averaged transmission for slab (½ ⟨...⟩ ≪ ½)
    @test s_thick.exit_inward_per_yr / s_thick.produced_per_yr < 0.30
    @test s_thick.exit_inward_per_yr / s_thick.produced_per_yr > 0.05
end

@testset "Per-source γ-rate documentation" begin
    sources = build_individual_sources(cryo, mat_Ti)
    println()
    println("── Self-shielded γ rates at source inner surface ──")
    @printf("  %-30s %-7s %12s %12s %6s\n",
            "source", "iso", "produced/yr", "exit_inw/yr", "frac")
    println("  ", "─"^77)
    for s in sources
        frac = s.exit_inward_per_yr / s.produced_per_yr
        @printf("  %-30s %-7s %12.3e %12.3e %6.3f\n",
                s.producer.name, string(s.isotope),
                s.produced_per_yr, s.exit_inward_per_yr, frac)
    end
    println()

    # All MC-active Bi-214 produced summed: should be roughly 84% of the
    # total full-cryostat Bi-214 rate (since we use mc_only=true)
    γ_active = sum(s.produced_per_yr for s in sources if s.isotope === :Bi214)
    # MC-active Ti mass = 2143.7 kg. Bi-214 γ/yr = 2143.7×0.08×1e-3×0.0155×3.1557e7
    # ≈ 8.39e4 γ/yr.
    @test isapprox(γ_active, 8.39e4; rtol=0.05)

    # Total Tl-208 companion produced rate = 0.99 × total Tl-208 produced rate
    γ_main = sum(s.produced_per_yr for s in sources if s.isotope === :Tl208)
    γ_comp = sum(s.produced_per_yr for s in sources if s.isotope === :Tl208c)
    @test isapprox(γ_comp, γ_main * BR_TL208_COMPANION; rtol=1e-12)
end
