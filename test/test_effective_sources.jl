# test/test_effective_sources.jl — Build the 6 LZ EffectiveSource objects
# (CB, CTH, CBH × Bi-214, Tl-208) and verify aggregation and attenuation.

using Test
using Printf

include("../src2/XLZD2.jl")
using .XLZD2

ti_path        = joinpath(@__DIR__, "..", "data", "nist_ti.csv")
geom_csv       = joinpath(@__DIR__, "..", "data", "lz_cryo_geometry.csv")
extras_csv     = joinpath(@__DIR__, "..", "data", "lz_cryo_extras.csv")
surfaces_csv   = joinpath(@__DIR__, "..", "data", "lz_cryo_surface_sources.csv")

mat_Ti = load_material("Ti", 4.510, ti_path)
cryo   = build_cryostat(geom_csv, extras_csv, surfaces_csv)
indiv  = build_individual_sources(cryo, mat_Ti)
effs   = build_effective_sources(indiv, cryo, mat_Ti)

# ---------------------------------------------------------------------------

# Trapezoidal helper for test-side spectrum integral checks
_trapz_t(y, x) = sum(0.5 * (y[i]+y[i+1]) * (x[i+1]-x[i]) for i in 1:length(x)-1)

# Pull individual sources by name
indiv_by_name = Dict(s.name => s for s in indiv)

@testset "Six effective sources built" begin
    @test length(effs) == 6
    names = sort([(e.name, e.isotope) for e in effs])
    @test names == [
        ("CB",  :Bi214), ("CB",  :Tl208),
        ("CBH", :Bi214), ("CBH", :Tl208),
        ("CTH", :Bi214), ("CTH", :Tl208),
    ]

    for e in effs
        @test e.region in (:barrel, :endcap_top, :endcap_bottom)
        @test e.E_MeV in (E_BI214_MEV, E_TL208_MEV)
        @test e.total_per_yr > 0
        # Spectrum integrates to total_per_yr
        @test isapprox(_trapz_t(e.dNdu, e.u_bins), e.total_per_yr; rtol=1e-10)
    end
end

@testset "Contribution counts per effective source (with MLI)" begin
    # MLI surface source attached to ICV_body adds one CB contribution
    by = Dict((e.name, e.isotope) => e for e in effs)
    @test length(by[("CB",  :Bi214)].contributions) == 8   # 7 Ti + 1 MLI
    @test length(by[("CTH", :Bi214)].contributions) == 6
    @test length(by[("CBH", :Bi214)].contributions) == 4
    # 8 + 6 + 4 = 18 individual sources accounted per isotope (17 Ti + 1 MLI)
    @test (length(by[("CB",  :Bi214)].contributions) +
           length(by[("CTH", :Bi214)].contributions) +
           length(by[("CBH", :Bi214)].contributions)) == 18

    # The MLI contribution sits in CB with one ICV-body slab downstream
    cb = by[("CB", :Bi214)]
    mli_contribs = [c for c in cb.contributions if c.source.producer isa PSurface]
    @test length(mli_contribs) == 1
    @test mli_contribs[1].source.producer.name == "Cryostat_insulation_MLI"
    @test length(mli_contribs[1].downstream_slabs) == 1
end

@testset "Downstream attenuation never increases flux" begin
    # For each contribution, attenuated dNdu(u) ≤ source.dNdu(u) for all u.
    for e in effs
        for c in e.contributions
            T = transmission_factor(c.downstream_slabs, c.source.u_bins)
            atten = c.source.dNdu .* T
            @test all(atten .<= c.source.dNdu .+ 1e-15)
        end
    end
end

@testset "Effective total ≤ sum of contribution exit_inward" begin
    # Each effective source's total must be ≤ Σ exit_inward over its
    # contributions (downstream attenuation only reduces flux).
    for e in effs
        upper = sum(c.source.exit_inward_per_yr for c in e.contributions)
        @test e.total_per_yr <= upper + 1e-9
    end
end

@testset "ICV-only contributions (no downstream) pass through unchanged" begin
    # CB Bi-214 contains ICV_barrel_Bi214 and ICV_main_flange_pair_Bi214 with
    # empty downstream lists; their attenuated dNdu must equal their raw dNdu.
    cb = first(e for e in effs if e.name=="CB" && e.isotope===:Bi214)
    icv_barrel  = first(c for c in cb.contributions
                          if c.source.name == "ICV_barrel_Bi214")
    icv_flange  = first(c for c in cb.contributions
                          if c.source.name == "ICV_main_flange_pair_Bi214")
    @test isempty(icv_barrel.downstream_slabs)
    @test isempty(icv_flange.downstream_slabs)
    T = transmission_factor(icv_barrel.downstream_slabs, icv_barrel.source.u_bins)
    @test all(T .== 1.0)
end

@testset "Documentation: effective source rates" begin
    println()
    println("── Effective source γ rates at ICV inner surfaces ──")
    @printf("  %-4s %-6s %12s %12s   contribs\n",
            "name", "iso", "total/yr", "Σ exit_in")
    println("  ", "─"^60)
    for e in effs
        upper = sum(c.source.exit_inward_per_yr for c in e.contributions)
        @printf("  %-4s %-6s %12.3e %12.3e   %d\n",
                e.name, string(e.isotope), e.total_per_yr, upper,
                length(e.contributions))
    end
    println()

    # Sanity: total Bi-214 across the three regions vs the sum of all
    # 17 individual sources' exit_inward.
    total_eff_Bi = sum(e.total_per_yr for e in effs if e.isotope===:Bi214)
    total_ind_Bi = sum(s.exit_inward_per_yr for s in indiv if s.isotope===:Bi214)
    @printf("  Total Bi-214 effective γ/yr = %.3e\n", total_eff_Bi)
    @printf("  Total Bi-214 individual exit_inward γ/yr = %.3e\n", total_ind_Bi)
    @printf("  ratio (effective/individual) = %.3f\n",
            total_eff_Bi / total_ind_Bi)
    println()

    # Effective should be smaller (downstream attenuation), so ratio < 1.
    @test total_eff_Bi < total_ind_Bi
    # And not too small (most sources are ICV-attached or have only one
    # ICV slab downstream): ratio > 0.4.
    @test total_eff_Bi / total_ind_Bi > 0.4
end

@testset "Sanity vs lz_cryo_linear_fit.csv (current reference)" begin
    # The Python pre-processing in py/lz_cryo.py produces
    # data/lz_cryo_linear_fit.csv with one row per (isotope, region).
    # We compare our effective-source totals to those reference values.
    # NOTE: not strict parity — our model includes flanges and main ports
    # absent from lz_cryo.py, so we expect comparable order of magnitude.
    by = Dict((e.name, e.isotope) => e.total_per_yr for e in effs)
    println("── Comparison vs Python reference (lz_cryo_linear_fit.csv) ──")
    println("  Region   Iso    Julia (γ/yr)   Python ref (γ/yr)")
    # Reference values from data/lz_cryo_linear_fit.csv (existing file)
    ref = Dict(
        ("CB",  :Bi214) => 3.382267e4,
        ("CB",  :Tl208) => 1.046081e6,
        ("CTH", :Bi214) => 1.493e4,    # endcap is total; we have CTH+CBH
        ("CTH", :Tl208) => 4.80e5,
        ("CBH", :Bi214) => 1.493e4,    # treat as half each by lz_cryo.py
        ("CBH", :Tl208) => 4.80e5,
    )
    for (key, jul) in by
        @printf("  %-7s %-6s %12.3e   %12.3e\n",
                key[1], string(key[2]), jul, get(ref, key, NaN))
    end
    println()

    # With MLI added, the barrel rates should be close to the Python
    # reference (within ~20%). The endcap rates remain ~half the Python
    # reference because Python's 1160 kg "extra Ti" surface source on
    # the bottom endcap is distributed in our model across explicit
    # geometric elements (legs, flanges, anchors) — most of which are
    # MC-inactive or end up CB-attached rather than head-attached.
    @test 0.7 < by[("CB", :Bi214)] / ref[("CB", :Bi214)] < 1.3
    @test 0.7 < by[("CB", :Tl208)] / ref[("CB", :Tl208)] < 1.3
    # Endcap residual: still within a factor of 5 (sanity-only)
    for region in ("CTH", "CBH")
        for iso in (:Bi214, :Tl208)
            jul = by[(region, iso)]
            @test 0.1 < jul / ref[(region, iso)] < 5.0
        end
    end
end
