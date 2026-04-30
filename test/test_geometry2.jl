# test/test_geometry2.jl — Verify that the new GCyl/GDisk primitives
# reproduce the cryostat masses recorded in data/lz_cryo_geometry.csv,
# and document the resulting Ti mass against the bb0nu paper's 2590 kg.

using Test
using DelimitedFiles
using Random
using Printf

include("../src2/XLZD2.jl")
using .XLZD2

const ρ_TI = 4.510  # g/cm³, Grade-1 titanium

# ---------------------------------------------------------------------------
# CSV loader
# ---------------------------------------------------------------------------
function read_cryo_geometry(path::AbstractString)
    raw, hdr = readdlm(path, ',', Any, '\n', header=true)
    header = vec(strip.(string.(hdr)))
    rows = Dict{String, Dict{String,Float64}}()
    for i in axes(raw, 1)
        elem = strip(string(raw[i, 1]))
        d = Dict{String,Float64}()
        for j in 2:length(header)
            d[header[j]] = Float64(raw[i, j])
        end
        rows[elem] = d
    end
    rows
end

# ---------------------------------------------------------------------------
# Build (barrel, top head, bottom head) from one CSV row
# ---------------------------------------------------------------------------
function build_vessel_objects(row::Dict{String,Float64})
    R   = row["R_cm"]
    H_c = row["H_cyl_cm"]
    t_w = row["t_wall_mm"] / 10.0
    t_t = row["t_top_mm"]  / 10.0
    t_b = row["t_bot_mm"]  / 10.0
    n_t = row["top_ar"]
    n_b = row["bot_ar"]

    # Place barrel between z=0 and z=H_c with bottom head below z=0 and
    # top head above z=H_c. Specific z values do not affect masses.
    z_eq_bot = 0.0
    z_eq_top = H_c

    barrel = GCyl(R, t_w, z_eq_bot, z_eq_top)
    top    = GDisk(R, t_t, z_eq_top, n_t, :up)
    bot    = GDisk(R, t_b, z_eq_bot, n_b, :down)
    (; barrel, top, bot)
end

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
csv_path = joinpath(@__DIR__, "..", "data", "lz_cryo_geometry.csv")
@assert isfile(csv_path) "Cryostat CSV not found at $csv_path"

geom = read_cryo_geometry(csv_path)

@testset "GDisk shape special cases" begin
    flat = GDisk(50.0, 0.5, 0.0, Inf, :up)
    @test depth(flat) == 0.0
    @test is_flat(flat)
    @test area_inner(flat) ≈ π * 50.0^2

    hemi = GDisk(50.0, 0.5, 0.0, 1.0, :up)
    @test depth(hemi) ≈ 50.0
    @test area_inner(hemi) ≈ 2π * 50.0^2

    e21 = GDisk(50.0, 0.5, 0.0, 2.0, :up)
    @test depth(e21) ≈ 25.0
    # Reference value 1.380 π R² for a 2:1 oblate head
    @test area_inner(e21) / (π * 50.0^2) ≈ 1.3802 rtol=1e-3

    e31 = GDisk(50.0, 0.5, 0.0, 3.0, :down)
    @test depth(e31) ≈ 50.0/3
    @test z_apex(e31) ≈ -50.0/3
    @test area_inner(e31) / (π * 50.0^2) ≈ 1.2078 rtol=1e-3
end

@testset "GCyl basic" begin
    g = GCyl(10.0, 0.5, 0.0, 100.0)
    @test R_outer(g) ≈ 10.5
    @test height(g) ≈ 100.0
    @test area_inner(g) ≈ 2π * 10.0 * 100.0
    @test volume_shell(g) ≈ π * (10.5^2 - 10.0^2) * 100.0
    @test mass(g, 4.510) ≈ volume_shell(g) * 4.510 / 1000.0
end

@testset "Cryostat mass from data/lz_cryo_geometry.csv" begin
    masses = Dict{String, Float64}()

    for elem in ("OCV", "ICV")
        row = geom[elem]
        v = build_vessel_objects(row)

        m_barrel = mass(v.barrel, ρ_TI)
        m_top    = mass(v.top,    ρ_TI)
        m_bot    = mass(v.bot,    ρ_TI)
        m_heads  = m_top + m_bot
        m_total  = m_barrel + m_heads

        # CSV uses thin-shell formula 2π·R_inner·H·t·ρ for the cylindrical
        # body, ours uses the exact π(R_o²−R_i²)·H formula. The two differ
        # by a factor 1 + t/(2R) ≈ 0.4–0.5% for LZ vessels. 1% rtol absorbs
        # this and any rounding in the CSV reference values.
        @testset "$elem" begin
            @test isapprox(m_barrel, row["mass_shell_kg"]; rtol=1e-2)
            @test isapprox(m_heads,  row["mass_heads_kg"]; rtol=1e-2)
            @test isapprox(m_total,  row["mass_total_kg"]; rtol=1e-2)
        end

        masses["$(elem)_barrel"] = m_barrel
        masses["$(elem)_top"]    = m_top
        masses["$(elem)_bot"]    = m_bot
        masses["$(elem)_total"]  = m_total
    end

    # Documentation summary
    println()
    println("── Cryostat mass from data/lz_cryo_geometry.csv ──")
    for elem in ("OCV", "ICV")
        row = geom[elem]
        @printf("  %-18s %8.1f kg\n",
                "$elem barrel",     masses["$(elem)_barrel"])
        @printf("  %-18s %8.1f kg   (R=%.1f, t=%.1f mm, %.0f:1)\n",
                "$elem top head",   masses["$(elem)_top"],
                row["R_cm"], row["t_top_mm"], row["top_ar"])
        @printf("  %-18s %8.1f kg   (R=%.1f, t=%.1f mm, %.0f:1)\n",
                "$elem bottom head", masses["$(elem)_bot"],
                row["R_cm"], row["t_bot_mm"], row["bot_ar"])
    end
    total_ti = masses["OCV_total"] + masses["ICV_total"]
    println("  ─────────────────────────────────────")
    @printf("  %-18s %8.1f kg\n", "Total Ti mass",        total_ti)
    @printf("  %-18s %8.1f kg\n", "bb0nu Table I",        2590.0)
    @printf("  %-18s %8.1f kg\n", "Unaccounted hardware", 2590.0 - total_ti)
    println()
end

# ---------------------------------------------------------------------------
# Full cryostat (barrels + heads + extras): mass budget against bb0nu
# ---------------------------------------------------------------------------
extras_path = joinpath(@__DIR__, "..", "data", "lz_cryo_extras.csv")
@assert isfile(extras_path) "Cryostat extras CSV not found at $extras_path"

@testset "Full cryostat Ti mass budget" begin
    cryo = build_cryostat(csv_path, extras_path)
    m_total = total_mass(cryo, ρ_TI)
    bd      = mass_breakdown(cryo, ρ_TI)
    m_mc    = mc_active_mass(cryo, ρ_TI)

    # bb0nu Table I full Ti cryostat assembly: 2590 kg
    @test isapprox(m_total, 2590.0; rtol=0.05)

    # Spot-check that every expected category is present
    for cat in (:barrel, :head, :flange, :port, :leg, :stiffener,
                :fin, :internal_anchor)
        @test haskey(bd, cat)
        @test bd[cat] > 0.0
    end

    # MC-active subset must be a strict subset (≥ barrels+heads, ≤ total)
    @test m_mc <= m_total
    @test m_mc >= bd[:barrel] + bd[:head]

    # Documentation block
    println()
    println("── Full LZ cryostat Ti mass budget ──")
    @printf("  %-30s %8.1f kg\n", "Barrels (OCV + ICV)",     bd[:barrel])
    @printf("  %-30s %8.1f kg\n", "Heads (4 ellipsoidal)",   bd[:head])
    @printf("  %-30s %8.1f kg\n", "Flanges",                 bd[:flange])
    @printf("  %-30s %8.1f kg\n", "Ports",                   bd[:port])
    @printf("  %-30s %8.1f kg\n", "Legs",                    bd[:leg])
    @printf("  %-30s %8.1f kg\n", "Stiffener",               bd[:stiffener])
    @printf("  %-30s %8.1f kg\n", "Coldhead fins",           bd[:fin])
    @printf("  %-30s %8.1f kg\n", "Internal anchor rings",   bd[:internal_anchor])
    println("  ──────────────────────────────────────")
    @printf("  %-30s %8.1f kg\n", "Total Ti mass (model)",   m_total)
    @printf("  %-30s %8.1f kg\n", "bb0nu Table I",           2590.0)
    @printf("  %-30s %+8.1f kg  (%+.1f%%)\n", "Residual (model − bb0nu)",
            m_total - 2590.0, 100.0 * (m_total - 2590.0) / 2590.0)
    println()
    println("── MC-active fraction ──")
    n_active = length(mc_active_extras(cryo))
    n_total  = length(cryo.extras)
    @printf("  %-30s %8.1f kg\n", "Barrels + heads (always)",
            bd[:barrel] + bd[:head])
    @printf("  %-30s %8.1f kg\n", "Active extras",
            m_mc - bd[:barrel] - bd[:head])
    println("  ──────────────────────────────────────")
    @printf("  %-30s %8.1f kg  (%.0f%% of total)\n",
            "MC-active Ti mass", m_mc, 100.0 * m_mc / m_total)
    @printf("  %-30s %8.1f kg  (%.0f%% of total)\n",
            "MC-inactive Ti", m_total - m_mc,
            100.0 * (m_total - m_mc) / m_total)
    @printf("  %-30s %d / %d extras flagged on\n",
            "Active extras count", n_active, n_total)
    println()
end
