# scripts/cryo_source.jl — Compute the dN/du angular spectra for the LZ
# cryostat effective γ sources at the ICV inner exit surfaces, and plot
# Bi-214, Tl-208 main and Tl-208 cascade-companion γ on the same axes
# for each of the three regions (CB, CTH, CBH).
#
# Run from the project root:
#     julia --project=. scripts/cryo_source.jl

using Plots
using Printf

include(joinpath(@__DIR__, "..", "src2", "XLZD2.jl"))
using .XLZD2

const PROJECT_ROOT = abspath(joinpath(@__DIR__, ".."))
const OUT_DIR      = joinpath(PROJECT_ROOT, "output", "cryo_sources")

const TI_PATH      = joinpath(PROJECT_ROOT, "data", "nist_ti.csv")
const GEOM_CSV     = joinpath(PROJECT_ROOT, "data", "lz_cryo_geometry.csv")
const EXTRAS_CSV   = joinpath(PROJECT_ROOT, "data", "lz_cryo_extras.csv")
const SURFACES_CSV = joinpath(PROJECT_ROOT, "data", "lz_cryo_surface_sources.csv")

# Discrete trapezoidal integral (matches the helper in sources.jl)
_trapz(y, x) = sum(0.5 * (y[i]+y[i+1]) * (x[i+1]-x[i]) for i in 1:length(x)-1)

# ---------------------------------------------------------------------------
# Load + build
# ---------------------------------------------------------------------------
mat_Ti = load_material("Ti", 4.510, TI_PATH)
cryo   = build_cryostat(GEOM_CSV, EXTRAS_CSV, SURFACES_CSV)
indiv  = build_individual_sources(cryo, mat_Ti)
effs   = build_effective_sources(indiv, cryo, mat_Ti)
by     = Dict(e.name => e for e in effs)

mkpath(OUT_DIR)

# ---------------------------------------------------------------------------
# Per-region plots: three isotopes overlaid
# ---------------------------------------------------------------------------
const REGION_TITLES = Dict(
    "CB"  => "CB — Cryostat Barrel (ICV inner cylindrical wall)",
    "CTH" => "CTH — Cryostat Top Head (ICV inner top head)",
    "CBH" => "CBH — Cryostat Bottom Head (ICV inner bottom head)",
)

for region in ("CB", "CTH", "CBH")
    e1 = by["$(region)_Bi214"]
    e2 = by["$(region)_Tl208"]
    e3 = by["$(region)_Tl208c"]

    plt = plot(
        xlabel = "u = cos θ",
        ylabel = "dN/du  (γ/yr per unit u)",
        yscale = :log10,
        title  = REGION_TITLES[region],
        legend = :bottomright,
        size   = (800, 500),
        dpi    = 150,
    )
    plot!(plt, e1.u_bins, e1.dNdu,
          label = "²¹⁴Bi  (2.448 MeV)", lw=2, color=:steelblue)
    plot!(plt, e2.u_bins, e2.dNdu,
          label = "²⁰⁸Tl  (2.615 MeV)", lw=2, color=:firebrick)
    plot!(plt, e3.u_bins, e3.dNdu,
          label = "²⁰⁸Tl companion (0.583 MeV)", lw=2, color=:darkorange,
          linestyle=:dash)

    fpath = joinpath(OUT_DIR, "$(region).png")
    savefig(plt, fpath)
    println("  wrote ", fpath)
end

# ---------------------------------------------------------------------------
# 3-panel composite
# ---------------------------------------------------------------------------
combined = plot(layout=(1,3), size=(1500, 450), dpi=150,
                left_margin=4Plots.mm, bottom_margin=4Plots.mm)
for (i, region) in enumerate(("CB", "CTH", "CBH"))
    e1 = by["$(region)_Bi214"]
    e2 = by["$(region)_Tl208"]
    e3 = by["$(region)_Tl208c"]
    plot!(combined, e1.u_bins, e1.dNdu, subplot=i, yscale=:log10,
          label="²¹⁴Bi",       lw=2, color=:steelblue,
          xlabel="u = cos θ", ylabel=(i==1 ? "dN/du (γ/yr/u)" : ""),
          title=region, legend=(i==3 ? :bottomright : false))
    plot!(combined, e2.u_bins, e2.dNdu, subplot=i,
          label="²⁰⁸Tl",       lw=2, color=:firebrick)
    plot!(combined, e3.u_bins, e3.dNdu, subplot=i,
          label="²⁰⁸Tl comp.", lw=2, color=:darkorange, linestyle=:dash)
end
fpath = joinpath(OUT_DIR, "all_regions.png")
savefig(combined, fpath)
println("  wrote ", fpath)

# ---------------------------------------------------------------------------
# Summary table to stdout
# ---------------------------------------------------------------------------
println()
println("── Effective source summary (at ICV inner exit surface) ──")
@printf("  %-12s %-7s %8s %14s\n", "name", "E (MeV)", "<u>", "total γ/yr")
println("  ", "─"^48)
for region in ("CB", "CTH", "CBH")
    for iso_suffix in ("Bi214", "Tl208", "Tl208c")
        e = by["$(region)_$iso_suffix"]
        # Mean u weighted by dN/du (trapezoidal integrals)
        num   = _trapz(e.u_bins .* e.dNdu, e.u_bins)
        denom = _trapz(e.dNdu, e.u_bins)
        mean_u = denom > 0 ? num / denom : NaN
        @printf("  %-12s %7.3f %8.3f %14.3e\n",
                e.name, e.E_MeV, mean_u, e.total_per_yr)
    end
end
println()
