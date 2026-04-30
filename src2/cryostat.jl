# src2/cryostat.jl — Composite Cryostat object built from CSV data files.
#
# Two input files:
#   data/lz_cryo_geometry.csv  → barrels (OCV, ICV) and four ellipsoidal heads
#   data/lz_cryo_extras.csv    → flanges, ports, legs, stiffeners, fins,
#                                internal anchor rings (each row is a GCyl).
#
# The full cryostat is the union of these objects. `include_in_mc` flags
# which extras are exposed as Monte Carlo source rings; barrels and heads
# are always MC-active.

# ---------------------------------------------------------------------------
# Types
# ---------------------------------------------------------------------------

"""
    CryostatExtra

One row of `lz_cryo_extras.csv`: an annular Ti shell with a category label
and a flag for whether it serves as an MC source. `count` lets the same
geometry represent N identical pieces (e.g. 6 coldhead fins, 3 tie-bar
ports), with the total mass scaled by `count`.
"""
struct CryostatExtra
    name::String
    category::Symbol
    count::Int
    include_in_mc::Bool
    shell::GCyl
end

"""
    Cryostat

Full LZ cryostat built from the two CSVs:
  * `barrels`  — OCV and ICV cylindrical bodies
  * `heads`    — four ellipsoidal heads (OCV top/bot, ICV top/bot)
  * `extras`   — every additional Ti element (flanges, ports, …)
"""
struct Cryostat
    barrels::Vector{GCyl}
    heads::Vector{GDisk}
    extras::Vector{CryostatExtra}
end

# ---------------------------------------------------------------------------
# CSV parsing helpers
# ---------------------------------------------------------------------------

function _read_csv_rows(path::AbstractString)
    raw, hdr = readdlm(path, ',', Any, '\n', header=true)
    header = vec(strip.(string.(hdr)))
    rows = Vector{Dict{String,Any}}()
    for i in axes(raw, 1)
        d = Dict{String,Any}()
        for j in 1:length(header)
            d[header[j]] = raw[i, j]
        end
        push!(rows, d)
    end
    rows
end

function _parse_bool(v)::Bool
    s = lowercase(strip(string(v)))
    s in ("true", "1", "yes", "t")  && return true
    s in ("false", "0", "no", "f")  && return false
    error("Cannot parse '$v' as Bool")
end

# ---------------------------------------------------------------------------
# Builder
# ---------------------------------------------------------------------------

"""
    build_cryostat(geom_csv, extras_csv) -> Cryostat

Construct the LZ cryostat object from `data/lz_cryo_geometry.csv` and
`data/lz_cryo_extras.csv`. Specific axial positions in the geometry CSV
do not affect mass; they are placeholders used by the radioactivity step
later.
"""
function build_cryostat(geom_csv::AbstractString, extras_csv::AbstractString)
    geom_rows = _read_csv_rows(geom_csv)
    barrels = GCyl[]
    heads   = GDisk[]
    for row in geom_rows
        R   = Float64(row["R_cm"])
        H_c = Float64(row["H_cyl_cm"])
        t_w = Float64(row["t_wall_mm"]) / 10.0
        t_t = Float64(row["t_top_mm"])  / 10.0
        t_b = Float64(row["t_bot_mm"])  / 10.0
        n_t = Float64(row["top_ar"])
        n_b = Float64(row["bot_ar"])

        z_eq_bot = 0.0
        z_eq_top = H_c

        push!(barrels, GCyl(R, t_w, z_eq_bot, z_eq_top))
        push!(heads,   GDisk(R, t_t, z_eq_top, n_t, :up))
        push!(heads,   GDisk(R, t_b, z_eq_bot, n_b, :down))
    end

    extras_rows = _read_csv_rows(extras_csv)
    extras = CryostatExtra[]
    for row in extras_rows
        name     = strip(string(row["name"]))
        category = Symbol(strip(string(row["category"])))
        R_i      = Float64(row["R_inner_cm"])
        R_o      = Float64(row["R_outer_cm"])
        z_min    = Float64(row["z_min_cm"])
        z_max    = Float64(row["z_max_cm"])
        cnt      = Int(row["count"])
        in_mc    = _parse_bool(row["include_in_mc"])

        shell = GCyl(R_i, R_o - R_i, z_min, z_max)
        push!(extras, CryostatExtra(name, category, cnt, in_mc, shell))
    end

    Cryostat(barrels, heads, extras)
end

# ---------------------------------------------------------------------------
# Mass aggregation
# ---------------------------------------------------------------------------

"""
    total_mass(c::Cryostat, ρ) -> Float64

Total Ti mass (kg) summed over barrels, heads, and all extras (active and
inactive). Each extra contributes `mass(shell, ρ) × count`.
"""
function total_mass(c::Cryostat, ρ::Real)::Float64
    m  = sum(mass(g, ρ) for g in c.barrels; init=0.0)
    m += sum(mass(h, ρ) for h in c.heads;   init=0.0)
    m += sum(mass(e.shell, ρ) * e.count for e in c.extras; init=0.0)
    m
end

"""
    mass_breakdown(c::Cryostat, ρ) -> Dict{Symbol, Float64}

Mass (kg) per category. Keys: `:barrel`, `:head`, plus every
category appearing in the extras (`:flange`, `:port`, `:leg`,
`:stiffener`, `:fin`, `:internal_anchor`, ...).
"""
function mass_breakdown(c::Cryostat, ρ::Real)::Dict{Symbol, Float64}
    out = Dict{Symbol, Float64}()
    out[:barrel] = sum(mass(g, ρ) for g in c.barrels; init=0.0)
    out[:head]   = sum(mass(h, ρ) for h in c.heads;   init=0.0)
    for e in c.extras
        out[e.category] = get(out, e.category, 0.0) + mass(e.shell, ρ) * e.count
    end
    out
end

"""
    mc_active_extras(c::Cryostat) -> Vector{CryostatExtra}

Subset of extras flagged `include_in_mc = true`.
"""
mc_active_extras(c::Cryostat)::Vector{CryostatExtra} =
    filter(e -> e.include_in_mc, c.extras)

"""
    mc_active_mass(c::Cryostat, ρ) -> Float64

Ti mass (kg) currently flagged as MC-active: barrels + heads + active extras.
"""
function mc_active_mass(c::Cryostat, ρ::Real)::Float64
    m  = sum(mass(g, ρ) for g in c.barrels; init=0.0)
    m += sum(mass(h, ρ) for h in c.heads;   init=0.0)
    m += sum(mass(e.shell, ρ) * e.count for e in mc_active_extras(c); init=0.0)
    m
end
