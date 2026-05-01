# src3/physics.jl — Photon-physics primitives for the per-photon MC.
#
# Ported from src/physics.jl with no logic change:
#   - XCOMTable: per-component cross sections (cm²/g) for photoelectric,
#     incoherent (Compton), and pair production (nuclear + electronic).
#   - load_xcom: reader for NIST XCOM tables in the 8-column format used by
#     data/nist.csv (Energy, Coherent, Incoherent, Photoelectric,
#     Nuclear PP, Electron PP, Tot. w/ Coherent, Tot. wo/ Coherent).
#   - σ_photo / σ_Compton / σ_pair / μ_total_lin: log-log interpolated
#     accessors plus the linear-attenuation total at a given density.
#   - sample_klein_nishina: Compton scatter via Kahn's rejection method,
#     returning (E_scattered, cos θ).
#   - rotate_direction: rotate a unit vector by polar angle (cos θ) and
#     azimuth φ; numerically stable near the poles.

const ME_C2_MEV = 0.510999  # electron rest mass

"""
    XCOMTable

NIST XCOM data with per-component (photoelectric, Compton, pair) mass
attenuation coefficients (cm²/g) plus their log-log grids for fast
interpolation. One table per material (typically LXe at the relevant
γ energies).
"""
struct XCOMTable
    energy_MeV::Vector{Float64}
    σ_photo::Vector{Float64}
    σ_Compton::Vector{Float64}
    σ_pair::Vector{Float64}
    log_E::Vector{Float64}
    log_σ_photo::Vector{Float64}
    log_σ_Compton::Vector{Float64}
    log_σ_pair::Vector{Float64}
end

"""
    load_xcom(path) -> XCOMTable

Read a NIST XCOM table (8-column format with Coherent / Incoherent /
Photoelectric / Nuclear PP / Electronic PP / Tot. w/ / Tot. wo/).
Drops the K-edge duplicate row by keeping the second (above-edge) entry
when two consecutive rows share the same energy. Pair production is
floored to 1e-30 below threshold so the log interpolation stays defined.
"""
function load_xcom(path::AbstractString)::XCOMTable
    energies = Float64[]
    photos   = Float64[]
    comptons = Float64[]
    pairs    = Float64[]

    for raw in eachline(path)
        line = strip(raw)
        isempty(line) && continue
        cols = split(line)
        length(cols) < 8 && continue
        e1 = tryparse(Float64, cols[1])
        e1 === nothing && continue   # header rows

        E      = e1
        comp   = parse(Float64, cols[3])  # incoherent (Compton)
        pe     = parse(Float64, cols[4])  # photoelectric
        pp_nuc = parse(Float64, cols[5])  # nuclear pair production
        pp_ele = parse(Float64, cols[6])  # electronic pair production
        pp     = pp_nuc + pp_ele

        # K-edge handling: if energy duplicates the previous row, keep
        # the second (above-edge) entry.
        if !isempty(energies) && E == energies[end]
            energies[end] = E
            photos[end]   = pe
            comptons[end] = comp
            pairs[end]    = pp
        else
            push!(energies, E)
            push!(photos,   pe)
            push!(comptons, comp)
            push!(pairs,    pp)
        end
    end
    isempty(energies) && error("load_xcom: no data parsed from $path")

    pp_floor = 1.0e-30
    pairs_safe = [max(p, pp_floor) for p in pairs]

    log_E  = log.(energies)
    log_pe = log.(photos)
    log_co = log.(comptons)
    log_pp = log.(pairs_safe)

    XCOMTable(energies, photos, comptons, pairs,
              log_E, log_pe, log_co, log_pp)
end

# ---------------------------------------------------------------------------
# Log-log interpolation
# ---------------------------------------------------------------------------

@inline function _log_log_interp(log_x::Vector{Float64},
                                 log_y::Vector{Float64},
                                 x_query::Float64)::Float64
    lq = log(x_query)
    idx = searchsortedlast(log_x, lq)
    idx < 1 && error("Energy $(x_query) MeV below table range")
    idx >= length(log_x) && (idx = length(log_x) - 1)
    t = (lq - log_x[idx]) / (log_x[idx+1] - log_x[idx])
    return exp(log_y[idx] + t * (log_y[idx+1] - log_y[idx]))
end

"σ_photo(table, E_MeV) — photoelectric mass attenuation (cm²/g)."
σ_photo(t::XCOMTable, E::Float64)   = _log_log_interp(t.log_E, t.log_σ_photo,   E)

"σ_Compton(table, E_MeV) — Compton mass attenuation (cm²/g)."
σ_Compton(t::XCOMTable, E::Float64) = _log_log_interp(t.log_E, t.log_σ_Compton, E)

"σ_pair(table, E_MeV) — pair-production mass attenuation (cm²/g). Zero below 1.022 MeV."
function σ_pair(t::XCOMTable, E::Float64)
    E < 1.022 && return 0.0
    _log_log_interp(t.log_E, t.log_σ_pair, E)
end

"""
    μ_total_lin(table, E_MeV, ρ) -> Float64

Total linear attenuation coefficient (cm⁻¹). Sums all three components
and multiplies by density. Used to draw photon step lengths.
"""
function μ_total_lin(t::XCOMTable, E::Float64, ρ::Float64)
    (σ_photo(t, E) + σ_Compton(t, E) + σ_pair(t, E)) * ρ
end

# ---------------------------------------------------------------------------
# Klein-Nishina sampling (Compton scatter)
# ---------------------------------------------------------------------------

"""
    sample_klein_nishina(rng, E_MeV) -> (E_scattered_MeV, cos_θ)

Sample the scattered-photon energy and polar angle from the
Klein-Nishina cross section using Kahn's rejection method.
The two values satisfy the kinematic constraint
    1 / E' − 1 / E = (1 − cos θ) / m_e c².
"""
function sample_klein_nishina(rng::AbstractRNG, E_MeV::Float64)
    κ = E_MeV / ME_C2_MEV
    while true
        r1 = rand(rng)
        r2 = rand(rng)
        r3 = rand(rng)

        if r1 ≤ (2κ + 1) / (2κ + 9)
            α = 1.0 + 2κ * r2
            if r3 ≤ 4.0 * (1.0 / α - 1.0 / (α * α))
                E_prime = E_MeV / α
                cos_θ = 1.0 - ME_C2_MEV * (1.0 / E_prime - 1.0 / E_MeV)
                return (E_prime, clamp(cos_θ, -1.0, 1.0))
            end
        else
            α = (1.0 + 2κ) / (1.0 + 2κ * r2)
            cos_θ_cand = 1.0 - (α - 1.0) / κ
            if r3 ≤ 0.5 * (cos_θ_cand * cos_θ_cand + 1.0 / α)
                E_prime = E_MeV / α
                cos_θ = 1.0 - ME_C2_MEV * (1.0 / E_prime - 1.0 / E_MeV)
                return (E_prime, clamp(cos_θ, -1.0, 1.0))
            end
        end
    end
end

# ---------------------------------------------------------------------------
# Direction rotation
# ---------------------------------------------------------------------------

"""
    rotate_direction(dx, dy, dz, cos_θ, φ) -> (dx', dy', dz')

Rotate the unit vector (dx, dy, dz) by polar angle θ (given as cos θ
from the input direction) and azimuth φ around it. Builds an orthogonal
basis perpendicular to the input direction; numerically stable near the
poles. Returns a unit vector.
"""
function rotate_direction(dx::Float64, dy::Float64, dz::Float64,
                          cos_θ::Float64, φ::Float64)
    sin_θ = sqrt(max(0.0, 1.0 - cos_θ * cos_θ))
    cos_φ = cos(φ)
    sin_φ = sin(φ)

    # Build a basis (u, v) ⟂ d. Pick u = ẑ × d unless d ≈ ẑ.
    if abs(dz) < 0.9
        ux = -dy;  uy = dx;  uz = 0.0
    else
        ux = 0.0;  uy = -dz;  uz = dy
    end
    norm_u = sqrt(ux*ux + uy*uy + uz*uz)
    ux /= norm_u;  uy /= norm_u;  uz /= norm_u

    vx = dy * uz - dz * uy
    vy = dz * ux - dx * uz
    vz = dx * uy - dy * ux

    dx_new = sin_θ * cos_φ * ux + sin_θ * sin_φ * vx + cos_θ * dx
    dy_new = sin_θ * cos_φ * uy + sin_θ * sin_φ * vy + cos_θ * dy
    dz_new = sin_θ * cos_φ * uz + sin_θ * sin_φ * vz + cos_θ * dz

    norm = sqrt(dx_new^2 + dy_new^2 + dz_new^2)
    return (dx_new / norm, dy_new / norm, dz_new / norm)
end
