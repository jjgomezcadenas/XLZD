# physics.jl — Cross sections, Klein-Nishina, direction rotation

const ME_C2_MEV = 0.510999  # electron rest mass in MeV

"""
    XCOMTable

NIST XCOM cross-section data for xenon, loaded once at startup.
Stores log-log interpolation grids for photoelectric and Compton cross sections.
"""
struct XCOMTable
    energy_MeV::Vector{Float64}
    σ_photo::Vector{Float64}       # cm²/g
    σ_Compton::Vector{Float64}     # cm²/g
    σ_pair::Vector{Float64}        # cm²/g (nuclear + electronic pair production)
    log_E::Vector{Float64}
    log_σ_photo::Vector{Float64}
    log_σ_Compton::Vector{Float64}
    log_σ_pair::Vector{Float64}
end

"""
    load_xcom(path::String) -> XCOMTable

Read `data/nist.csv`, drop the below-K-edge duplicate row at 34.56 keV,
and set up log-log interpolation grids.
"""
function load_xcom(path::String)::XCOMTable
    # Parse whitespace-delimited data, skipping 2 header lines
    lines = readlines(path)
    energies  = Float64[]
    photos    = Float64[]
    comptons  = Float64[]
    pairs     = Float64[]

    for line in lines[3:end]
        stripped = strip(line)
        isempty(stripped) && continue
        cols = split(stripped)
        length(cols) < 8 && continue
        E    = parse(Float64, cols[1])
        comp = parse(Float64, cols[3])  # incoherent scatter = Compton
        pe   = parse(Float64, cols[4])  # photoelectric absorption
        pp_nuc = parse(Float64, cols[5])  # nuclear pair production
        pp_ele = parse(Float64, cols[6])  # electronic pair production
        pp = pp_nuc + pp_ele              # total pair production

        # K-edge handling: if energy duplicates previous, keep the later (above-edge) row
        if !isempty(energies) && E == energies[end]
            energies[end] = E
            photos[end]   = pe
            comptons[end]  = comp
            pairs[end]     = pp
        else
            push!(energies, E)
            push!(photos, pe)
            push!(comptons, comp)
            push!(pairs, pp)
        end
    end

    # Pair production is zero below threshold (~1.022 MeV). To avoid log(0)
    # in log-log interpolation, clamp to a tiny floor value.
    pp_floor = 1.0e-30
    pairs_safe = [max(p, pp_floor) for p in pairs]

    log_E  = log.(energies)
    log_pe = log.(photos)
    log_co = log.(comptons)
    log_pp = log.(pairs_safe)

    XCOMTable(energies, photos, comptons, pairs, log_E, log_pe, log_co, log_pp)
end

"""
    _log_log_interp(log_x, log_y, x_query) -> Float64

Log-log linear interpolation. Assumes `log_x` is sorted ascending.
"""
function _log_log_interp(log_x::Vector{Float64}, log_y::Vector{Float64},
                         x_query::Float64)::Float64
    lq = log(x_query)
    # Find bracketing interval via binary search
    idx = searchsortedlast(log_x, lq)
    idx < 1 && error("Energy $(x_query) MeV below table range")
    idx >= length(log_x) && (idx = length(log_x) - 1)
    t = (lq - log_x[idx]) / (log_x[idx+1] - log_x[idx])
    return exp(log_y[idx] + t * (log_y[idx+1] - log_y[idx]))
end

"""
    σ_photo(table::XCOMTable, E_MeV::Float64) -> Float64

Log-log interpolated photoelectric cross section (cm²/g).
"""
function σ_photo(table::XCOMTable, E_MeV::Float64)::Float64
    _log_log_interp(table.log_E, table.log_σ_photo, E_MeV)
end

"""
    σ_Compton(table::XCOMTable, E_MeV::Float64) -> Float64

Log-log interpolated Compton (incoherent) cross section (cm²/g).
"""
function σ_Compton(table::XCOMTable, E_MeV::Float64)::Float64
    _log_log_interp(table.log_E, table.log_σ_Compton, E_MeV)
end

"""
    σ_pair(table::XCOMTable, E_MeV::Float64) -> Float64

Log-log interpolated pair production cross section (cm²/g).
Sum of nuclear and electronic pair production. Zero below ~1.022 MeV.
"""
function σ_pair(table::XCOMTable, E_MeV::Float64)::Float64
    # Below pair production threshold, return zero (don't interpolate the floor)
    E_MeV < 1.022 && return 0.0
    _log_log_interp(table.log_E, table.log_σ_pair, E_MeV)
end

"""
    μ_total_lin(table::XCOMTable, E_MeV::Float64, ρ_LXe::Float64) -> Float64

Total linear attenuation coefficient (cm⁻¹) = (σ_photo + σ_Compton + σ_pair) × ρ_LXe.
Includes pair production so that the mean free path and interaction branching
are physically correct. When pair production is sampled, the event is
immediately classified as MS-rejected (the two 511 keV annihilation gammas
travel ~6 cm in LXe and virtually always produce spatially separated deposits).
"""
function μ_total_lin(table::XCOMTable, E_MeV::Float64, ρ_LXe::Float64)::Float64
    (σ_photo(table, E_MeV) + σ_Compton(table, E_MeV) + σ_pair(table, E_MeV)) * ρ_LXe
end

"""
    sample_klein_nishina(rng, E_MeV::Float64) -> (E_scattered_MeV, cos_theta)

Sample scattered photon energy from Klein-Nishina via Kahn's rejection method.
Returns scattered energy and cos(θ) from the Compton kinematic relation.
"""
function sample_klein_nishina(rng::AbstractRNG, E_MeV::Float64)
    κ = E_MeV / ME_C2_MEV
    # Kahn's algorithm
    while true
        r1 = rand(rng)
        r2 = rand(rng)
        r3 = rand(rng)

        if r1 ≤ (2κ + 1) / (2κ + 9)
            # Branch 1
            α = 1.0 + 2κ * r2
            if r3 ≤ 4.0 * (1.0 / α - 1.0 / (α * α))
                E_prime = E_MeV / α
                cos_θ = 1.0 - ME_C2_MEV * (1.0 / E_prime - 1.0 / E_MeV)
                return (E_prime, clamp(cos_θ, -1.0, 1.0))
            end
        else
            # Branch 2
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

"""
    rotate_direction(dx, dy, dz, cos_theta, phi) -> (dx', dy', dz')

Rotate a unit vector by polar angle θ (given as cos θ) and azimuth φ.
Builds an orthogonal basis perpendicular to the input direction.
Numerically stable for directions near the poles.
"""
function rotate_direction(dx::Float64, dy::Float64, dz::Float64,
                          cos_theta::Float64, phi::Float64)
    sin_theta = sqrt(max(0.0, 1.0 - cos_theta * cos_theta))
    cos_phi = cos(phi)
    sin_phi = sin(phi)

    # Build orthogonal basis (u, v) perpendicular to d
    if abs(dz) < 0.9
        # u = normalize(ẑ × d)
        ux = -dy;  uy = dx;  uz = 0.0
    else
        # u = normalize(x̂ × d)
        ux = 0.0;  uy = -dz;  uz = dy
    end
    norm_u = sqrt(ux * ux + uy * uy + uz * uz)
    ux /= norm_u;  uy /= norm_u;  uz /= norm_u

    # v = d × u
    vx = dy * uz - dz * uy
    vy = dz * ux - dx * uz
    vz = dx * uy - dy * ux

    # New direction
    dx_new = sin_theta * cos_phi * ux + sin_theta * sin_phi * vx + cos_theta * dx
    dy_new = sin_theta * cos_phi * uy + sin_theta * sin_phi * vy + cos_theta * dy
    dz_new = sin_theta * cos_phi * uz + sin_theta * sin_phi * vz + cos_theta * dz

    # Renormalize for safety
    norm = sqrt(dx_new^2 + dy_new^2 + dz_new^2)
    return (dx_new / norm, dy_new / norm, dz_new / norm)
end
