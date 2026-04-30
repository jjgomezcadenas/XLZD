# src2/material.jl — Material with energy-dependent linear attenuation.
#
# Loads a NIST XCOM table and exposes μ_lin(E) as a callable closure.
# The NIST file format (whitespace-delimited):
#   - some header lines (skipped: any line that does not begin with a number)
#   - data: either 2 columns (energy, total_mass_attn_w_coherent) or 8 columns
#     (energy, coherent, incoherent, photoelectric, nuclear_pp, electron_pp,
#      tot_w_coherent, tot_wo_coherent). We always read column 1 as energy
#     (MeV) and "Tot. w/ Coherent" (col 2 in the 2-col file, col 7 in the 8-col file).

"""
    Material(name, density_g_cm3, μ_lin)

A material with name (e.g. "Ti", "LXe"), bulk density in g/cm³, and a
callable `μ_lin(E_MeV) → cm⁻¹` that interpolates the NIST XCOM table
in log-log space.
"""
struct Material
    name::String
    density::Float64       # g/cm³
    μ_lin::Function        # E (MeV) → cm⁻¹
end

"""
    _read_nist_table(path) -> (energies, mu_over_rho)

Read a NIST XCOM table. Returns parallel arrays of energies (MeV) and
mass attenuation coefficients with coherent (cm²/g).
"""
function _read_nist_table(path::AbstractString)
    energies = Float64[]
    mu_rho   = Float64[]
    for raw in eachline(path)
        line = strip(raw)
        isempty(line) && continue
        parts = split(line)
        # Skip header rows (first token must parse as Float64)
        e = tryparse(Float64, parts[1])
        e === nothing && continue
        local val::Float64
        if length(parts) == 2
            val_p = tryparse(Float64, parts[2])
        elseif length(parts) >= 7
            val_p = tryparse(Float64, parts[7])
        else
            continue
        end
        val_p === nothing && continue
        push!(energies, e)
        push!(mu_rho,   val_p)
    end
    isempty(energies) && error("No data rows parsed from $path")
    energies, mu_rho
end

"""
    load_material(name, density, nist_path) -> Material

Load a NIST mass-attenuation table and build a `Material` with a log-log
linear μ_lin(E) closure.
"""
function load_material(name::AbstractString, density::Real,
                       nist_path::AbstractString)::Material
    energies, mu_rho = _read_nist_table(nist_path)
    log_E  = log.(energies)
    log_mu = log.(mu_rho)
    E_min, E_max = extrema(energies)

    function μ_lin(E_MeV::Real)
        E_MeV >= E_min ||
            error("μ_lin($name): E=$(E_MeV) MeV below table minimum $(E_min)")
        E_MeV <= E_max ||
            error("μ_lin($name): E=$(E_MeV) MeV above table maximum $(E_max)")
        lE = log(Float64(E_MeV))
        idx = searchsortedlast(log_E, lE)
        idx == length(log_E) && (idx -= 1)
        idx == 0 && (idx = 1)
        t = (lE - log_E[idx]) / (log_E[idx+1] - log_E[idx])
        mu_rho_interp = exp(log_mu[idx] + t * (log_mu[idx+1] - log_mu[idx]))
        return mu_rho_interp * density
    end

    Material(String(name), Float64(density), μ_lin)
end
