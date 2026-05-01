# src3/transport.jl — Plane-parallel slab transmission for the analytic
# downstream attenuation in the angular spectrum.
#
# A `Slab` is a homogeneous Ti (or other) layer with linear attenuation
# coefficient μ (cm⁻¹) and thickness t (cm). The inward transmission
# at angle θ from the slab normal (u = cos θ) is exp(−μ·t / u).
# A list of slabs adds optical depths: τ = Σᵢ μᵢ·tᵢ.

"""
    Slab(μ, t_cm, label="")

Plane-parallel slab with linear attenuation `μ` (cm⁻¹) and thickness
`t_cm` (cm). `label` is an optional human tag for debugging /
introspection.
"""
struct Slab
    μ::Float64
    t_cm::Float64
    label::String
end
Slab(μ::Real, t_cm::Real, label::AbstractString="") =
    Slab(Float64(μ), Float64(t_cm), String(label))

"""
    optical_depth(slabs::Vector{Slab}) -> Float64

Total normal-incidence optical depth Σ μᵢ·tᵢ (dimensionless).
"""
optical_depth(slabs::Vector{Slab})::Float64 =
    sum(s.μ * s.t_cm for s in slabs; init=0.0)

"""
    transmission_factor(slabs, u_bins) -> Vector{Float64}

Inward-transmission factor exp(−τ / u) at each value of u in `u_bins`,
where τ is the total optical depth of the slab list. With an empty
list (no downstream Ti) returns 1 for all u.
"""
function transmission_factor(slabs::Vector{Slab},
                              u_bins::Vector{Float64})::Vector{Float64}
    τ = optical_depth(slabs)
    if τ == 0.0
        return ones(Float64, length(u_bins))
    end
    [exp(-τ / u) for u in u_bins]
end
