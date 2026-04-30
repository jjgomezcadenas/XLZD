# src2/pobjects.jl — Physics-aware geometric objects.
#
# `PCyl` and `PDisk` wrap a `GCyl`/`GDisk` together with a `Material` and
# late-chain specific activities. They expose mass, total activities,
# γ production rates, and a self-shielded angular spectrum at the source's
# own inner surface.
#
# Constants for chain → γ conversion:
#   ²¹⁴Bi 2.448 MeV γ: BR = 1.55% per Bi-214 decay
#   ²⁰⁸Tl 2.615 MeV γ: produced in 35.9% of ²³²Th-late decays (²⁰⁸Tl
#     branching from ²¹²Bi); subsequent γ BR is ~100%, so the per-decay
#     factor is 0.359.

const BR_BI214_GAMMA       = 0.0155
const BR_TL208_FROM_CHAIN  = 0.359
const SEC_PER_YEAR         = 3.1557e7   # Julian year
const E_BI214_MEV          = 2.448
const E_TL208_MEV          = 2.615

# ---------------------------------------------------------------------------
# Types
# ---------------------------------------------------------------------------

"""
    PCyl(geom, material, act_U238_late_mBqkg, act_Th232_late_mBqkg, count=1, name="")

Cylindrical-shell physical object: a `GCyl` plus a `Material` and the
late-chain specific activities (mBq/kg) of ²³⁸U and ²³²Th. `count` is
the multiplicity (number of identical pieces sharing this geometry).
"""
struct PCyl
    geom::GCyl
    material::Material
    act_U238_late_mBqkg::Float64
    act_Th232_late_mBqkg::Float64
    count::Int
    name::String
end
PCyl(geom, material, aU, aTh; count=1, name="") =
    PCyl(geom, material, Float64(aU), Float64(aTh), Int(count), String(name))

"""
    PDisk(geom, material, act_U238_late_mBqkg, act_Th232_late_mBqkg, count=1, name="")

Head/disc physical object: a `GDisk` plus a `Material` and late-chain
specific activities. `count` is the multiplicity.
"""
struct PDisk
    geom::GDisk
    material::Material
    act_U238_late_mBqkg::Float64
    act_Th232_late_mBqkg::Float64
    count::Int
    name::String
end
PDisk(geom, material, aU, aTh; count=1, name="") =
    PDisk(geom, material, Float64(aU), Float64(aTh), Int(count), String(name))

"""
    PSurface(name, attached_to, material_name, total_mass_kg,
             act_U238_late_mBqkg, act_Th232_late_mBqkg)

Surface (zero-thickness) source. The activity is given as
specific-activity × explicit `total_mass_kg` (not derived from a volume).
`attached_to` is a Symbol identifying which Ti wall the surface sits on
(`:ICV_body`, `:ICV_top`, `:ICV_bot`, …); the effective-source registry
uses it to select the right downstream slab list.
"""
struct PSurface
    name::String
    attached_to::Symbol
    material_name::String
    total_mass_kg::Float64
    act_U238_late_mBqkg::Float64
    act_Th232_late_mBqkg::Float64
end

const PObject = Union{PCyl, PDisk, PSurface}

# ---------------------------------------------------------------------------
# Mass and activity
# ---------------------------------------------------------------------------

"Total Ti mass (kg) for this object: shell volume × density × count."
mass(p::PCyl)::Float64 = mass(p.geom, p.material.density) * p.count
mass(p::PDisk)::Float64 = mass(p.geom, p.material.density) * p.count

"Mass (kg) for a surface source: explicit field, no geometry derivation."
mass(p::PSurface)::Float64 = p.total_mass_kg

"Total ²³⁸U-late activity (mBq)."
activity_U238_late(p::PObject)::Float64 = mass(p) * p.act_U238_late_mBqkg

"Total ²³²Th-late activity (mBq)."
activity_Th232_late(p::PObject)::Float64 = mass(p) * p.act_Th232_late_mBqkg

"Bi-214 2.448 MeV γ production rate (γ/yr)."
gamma_rate_Bi214(p::PObject)::Float64 =
    activity_U238_late(p) * 1.0e-3 * BR_BI214_GAMMA * SEC_PER_YEAR

"Tl-208 2.615 MeV γ production rate (γ/yr)."
gamma_rate_Tl208(p::PObject)::Float64 =
    activity_Th232_late(p) * 1.0e-3 * BR_TL208_FROM_CHAIN * SEC_PER_YEAR

# ---------------------------------------------------------------------------
# Self-shielding inward angular spectrum at the source's own inner surface
# ---------------------------------------------------------------------------

"""
    source_slab_thickness(p::PObject) -> Float64 (cm)

The slab thickness used in the self-shielding integral. For both `PCyl`
and `PDisk` this is the wall thickness in the inward direction (radial
for cylinders, normal-to-surface for discs).
"""
source_slab_thickness(p::PCyl)::Float64    = p.geom.wall_thickness
source_slab_thickness(p::PDisk)::Float64   = p.geom.wall_thickness
source_slab_thickness(::PSurface)::Float64 = 0.0

"""
    self_shielded_spectrum(p, gamma_rate_per_yr, E_MeV, u_bins) -> dNdu

Compute the angular spectrum dN/du(u) (γ/yr per unit u = cos θ, with
θ measured from the local inward normal) at the source's own inner
surface, after self-shielding within the source's Ti slab. Birth depth
is uniform in [0, t_src]; only the inward hemisphere contributes
(factor 1/2). Plane-parallel slab approximation:

    dN/du(u) = (R / 2 t_src) · (u / μ_src) · (1 − exp(−μ_src·t_src/u))

where R is the total γ-production rate (γ/yr) and t_src the source
slab thickness. The integral over u ∈ (0, 1] gives the inward-going
γ/yr at the inner surface.
"""
function self_shielded_spectrum(p::Union{PCyl,PDisk}, gamma_rate_per_yr::Real,
                                E_MeV::Real, u_bins::Vector{Float64})
    μ_src = p.material.μ_lin(E_MeV)
    t_src = source_slab_thickness(p)
    R = Float64(gamma_rate_per_yr)
    dNdu = similar(u_bins)
    @inbounds for i in eachindex(u_bins)
        u = u_bins[i]
        if u <= 0
            dNdu[i] = 0.0
        else
            dNdu[i] = (R / (2.0 * t_src)) * (u / μ_src) *
                      (1.0 - exp(-μ_src * t_src / u))
        end
    end
    dNdu
end

"""
    self_shielded_spectrum(p::PSurface, R, E_MeV, u_bins) -> dN/du

Surface source: emission is isotropic into the inward hemisphere with
no self-shielding. Per inward hemisphere `dN/du(u) = R/2` (constant in u);
integrating over u ∈ [0, 1] recovers `R/2` (the inward fraction of the
total γ rate). The `E_MeV` argument is ignored; the surface has no
material-thickness dependence.
"""
function self_shielded_spectrum(p::PSurface, gamma_rate_per_yr::Real,
                                E_MeV::Real, u_bins::Vector{Float64})
    R = Float64(gamma_rate_per_yr)
    dNdu = similar(u_bins)
    val = R / 2.0
    @inbounds for i in eachindex(u_bins)
        dNdu[i] = u_bins[i] > 0 ? val : 0.0
    end
    dNdu
end
