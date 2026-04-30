# src2/XLZD2.jl — Module root for the geometry-driven refactor.

module XLZD2

using Random
using Printf
using DelimitedFiles

include("geometry.jl")
include("cryostat.jl")
include("material.jl")
include("pobjects.jl")
include("sources.jl")
include("transport.jl")
include("effective_sources.jl")
include("lxe_detector.jl")
include("physics.jl")

# Geometry primitives
export GCyl, GDisk

# GCyl-specific
export R_outer, height, area_outer, volume_inner

# GDisk-specific
export depth, z_apex, is_flat

# Shared (multi-method)
export area_inner, volume_shell, mass
export sample_inner_surface, inward_normal, path_through_shell

# Cryostat composite
export CryostatExtra, CryostatSurface, Cryostat, build_cryostat
export total_mass, mass_breakdown, mc_active_extras, mc_active_mass

# Materials
export Material, load_material

# Physics-aware geometry
export PCyl, PDisk, PSurface, PObject
export activity_U238_late, activity_Th232_late
export gamma_rate_Bi214, gamma_rate_Tl208, gamma_rate_Tl208_companion
export source_slab_thickness, self_shielded_spectrum

# Constants
export BR_BI214_GAMMA, BR_TL208_FROM_CHAIN, BR_TL208_COMPANION, SEC_PER_YEAR
export E_BI214_MEV, E_TL208_MEV, E_TL208_COMPANION_MEV
export TI_BB0NU_U238_LATE_MBQKG, TI_BB0NU_TH232_LATE_MBQKG

# Sources
export GammaSource, make_gamma_source
export pobjects_from_cryostat, build_individual_sources
export DEFAULT_U_BINS

# Transport
export Slab, optical_depth, transmission_factor

# Effective sources
export SourceContribution, EffectiveSource
export aggregate_dNdu, build_effective_source, build_effective_sources

# LXe detector
export LXeDetector, build_lxe_detector
export region_at, μ_LXe
export active_volume_cm3, active_mass_kg, skin_volume_cm3, skin_mass_kg

# Photon physics
export XCOMTable, load_xcom
export σ_photo, σ_Compton, σ_pair, μ_total_lin
export sample_klein_nishina, rotate_direction
export ME_C2_MEV

end # module XLZD2
