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
export CryostatExtra, Cryostat, build_cryostat
export total_mass, mass_breakdown, mc_active_extras, mc_active_mass

# Materials
export Material, load_material

# Physics-aware geometry
export PCyl, PDisk, PObject
export activity_U238_late, activity_Th232_late
export gamma_rate_Bi214, gamma_rate_Tl208
export source_slab_thickness, self_shielded_spectrum

# Constants
export BR_BI214_GAMMA, BR_TL208_FROM_CHAIN, SEC_PER_YEAR
export E_BI214_MEV, E_TL208_MEV
export TI_BB0NU_U238_LATE_MBQKG, TI_BB0NU_TH232_LATE_MBQKG

# Sources
export GammaSource, make_gamma_source
export pobjects_from_cryostat, build_individual_sources
export DEFAULT_U_BINS

end # module XLZD2
