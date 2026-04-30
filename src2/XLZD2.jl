# src2/XLZD2.jl — Module root for the geometry-driven refactor.

module XLZD2

using Random
using Printf
using DelimitedFiles

include("geometry.jl")
include("cryostat.jl")

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

end # module XLZD2
