module XLZD

using Random
using Base.Threads
using HDF5
using JSON3
using Plots
using StatsBase
using DelimitedFiles
using Printf

include("geometry.jl")
include("physics.jl")
include("trajectory.jl")
include("mc.jl")
include("output.jl")

export Params, Geometry, BFV, SourceConfig, XCOMTable, Result, copy_params
export build_geometry, build_bfv, in_bfv, in_fv, in_skin, in_active_lxe
export sample_entry_point, sample_entry_direction, sample_u
export sample_companion_energy, track_companion_gamma
export path_to_cylinder_exit
export load_xcom, σ_photo, σ_Compton, σ_pair, μ_total_lin
export sample_klein_nishina, rotate_direction
export run_mc, write_outputs

end # module XLZD
