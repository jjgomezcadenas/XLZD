# src3/XLZD3.jl — Module root for the geometry-driven refactor.

module XLZD3

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
include("sampling.jl")
include("mc_params.jl")
include("stack.jl")
include("clusters.jl")
include("histograms.jl")
include("select.jl")
include("classify.jl")
include("tracker.jl")
include("mc.jl")
include("run.jl")

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

# Source sampling
export build_cdf, sample_u
export sample_barrel_entry, sample_endcap_entry, sample_entry
export icv_top_inner_disk, icv_bot_inner_disk
export ICV_TOP_ASPECT, ICV_BOT_ASPECT

# MC params and helper functions
export MCParams, in_fv
export E_tracking_cutoff_MeV, Δz_threshold_cm
export path_to_next_region, companion_visible!, companion_reach_prob

# Run driver
export MCResult, run_mc, run_mc_all

# Cluster type and builder
export Cluster, build_clusters

# Cluster selection predicates
export select_SC, select_ROI, select_FV, select_skin

# Event-outcome classifier
export classify_event, CLASSIFY_EVENT_OUTCOMES, TRACK_STATUSES

# Stack tracker (Phase 2; src3/tracker.jl)
export track_photon_stack, fast_veto

# Stack & cluster diagnostic histograms
export StackHistogramSet, update_stack_histograms!, merge_stack_histograms!
export ClusterHistogramSet, update_cluster_histograms!, merge_cluster_histograms!

# Cut-flow histograms
export CutHistograms, merge_cut_histograms!
export fill_cut1_u!, fill_cut2_first_interaction!, fill_cut3!, fill_cut4!

# Stack tracker (Phase 1 of src3 refactor)
export StackRow, PhotonStack, push_row!
export INT_PHOTO, INT_COMPTON, INT_PAIR, INT_BELOW_THRESH

end # module XLZD3
