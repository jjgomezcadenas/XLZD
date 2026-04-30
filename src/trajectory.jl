# trajectory.jl — Per-photon trajectory recording for the 3D viewer

const OUTCOME_CATEGORIES = [:SS_in_ROI, :SS_outside_ROI, :SS_outside_FV,
                             :MS_rejected, :escaped, :outside_bfv,
                             :companion_vetoed, :skin_vetoed]

"""
    InteractionPoint

One interaction along a photon's path.
`type` is one of :entry, :compton, :photo, :pair, :cutoff.
"""
struct InteractionPoint
    x::Float64
    y::Float64
    z::Float64
    E_after::Float64
    type::Symbol
    E_dep::Float64
end

"""
    Trajectory

Complete record of one photon: entry point, interactions, outcome, cluster info.
"""
struct Trajectory
    entry_x::Float64
    entry_y::Float64
    entry_z::Float64
    source_label::String
    dir_x::Float64
    dir_y::Float64
    dir_z::Float64
    interactions::Vector{InteractionPoint}
    outcome::Symbol
    cluster_z::Float64
    cluster_x::Float64
    cluster_y::Float64
    cluster_E::Float64
    cluster_started::Bool
    n_interactions::Int
end

"""
    TrajectoryBuffer

Per-thread storage for trajectory recording.
Fixed-capacity buffer per outcome category.
"""
mutable struct TrajectoryBuffer
    buffers::Dict{Symbol, Vector{Trajectory}}
    max_per_outcome::Int
end

"""
    TrajectoryBuffer(max_per_outcome::Int) -> TrajectoryBuffer

Create an empty buffer with slots for each outcome category.
"""
function TrajectoryBuffer(max_per_outcome::Int)
    bufs = Dict{Symbol, Vector{Trajectory}}()
    for cat in OUTCOME_CATEGORIES
        bufs[cat] = Trajectory[]
    end
    TrajectoryBuffer(bufs, max_per_outcome)
end

"""
    buffer_space_available(buffer::TrajectoryBuffer, outcome::Symbol) -> Bool

Return true if this thread still has space for `outcome` trajectories.
"""
function buffer_space_available(buffer::TrajectoryBuffer, outcome::Symbol)::Bool
    v = get(buffer.buffers, outcome, nothing)
    v === nothing && return false
    return length(v) < buffer.max_per_outcome
end

"""
    commit!(buffer::TrajectoryBuffer, traj::Trajectory)

Store a completed trajectory in the buffer slot for its outcome, if there's room.
"""
function commit!(buffer::TrajectoryBuffer, traj::Trajectory)
    v = get(buffer.buffers, traj.outcome, nothing)
    v === nothing && return
    if length(v) < buffer.max_per_outcome
        push!(v, traj)
    end
    nothing
end

"""
    merge_buffers(buffers::Vector{TrajectoryBuffer}, max_per_outcome::Int) -> Vector{Trajectory}

Merge per-thread buffers into a single list, capped at `max_per_outcome` per outcome.
"""
function merge_buffers(buffers::Vector{TrajectoryBuffer}, max_per_outcome::Int)::Vector{Trajectory}
    merged = Trajectory[]
    for cat in OUTCOME_CATEGORIES
        pool = Trajectory[]
        for buf in buffers
            append!(pool, buf.buffers[cat])
        end
        append!(merged, pool[1:min(length(pool), max_per_outcome)])
    end
    merged
end

"""
    _trajectories_to_list(trajectories::Vector{Trajectory}) -> Vector{Dict}

Convert a vector of Trajectory structs to a list of Dicts for JSON serialization.
"""
function _trajectories_to_list(trajectories::Vector{Trajectory})
    traj_list = []
    for t in trajectories
        ips = [Dict(
            "x" => ip.x, "y" => ip.y, "z" => ip.z,
            "E_after" => ip.E_after, "type" => string(ip.type),
            "E_dep" => ip.E_dep
        ) for ip in t.interactions]

        cluster = t.cluster_started ?
            Dict("x" => t.cluster_x, "y" => t.cluster_y,
                 "z" => t.cluster_z, "E_total" => t.cluster_E) :
            nothing

        push!(traj_list, Dict(
            "entry" => Dict("x" => t.entry_x, "y" => t.entry_y,
                            "z" => t.entry_z, "source" => t.source_label),
            "direction" => Dict("dx" => t.dir_x, "dy" => t.dir_y, "dz" => t.dir_z),
            "outcome" => string(t.outcome),
            "interactions" => ips,
            "cluster" => cluster,
            "n_interactions" => t.n_interactions
        ))
    end
    traj_list
end

"""
    to_json(trajectories, trajectories_fv, geom, params) -> String

Serialize scene metadata + trajectories (all + FV-only) into JSON for `viewer.html`.
"""
function to_json(trajectories::Vector{Trajectory}, trajectories_fv::Vector{Trajectory},
                 geom::Geometry, params::Params)::String
    data = Dict(
        "geometry" => Dict(
            "R_lxe" => geom.R_lxe, "L_lxe" => geom.L_lxe,
            "fv_z_min" => params.fv_z_min_cm, "fv_z_max" => params.fv_z_max_cm,
            "fv_r" => sqrt(params.fv_r2_max_cm2)
        ),
        "trajectories" => _trajectories_to_list(trajectories),
        "trajectories_fv" => _trajectories_to_list(trajectories_fv)
    )

    JSON3.write(data)
end
