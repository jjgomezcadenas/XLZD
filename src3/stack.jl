# src3/stack.jl — Per-photon interaction stack for the new tracker.
#
# A `StackRow` records ONE physics interaction during photon transport.
# Multiple rows assembled in time order form the full event history,
# including pair-production secondaries (each child photon produces its
# own rows, with `nm` pointing to the parent interaction's `ng`).
#
# Field semantics:
#   ng              — global interaction index in the event (1, 2, 3, …)
#   nm              — `ng` of the parent interaction (the one that produced
#                     the photon being tracked here). `nm = 0` for the
#                     original source photon's first interaction.
#   parent_region   — region of the photon's PREVIOUS interaction (or, for a
#                     fresh source photon's first row, the source region
#                     such as :CB / :CTH / :CBH). Causal, NOT spatial:
#                     transparent regions traversed between two interactions
#                     are not recorded here. See design/tracking.md.
#   region          — region where THIS interaction happened
#                     (:TPC, :Skin, :Inert, :Gas, …).
#   interaction     — :PHOTO, :COMPTON, :PAIR, or :BELOW_THRESH.
#   x, y, z         — interaction position (cm).
#   epre            — photon energy ARRIVING at this interaction (MeV).
#   edep            — energy deposited locally by this interaction (MeV).
#                     For Compton: electron kinetic energy = epre - epost.
#                     For Photo:   = epre.
#                     For Pair:    = epre - 2·m_e c²  (kinetic of e⁺ + e⁻).
#                     For BELOW_THRESH: residual photon energy dumped where
#                                       the photon falls below cutoff.
#
# All units: cm, MeV.

struct StackRow
    ng::Int
    nm::Int
    parent_region::Symbol
    region::Symbol
    interaction::Symbol
    x::Float64
    y::Float64
    z::Float64
    epre::Float64
    edep::Float64
end

# Interaction symbol constants. Use these in callers instead of bare
# Symbol literals so a typo becomes a compile-time UndefVarError.
const INT_PHOTO        = :PHOTO
const INT_COMPTON      = :COMPTON
const INT_PAIR         = :PAIR
const INT_BELOW_THRESH = :BELOW_THRESH

"""
    PhotonStack()

Container of `StackRow` records for ONE source-photon event, plus a
running counter `next_ng` used to assign sequential `ng` values, plus a
running `path_length_LXe` accumulator (cm) for the photon's total
distance through any LXe region (`:TPC`, `:Skin`).

Reuse across events by calling `empty!(stack)` between events; this
clears the rows, resets `next_ng = 1`, and resets `path_length_LXe = 0`.
"""
mutable struct PhotonStack
    rows::Vector{StackRow}
    next_ng::Int
    path_length_LXe::Float64
end

PhotonStack() = PhotonStack(StackRow[], 1, 0.0)

Base.length(stack::PhotonStack) = length(stack.rows)

function Base.empty!(stack::PhotonStack)
    empty!(stack.rows)
    stack.next_ng = 1
    stack.path_length_LXe = 0.0
    stack
end

"""
    push_row!(stack; nm, parent_region, region, interaction,
                     x, y, z, epre, edep) -> Int

Append a `StackRow` to `stack`, assigning the next sequential `ng`.
Returns the assigned `ng` so the caller can use it as `nm` for any
child interactions spawned by this row (Compton outgoing γ, pair γγ).
"""
function push_row!(stack::PhotonStack;
                   nm::Integer,
                   parent_region::Symbol,
                   region::Symbol,
                   interaction::Symbol,
                   x::Real, y::Real, z::Real,
                   epre::Real, edep::Real)::Int
    ng = stack.next_ng
    push!(stack.rows, StackRow(ng, Int(nm), parent_region, region, interaction,
                               Float64(x), Float64(y), Float64(z),
                               Float64(epre), Float64(edep)))
    stack.next_ng = ng + 1
    return ng
end
