# test/test_tracker3.jl — Smoke tests for the new stack-based tracker.
#
# Step 11a: only the scaffold exists. `track_photon_stack` is a stub
# that returns `:escaped` and leaves the stack untouched. These tests
# confirm the scaffold is in place; richer physics tests land in
# steps 11b–11f as the implementation grows.

using Test
using Random

isdefined(Main, :XLZD3) || include("../src3/XLZD3.jl")
using .XLZD3

# ---------------------------------------------------------------------------
# Setup: real det / eff / xcom / params from data files
# ---------------------------------------------------------------------------
const lxe_csv      = joinpath(@__DIR__, "..", "data", "lxe_detector.csv")
const lxe_nist     = joinpath(@__DIR__, "..", "data", "nist_lxe.csv")
const ti_path      = joinpath(@__DIR__, "..", "data", "nist_ti.csv")
const geom_csv     = joinpath(@__DIR__, "..", "data", "lz_cryo_geometry.csv")
const extras_csv   = joinpath(@__DIR__, "..", "data", "lz_cryo_extras.csv")
const surfaces_csv = joinpath(@__DIR__, "..", "data", "lz_cryo_surface_sources.csv")
const xcom_path    = joinpath(@__DIR__, "..", "data", "nist.csv")

mat_LXe = load_material("LXe", 2.953, lxe_nist)
mat_Ti  = load_material("Ti",  4.510, ti_path)
det     = build_lxe_detector(lxe_csv, mat_LXe)
cryo    = build_cryostat(geom_csv, extras_csv, surfaces_csv)
indiv   = build_individual_sources(cryo, mat_Ti)
effs    = build_effective_sources(indiv, cryo, mat_Ti)
xcom    = load_xcom(xcom_path)
params  = MCParams()
by_name = Dict(e.name => e for e in effs)
const eff_test = by_name["CB_Bi214"]

# ---------------------------------------------------------------------------
# 1. Scaffold: callable, returns a TRACK_STATUSES symbol, stack untouched
# ---------------------------------------------------------------------------

@testset "11a scaffold: track_photon_stack is callable and well-typed" begin
    rng   = MersenneTwister(0xABCD)
    stack = PhotonStack()

    status = track_photon_stack(rng, det, eff_test, xcom, params, stack)

    @test status isa Symbol
    @test status in TRACK_STATUSES
    @test length(stack)   == 0
    @test stack.next_ng   == 1
end

# ---------------------------------------------------------------------------
# 2. empty! reuse pattern (the production caller will use this)
# ---------------------------------------------------------------------------

@testset "11a scaffold: stack reuse via empty!" begin
    rng   = MersenneTwister(0xBEEF)
    stack = PhotonStack()

    for _ in 1:5
        empty!(stack)
        s = track_photon_stack(rng, det, eff_test, xcom, params, stack)
        @test s in TRACK_STATUSES
        @test length(stack) == 0
        @test stack.next_ng == 1
    end
end

# ---------------------------------------------------------------------------
# 3. rej_hist=nothing path doesn't error (default keyword)
# ---------------------------------------------------------------------------

@testset "11a scaffold: rej_hist=nothing default" begin
    rng   = MersenneTwister(0x1234)
    stack = PhotonStack()
    @test track_photon_stack(rng, det, eff_test, xcom, params, stack;
                              rej_hist=nothing) in TRACK_STATUSES
end

# ---------------------------------------------------------------------------
# 4. rej_hist=RejectionHistograms path also accepted (signature smoke)
# ---------------------------------------------------------------------------

@testset "11a scaffold: rej_hist accepts RejectionHistograms" begin
    rng   = MersenneTwister(0xCAFE)
    stack = PhotonStack()
    rej   = RejectionHistograms(r2_max_cm2 = det.R_ICV_inner^2,
                                 z_min_cm   = det.z_LXe_bottom,
                                 z_max_cm   = det.z_gate)
    @test track_photon_stack(rng, det, eff_test, xcom, params, stack;
                              rej_hist=rej) in TRACK_STATUSES
end

# ---------------------------------------------------------------------------
# 5. The 11a stub specifically returns :escaped (documented behavior)
# ---------------------------------------------------------------------------

@testset "11a stub returns :escaped (delete this test in 11b)" begin
    rng   = MersenneTwister(0x11A)
    stack = PhotonStack()
    @test track_photon_stack(rng, det, eff_test, xcom, params, stack) === :escaped
end

println("\n  ── test_tracker3.jl: track_photon_stack scaffold OK ──\n")
