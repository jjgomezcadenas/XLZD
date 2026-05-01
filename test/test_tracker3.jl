# test/test_tracker3.jl — Tests for the stack-based tracker.
#
# Step 11b: source sampling + transparent advance + first-interaction
# loop with forced PHOTO. Cross-section sampling, Compton recursion,
# pair production, and early-reject paths land in steps 11c-11f.

using Test
using Random
using Printf

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

const LXE_INTERACTIVE = (:TPC, :Skin, :Inert)

# ---------------------------------------------------------------------------
# 1. Status alphabet
# ---------------------------------------------------------------------------

@testset "1. status is in TRACK_STATUSES" begin
    rng   = MersenneTwister(0xABCD)
    stack = PhotonStack()
    for _ in 1:200
        empty!(stack)
        s = track_photon_stack(rng, det, eff_test, xcom, params, stack)
        @test s isa Symbol
        @test s in TRACK_STATUSES
        @test s in (:completed, :escaped)   # 11b only produces these two
    end
end

# ---------------------------------------------------------------------------
# 2/3. :completed -> 1 row;  :escaped -> 0 rows
# ---------------------------------------------------------------------------

@testset "2/3. row count matches return status" begin
    rng   = MersenneTwister(0xBEEF)
    stack = PhotonStack()
    n_completed = 0
    n_escaped   = 0
    for _ in 1:500
        empty!(stack)
        s = track_photon_stack(rng, det, eff_test, xcom, params, stack)
        if s === :completed
            n_completed += 1
            @test length(stack)  == 1
            @test stack.next_ng  == 2
        elseif s === :escaped
            n_escaped += 1
            @test length(stack)  == 0
            @test stack.next_ng  == 1
        end
    end
    @test n_completed > 0
    @test n_escaped   > 0
end

# ---------------------------------------------------------------------------
# 4. Forced PHOTO + full deposit
# ---------------------------------------------------------------------------

@testset "4. row interaction is one of {PHOTO, COMPTON, PAIR}; edep/epre invariants" begin
    # 11c: cross-section sampling — interaction is no longer forced PHOTO.
    rng   = MersenneTwister(0xCAFE)
    stack = PhotonStack()
    for _ in 1:300
        empty!(stack)
        s = track_photon_stack(rng, det, eff_test, xcom, params, stack)
        if s === :completed
            r = stack.rows[1]
            @test r.interaction in (INT_PHOTO, INT_COMPTON, INT_PAIR)
            @test r.epre ≈ eff_test.E_MeV
            @test 0.0 < r.edep <= r.epre + 1e-9
            if r.interaction === INT_PHOTO
                @test r.edep ≈ r.epre
            end
        end
    end
end

# ---------------------------------------------------------------------------
# 5. parent_region matches eff.region
# ---------------------------------------------------------------------------

@testset "5. row.parent_region == eff.region" begin
    rng   = MersenneTwister(0x1234)
    stack = PhotonStack()
    for _ in 1:300
        empty!(stack)
        s = track_photon_stack(rng, det, eff_test, xcom, params, stack)
        if s === :completed
            @test stack.rows[1].parent_region === eff_test.region
        end
    end
    # Also exercise the endcap source so we catch any source-region mishandling
    eff_cth = by_name["CTH_Bi214"]
    rng2 = MersenneTwister(0x5678)
    for _ in 1:300
        empty!(stack)
        s = track_photon_stack(rng2, det, eff_cth, xcom, params, stack)
        if s === :completed
            @test stack.rows[1].parent_region === eff_cth.region
        end
    end
end

# ---------------------------------------------------------------------------
# 6. Row's region is in the LXe interactor set (never :FC / :Gas / :Outside)
# ---------------------------------------------------------------------------

@testset "6. row.region in {:TPC, :Skin, :Inert}" begin
    rng   = MersenneTwister(0xDEADBEEF)
    stack = PhotonStack()
    for _ in 1:1000
        empty!(stack)
        s = track_photon_stack(rng, det, eff_test, xcom, params, stack)
        if s === :completed
            @test stack.rows[1].region in LXE_INTERACTIVE
        end
    end
end

# ---------------------------------------------------------------------------
# 7. Calibration: most CB_Bi214 photons interact (mfp ~10 cm in 80 cm detector)
# ---------------------------------------------------------------------------

@testset "7. CB_Bi214 :completed fraction in [0.50, 0.70] at N=5000" begin
    # The legacy tracker (track_one_photon!) shows ~0.58 non-escape
    # fraction at the same source/seed; our 11b tracker should match
    # within ~few percent since it samples the same source and uses
    # the same μ_LXe. Once cross-section sampling lands in 11c the
    # non-escape fraction may shift slightly (Compton outgoing γ may
    # escape on its own); we re-tune then.
    rng   = MersenneTwister(0x42)
    stack = PhotonStack()
    n_completed = 0
    N = 5000
    for _ in 1:N
        empty!(stack)
        s = track_photon_stack(rng, det, eff_test, xcom, params, stack)
        s === :completed && (n_completed += 1)
    end
    frac = n_completed / N
    @printf("\n     :completed fraction (CB_Bi214, N=%d): %.3f\n", N, frac)
    @test 0.50 <= frac <= 0.70
end

# ---------------------------------------------------------------------------
# 8. Reproducibility: same seed -> same status sequence
# ---------------------------------------------------------------------------

@testset "8. reproducibility on fixed seed" begin
    function run_seq(seed)
        rng   = MersenneTwister(seed)
        stack = PhotonStack()
        out = Symbol[]
        for _ in 1:200
            empty!(stack)
            push!(out, track_photon_stack(rng, det, eff_test, xcom, params, stack))
        end
        out
    end
    @test run_seq(0x9999) == run_seq(0x9999)
    @test run_seq(0x9999) != run_seq(0xAAAA)
end

# ---------------------------------------------------------------------------
# 9. Stack reuse via empty!
# ---------------------------------------------------------------------------

@testset "9. stack reuse via empty! between events" begin
    rng   = MersenneTwister(0x1111)
    stack = PhotonStack()
    for i in 1:50
        empty!(stack)
        @test length(stack)  == 0
        @test stack.next_ng  == 1
        track_photon_stack(rng, det, eff_test, xcom, params, stack)
        @test length(stack)  in (0, 1)        # 0 if escaped, 1 if completed
    end
end

# ---------------------------------------------------------------------------
# 10. rej_hist accepted; nothing written to it in 11b
# ---------------------------------------------------------------------------

@testset "10. rej_hist accepted; no fills in 11b" begin
    rng   = MersenneTwister(0xCC)
    stack = PhotonStack()
    rej   = RejectionHistograms(r2_max_cm2 = det.R_ICV_inner^2,
                                 z_min_cm   = det.z_LXe_bottom,
                                 z_max_cm   = det.z_gate)
    # Snapshot before
    skin_E_total_before = sum(rej.skin_E_counts)
    fv_E_total_before   = sum(rej.fv_E_counts)
    for _ in 1:300
        empty!(stack)
        track_photon_stack(rng, det, eff_test, xcom, params, stack;
                            rej_hist=rej)
    end
    # 11b does not fill rej_hist (skin/FV early-reject lands in 11f).
    @test sum(rej.skin_E_counts) == skin_E_total_before
    @test sum(rej.fv_E_counts)   == fv_E_total_before
end

# ---------------------------------------------------------------------------
# 11. All three interaction types appear at N=1000; Compton dominates
# ---------------------------------------------------------------------------

@testset "11. PHOTO + COMPTON + PAIR all reachable; Compton dominates" begin
    rng   = MersenneTwister(0x11CC)
    stack = PhotonStack()
    n_photo   = 0
    n_compton = 0
    n_pair    = 0
    for _ in 1:1000
        empty!(stack)
        s = track_photon_stack(rng, det, eff_test, xcom, params, stack)
        if s === :completed
            it = stack.rows[1].interaction
            it === INT_PHOTO   && (n_photo   += 1)
            it === INT_COMPTON && (n_compton += 1)
            it === INT_PAIR    && (n_pair    += 1)
        end
    end
    @printf("\n     interaction-type counts (CB_Bi214, N=1000):\n")
    @printf("       INT_PHOTO   = %d\n", n_photo)
    @printf("       INT_COMPTON = %d\n", n_compton)
    @printf("       INT_PAIR    = %d\n", n_pair)
    @test n_photo   > 0
    @test n_compton > 0
    @test n_pair    > 0
    @test n_compton > n_photo     # Compton dominates at 2.448 MeV in LXe
end

# ---------------------------------------------------------------------------
# 12. Compton edep is strictly between 0 and epre (electron KE, not full)
# ---------------------------------------------------------------------------

@testset "12. Compton: 0 < edep < epre" begin
    rng   = MersenneTwister(0x12CC)
    stack = PhotonStack()
    n_compton_seen = 0
    for _ in 1:1000
        empty!(stack)
        s = track_photon_stack(rng, det, eff_test, xcom, params, stack)
        if s === :completed && stack.rows[1].interaction === INT_COMPTON
            r = stack.rows[1]
            @test 0.0 < r.edep
            @test r.edep < r.epre   # electron KE strictly less than incoming
            n_compton_seen += 1
        end
    end
    @test n_compton_seen > 0    # we should have hit some Compton events
end

# ---------------------------------------------------------------------------
# 13. Pair edep ≈ epre − 2·m_e·c² (deterministic vertex deposit)
# ---------------------------------------------------------------------------

@testset "13. Pair: edep == epre - 2·ME_C2_MEV" begin
    rng   = MersenneTwister(0x13CC)
    stack = PhotonStack()
    n_pair_seen = 0
    for _ in 1:2000
        empty!(stack)
        s = track_photon_stack(rng, det, eff_test, xcom, params, stack)
        if s === :completed && stack.rows[1].interaction === INT_PAIR
            r = stack.rows[1]
            @test r.edep ≈ r.epre - 2.0 * ME_C2_MEV
            n_pair_seen += 1
        end
    end
    @test n_pair_seen > 0
end

# ---------------------------------------------------------------------------
# 14. No PAIR rows from sub-threshold sources (Tl-208 companion, 583 keV)
# ---------------------------------------------------------------------------

@testset "14. sub-threshold source produces no PAIR" begin
    eff_comp = by_name["CB_Tl208c"]   # E_MeV = 0.583, well below 1.022
    @test eff_comp.E_MeV < 2.0 * ME_C2_MEV
    rng   = MersenneTwister(0x14CC)
    stack = PhotonStack()
    n_pair = 0
    for _ in 1:2000
        empty!(stack)
        s = track_photon_stack(rng, det, eff_comp, xcom, params, stack)
        if s === :completed && stack.rows[1].interaction === INT_PAIR
            n_pair += 1
        end
    end
    @test n_pair == 0
end

println("\n  ── test_tracker3.jl: track_photon_stack 11c OK ──\n")
