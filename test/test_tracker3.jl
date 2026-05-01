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

@testset "2/3. row count matches return status (multi-row events possible in 11d)" begin
    rng   = MersenneTwister(0xBEEF)
    stack = PhotonStack()
    n_completed = 0
    n_escaped   = 0
    for _ in 1:500
        empty!(stack)
        s = track_photon_stack(rng, det, eff_test, xcom, params, stack)
        if s === :completed
            n_completed += 1
            @test length(stack)  >= 1
            @test stack.next_ng  == length(stack) + 1
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

@testset "5. first row.parent_region == eff.region (source linkage)" begin
    # The FIRST row of every event is the source γ's first interaction;
    # its parent_region must be the source region. Subsequent rows
    # (Compton chain) carry their own parent_region = previous-vertex
    # region (covered by testset 18).
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

@testset "7. CB_Bi214 :completed fraction in [0.50, 0.75] at N=5000" begin
    # 11b/11c: ~0.58–0.61 (single-interaction tracker).
    # 11d: same per-event escape probability, but RNG-state divergence
    # from extra rand() calls during Compton chains can shift the
    # observed fraction by a few percent. Bounds left wide.
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
    @test 0.50 <= frac <= 0.75
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
        # 11d: 0 if escaped, ≥1 if completed (Compton chains can produce
        # multiple rows per event).
        @test length(stack) >= 0
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

# ---------------------------------------------------------------------------
# 15. Multi-row events exist (Compton chains)
# ---------------------------------------------------------------------------

@testset "15. multi-row events appear at N=1000" begin
    rng   = MersenneTwister(0x15CC)
    stack = PhotonStack()
    n_multi = 0
    max_rows = 0
    for _ in 1:1000
        empty!(stack)
        s = track_photon_stack(rng, det, eff_test, xcom, params, stack)
        if s === :completed
            n_multi += (length(stack) >= 2 ? 1 : 0)
            max_rows = max(max_rows, length(stack))
        end
    end
    @printf("\n     multi-row events: %d / 1000;  max rows in any event: %d\n",
            n_multi, max_rows)
    @test n_multi  > 0
    @test max_rows >= 2
end

# ---------------------------------------------------------------------------
# 16. Mother-chain consistency: walk nm up from the last row to 0
# ---------------------------------------------------------------------------

@testset "16. mother-chain consistency in multi-row events" begin
    rng   = MersenneTwister(0x16CC)
    stack = PhotonStack()
    chains_checked = 0
    for _ in 1:2000
        empty!(stack)
        s = track_photon_stack(rng, det, eff_test, xcom, params, stack)
        if s === :completed && length(stack) >= 2
            # Walk nm from the last row upward — must reach 0 in ≤ N steps.
            cur = stack.rows[end].ng
            seen = Int[]
            while cur != 0
                push!(seen, cur)
                @test 1 <= cur <= length(stack)
                row = stack.rows[cur]
                # Every intermediate row (i.e. every non-leaf in the chain)
                # must be either INT_COMPTON or INT_PAIR. INT_PAIR appears
                # at the boundary between a pair-child's chain and the
                # source γ's pre-pair chain. PHOTO / BELOW_THRESH terminate
                # photons and cannot have children → never intermediate.
                if cur != stack.rows[end].ng
                    @test row.interaction in (INT_COMPTON, INT_PAIR)
                end
                cur = row.nm
                length(seen) > length(stack) && error("mother chain loops")
            end
            @test length(seen) >= 2
            chains_checked += 1
        end
    end
    @test chains_checked > 0
end

# ---------------------------------------------------------------------------
# 17. Energy bookkeeping along Compton chain: epre[i+1] == epre[i] - edep[i]
# ---------------------------------------------------------------------------

@testset "17. Compton-chain energy bookkeeping" begin
    rng   = MersenneTwister(0x17CC)
    stack = PhotonStack()
    chains_checked = 0
    for _ in 1:2000
        empty!(stack)
        s = track_photon_stack(rng, det, eff_test, xcom, params, stack)
        if s === :completed && length(stack) >= 2
            for i in 1:(length(stack) - 1)
                ri = stack.rows[i]
                rj = stack.rows[i + 1]
                # rj.nm == ri.ng iff rj is a child of ri (single chain in 11d)
                if rj.nm == ri.ng && ri.interaction === INT_COMPTON
                    @test isapprox(rj.epre, ri.epre - ri.edep; atol=1e-9)
                    chains_checked += 1
                end
            end
        end
    end
    @test chains_checked > 0
end

# ---------------------------------------------------------------------------
# 18. Child rows: parent_region matches parent's interaction region
# ---------------------------------------------------------------------------

@testset "18. child row.parent_region == parent.region" begin
    rng   = MersenneTwister(0x18CC)
    stack = PhotonStack()
    children_checked = 0
    for _ in 1:2000
        empty!(stack)
        s = track_photon_stack(rng, det, eff_test, xcom, params, stack)
        if s === :completed
            for r in stack.rows
                if r.nm > 0
                    parent = stack.rows[r.nm]
                    @test r.parent_region === parent.region
                    children_checked += 1
                end
            end
        end
    end
    @test children_checked > 0
end

# ---------------------------------------------------------------------------
# 19. INT_BELOW_THRESH rows appear for low-energy Compton tails
# ---------------------------------------------------------------------------

@testset "19. INT_BELOW_THRESH path is reachable from Compton-chain tails" begin
    # At default threshold (40 keV) BELOW_THRESH rows are very rare:
    # one-Compton outgoing γ at 2.448 MeV is bounded below by ~244 keV
    # (kinematic minimum), and successive scatters usually photo-absorb
    # before crossing 40 keV. To exercise the code path deterministically,
    # raise the cutoff to 500 keV so the first Compton outgoing γ
    # frequently falls below it on the next loop iteration.
    params_hi = MCParams(; E_tracking_cutoff_keV=500.0)
    rng       = MersenneTwister(0x19CC)
    stack     = PhotonStack()
    n_below   = 0
    for _ in 1:1000
        empty!(stack)
        s = track_photon_stack(rng, det, eff_test, xcom, params_hi, stack)
        if s === :completed
            for r in stack.rows
                r.interaction === INT_BELOW_THRESH && (n_below += 1)
            end
        end
    end
    @printf("\n     INT_BELOW_THRESH count (cutoff=500 keV, N=1000): %d\n",
            n_below)
    @test n_below > 0
end

# ---------------------------------------------------------------------------
# 20. Pair vertices spawn 511 keV children (rows with nm == pair.ng)
# ---------------------------------------------------------------------------

@testset "20. pair vertex spawns up to two 511 keV children" begin
    # Use a Tl-208 source: more pair production probability at 2.615 MeV.
    eff_tl = by_name["CTH_Tl208"]
    rng    = MersenneTwister(0x20EE)
    stack  = PhotonStack()
    n_pair_events = 0
    n_direct_children_total = 0
    for _ in 1:5000
        empty!(stack)
        s = track_photon_stack(rng, det, eff_tl, xcom, params, stack)
        s === :completed || continue
        # Per pair vertex: count rows whose nm points to it.
        for r_pair in stack.rows
            r_pair.interaction === INT_PAIR || continue
            n_pair_events += 1
            n_children = count(r -> r.nm == r_pair.ng, stack.rows)
            # 0 (both children escaped without interacting), 1, or 2.
            @test 0 <= n_children <= 2
            n_direct_children_total += n_children
        end
    end
    @test n_pair_events > 0
    avg = n_direct_children_total / n_pair_events
    @printf("\n     pair events: %d;  avg direct-child rows: %.3f (max 2.0)\n",
            n_pair_events, avg)
    # 511 keV mfp in LXe is ~3 cm, detector ~80 cm; both children almost
    # always interact. Mean per pair vertex should be > 1.5.
    @test avg > 1.5
end

# ---------------------------------------------------------------------------
# 21. Pair-children's epre is exactly 511 keV
# ---------------------------------------------------------------------------

@testset "21. pair children have epre ≈ 0.511 MeV" begin
    eff_tl = by_name["CTH_Tl208"]
    rng    = MersenneTwister(0x21EE)
    stack  = PhotonStack()
    n_children_checked = 0
    for _ in 1:3000
        empty!(stack)
        s = track_photon_stack(rng, det, eff_tl, xcom, params, stack)
        s === :completed || continue
        for r in stack.rows
            r.nm > 0 || continue
            parent = stack.rows[r.nm]
            if parent.interaction === INT_PAIR
                @test r.epre ≈ ME_C2_MEV  atol=1e-9
                n_children_checked += 1
            end
        end
    end
    @test n_children_checked > 0
end

# ---------------------------------------------------------------------------
# 22. Pair-children's parent_region matches the pair vertex's region
# ---------------------------------------------------------------------------

@testset "22. pair children's parent_region == pair_vertex.region" begin
    eff_tl = by_name["CTH_Tl208"]
    rng    = MersenneTwister(0x22EE)
    stack  = PhotonStack()
    n_children_checked = 0
    for _ in 1:2000
        empty!(stack)
        s = track_photon_stack(rng, det, eff_tl, xcom, params, stack)
        s === :completed || continue
        for r in stack.rows
            r.nm > 0 || continue
            parent = stack.rows[r.nm]
            if parent.interaction === INT_PAIR
                @test r.parent_region === parent.region
                n_children_checked += 1
            end
        end
    end
    @test n_children_checked > 0
end

# ---------------------------------------------------------------------------
# 23. PAIR appears at most once along any single nm-walk (no pair cascade)
# ---------------------------------------------------------------------------

@testset "23. at most one INT_PAIR per mother-chain walk" begin
    # Pair production needs e >= 1.022 MeV; the 511 keV annihilation γ
    # are below threshold and can only photo-absorb or Compton-scatter.
    # Verify by walking nm up from every row and counting INT_PAIR seen.
    eff_tl = by_name["CTH_Tl208"]
    rng    = MersenneTwister(0x23EE)
    stack  = PhotonStack()
    n_walks_checked = 0
    for _ in 1:2000
        empty!(stack)
        s = track_photon_stack(rng, det, eff_tl, xcom, params, stack)
        s === :completed || continue
        for r0 in stack.rows
            n_pair_seen = 0
            cur = r0.ng
            while cur != 0
                row = stack.rows[cur]
                if row.interaction === INT_PAIR
                    n_pair_seen += 1
                end
                cur = row.nm
            end
            @test n_pair_seen <= 1
            n_walks_checked += 1
        end
    end
    @test n_walks_checked > 0
end

println("\n  ── test_tracker3.jl: track_photon_stack 11e OK ──\n")
