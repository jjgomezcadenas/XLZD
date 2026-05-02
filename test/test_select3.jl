# test/test_select3.jl — Pure-predicate tests for select_SC + select_ROI.
#
# Clusters are constructed by hand with explicit (ec, es) values so the
# tests are fully deterministic — no RNG, no smearing happens here.
# The smearing distribution is exercised in test_build_clusters3.jl.

using Test
using Printf

isdefined(Main, :XLZD3) || include("../src3/XLZD3.jl")
using .XLZD3

const params = MCParams()                   # Q_ββ = 2458 keV, halfwidth = 17.2 keV
const Q_MeV  = params.Q_betabeta_keV / 1000.0
const HW_keV = params.ROI_halfwidth_keV
const HW_MeV = HW_keV / 1000.0

# Helper: hand-build a cluster with given (ec, es).
mk(ec, es; x=0.0, y=0.0, z=50.0) = Cluster(x, y, z, ec, es)

# ---------------------------------------------------------------------------
# select_SC — visible-cluster-count predicate
# (sub-visible clusters with ec ≤ E_visible_keV/1000 are ignored — this is
#  what makes the SS/MS check geometry-invariant. See test_select3 cases
#  23-27 below for the visibility-filter behaviour.)
# ---------------------------------------------------------------------------

const E_VIS_MEV_SC = params.E_visible_keV / 1000.0

@testset "1. select_SC: empty list → false" begin
    @test select_SC(Cluster[], params) === false
end

@testset "2. select_SC: one visible cluster → true" begin
    @test select_SC([mk(Q_MeV, Q_MeV)], params) === true
end

@testset "3. select_SC: two visible clusters → false" begin
    @test select_SC([mk(1.0, 1.0), mk(0.5, 0.5)], params) === false
end

@testset "4. select_SC: three visible clusters → false" begin
    @test select_SC([mk(0.3, 0.3), mk(0.7, 0.7), mk(1.2, 1.2)], params) === false
end

@testset "4a. select_SC: 1 visible + N sub-visible → true (regression for the SS/MS bug)" begin
    # The bug we're guarding against: build_clusters keeps every :TPC row
    # with edep > 0, so a vertical photon descending through LXe creates
    # extra sub-visible clusters at z-spaced positions. Counting those in
    # select_SC penalised vertical trajectories (e.g. CTH) far more than
    # horizontal ones (CB), producing a spurious SS_in_ROI suppression.
    sub = E_VIS_MEV_SC * 0.5         # well below visible threshold
    cs = [mk(Q_MeV, Q_MeV; z=50.0),
          mk(sub,    sub;    z=20.0),
          mk(sub,    sub;    z=80.0)]
    @test select_SC(cs, params) === true
end

@testset "4b. select_SC: 0 visible + many sub-visible → false" begin
    sub = E_VIS_MEV_SC * 0.3
    cs = [mk(sub, sub; z=20.0), mk(sub, sub; z=50.0), mk(sub, sub; z=80.0)]
    @test select_SC(cs, params) === false
end

@testset "4c. select_SC: 2 visible + N sub-visible → false" begin
    sub = E_VIS_MEV_SC * 0.5
    cs = [mk(1.0, 1.0; z=30.0),
          mk(sub, sub; z=50.0),
          mk(0.7, 0.7; z=70.0)]
    @test select_SC(cs, params) === false
end

@testset "4d. select_SC: visibility threshold is strict (ec == thr → not visible)" begin
    cs = [mk(E_VIS_MEV_SC, E_VIS_MEV_SC; z=50.0)]    # ec exactly at threshold
    @test select_SC(cs, params) === false
end

# ---------------------------------------------------------------------------
# select_ROI — deterministic Bool from cluster.es
# ---------------------------------------------------------------------------

@testset "5. select_ROI: es == Q_ββ → true" begin
    c = mk(Q_MeV, Q_MeV)
    @test select_ROI(c, params) === true
end

@testset "6. select_ROI: es just inside +halfwidth → true" begin
    ε_keV = 0.01
    c = mk(Q_MeV, Q_MeV + (HW_keV - ε_keV) / 1000.0)
    @test select_ROI(c, params) === true
end

@testset "7. select_ROI: es just outside +halfwidth → false" begin
    ε_keV = 0.01
    c = mk(Q_MeV, Q_MeV + (HW_keV + ε_keV) / 1000.0)
    @test select_ROI(c, params) === false
end

@testset "8. select_ROI: es just outside -halfwidth → false" begin
    ε_keV = 0.01
    c = mk(Q_MeV, Q_MeV - (HW_keV + ε_keV) / 1000.0)
    @test select_ROI(c, params) === false
end

@testset "9. select_ROI: es far below Q (1.0 MeV) → false" begin
    @test select_ROI(mk(1.0, 1.0), params) === false
end

@testset "10. select_ROI: ec is irrelevant; only es matters" begin
    # ec deliberately at Q to confirm the function does NOT read ec.
    c_in  = mk(Q_MeV, Q_MeV)                                # es = Q  → in
    c_out = mk(Q_MeV, Q_MeV + 5 * HW_MeV)                  # es far  → out
    @test select_ROI(c_in,  params) === true
    @test select_ROI(c_out, params) === false
end

# ---------------------------------------------------------------------------
# select_FV — slow per-cluster FV check
# ---------------------------------------------------------------------------

# FV defaults (mirroring test_classify_event3): z ∈ [26, 96], r² ≤ 39²=1521.
const Z_MID_FV  = 0.5 * (params.fv_z_min_cm + params.fv_z_max_cm)
const Z_OUT_FV  = params.fv_z_max_cm + 5.0
const R_OUT_FV  = sqrt(params.fv_r2_max_cm2) + 5.0
const E_VIS_MEV = params.E_visible_keV / 1000.0

@testset "11. select_FV: empty → true" begin
    @test select_FV(Cluster[], params) === true
end

@testset "12. select_FV: one cluster in FV (ec=Q) → true" begin
    cs = [Cluster(0.0, 0.0, Z_MID_FV, Q_MeV, Q_MeV)]
    @test select_FV(cs, params) === true
end

@testset "13. select_FV: one cluster outside FV with ec > visible → false" begin
    cs = [Cluster(0.0, 0.0, Z_OUT_FV, Q_MeV, Q_MeV)]
    @test select_FV(cs, params) === false
end

@testset "14. select_FV: sub-visible cluster outside FV → true (ignored)" begin
    # ec just below E_visible_keV/1000: should NOT trigger reject.
    sub_visible = E_VIS_MEV * 0.5
    cs = [Cluster(0.0, 0.0, Z_OUT_FV, sub_visible, sub_visible)]
    @test select_FV(cs, params) === true
end

@testset "15. select_FV: two clusters both in FV → true" begin
    cs = [Cluster(0.0, 0.0, 30.0, Q_MeV, Q_MeV),
          Cluster(0.0, 0.0, 60.0, 0.5,    0.5)]
    @test select_FV(cs, params) === true
end

@testset "16. select_FV: two clusters, one outside FV (above visible) → false" begin
    cs = [Cluster(0.0, 0.0, 30.0,    Q_MeV, Q_MeV),
          Cluster(R_OUT_FV, 0.0, 60.0, 0.5, 0.5)]
    @test select_FV(cs, params) === false
end

# ---------------------------------------------------------------------------
# select_skin — slow cumulative skin-energy check
# ---------------------------------------------------------------------------

const E_SKIN_VETO_KEV = params.E_skin_veto_keV
const E_SKIN_VETO_MEV = E_SKIN_VETO_KEV / 1000.0

# Helper: build a stack with given (region, edep) pairs.
function _stack_from(rows::Vector{Tuple{Symbol,Float64}})
    s = PhotonStack()
    for (reg, ed) in rows
        push_row!(s; nm=0, parent_region=:CB, region=reg,
                      interaction=INT_PHOTO,
                      x=0.0, y=0.0, z=50.0, epre=ed, edep=ed)
    end
    s
end

@testset "17. select_skin: empty stack → true" begin
    @test select_skin(PhotonStack(), params) === true
end

@testset "18. select_skin: only :TPC rows → true" begin
    s = _stack_from([(:TPC, 1.0), (:TPC, 0.5)])
    @test select_skin(s, params) === true
end

@testset "19. select_skin: one :Skin row below threshold → true" begin
    s = _stack_from([(:Skin, 0.5 * E_SKIN_VETO_MEV)])
    @test select_skin(s, params) === true
end

@testset "20. select_skin: one :Skin row above threshold → false" begin
    s = _stack_from([(:Skin, 1.5 * E_SKIN_VETO_MEV)])
    @test select_skin(s, params) === false
end

@testset "21. select_skin: cumulative > threshold → false" begin
    # Two skin rows each below threshold but summing above.
    half = 0.6 * E_SKIN_VETO_MEV
    s = _stack_from([(:Skin, half), (:Skin, half)])
    @test select_skin(s, params) === false
end

@testset "22. select_skin: mixed regions, skin sum below → true" begin
    s = _stack_from([(:TPC, 1.0),
                     (:Skin, 0.3 * E_SKIN_VETO_MEV),
                     (:Inert, 0.5),
                     (:Skin, 0.4 * E_SKIN_VETO_MEV)])
    @test select_skin(s, params) === true
end

println("\n  ── test_select3.jl: select_SC + select_ROI + select_FV + select_skin OK ──\n")
