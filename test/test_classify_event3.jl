# test/test_classify_event3.jl — Pure-decision tests for classify_event.
#
# Hand-built (status, clusters) inputs, no MC, no tracking.
# Each testset asserts the resulting outcome symbol.

using Test

isdefined(Main, :XLZD3) || include("../src3/XLZD3.jl")
using .XLZD3

const params = MCParams()
const Q_MeV  = params.Q_betabeta_keV / 1000.0
const HW_keV = params.ROI_halfwidth_keV
const HW_MeV = HW_keV / 1000.0

# Most legacy testsets don't depend on stack contents — only the slow-skin
# path does. Keep one fresh-empty stack handy for the legacy calls.
empty_stack() = PhotonStack()

# FV-bracket helpers from defaults: z ∈ [26, 96] cm, r² ≤ 39² = 1521 cm²
const Z_MID  = 0.5 * (params.fv_z_min_cm + params.fv_z_max_cm)   # well inside FV
const Z_LOW  = params.fv_z_min_cm - 1.0                          # below FV
const Z_HIGH = params.fv_z_max_cm + 1.0                          # above FV
const R_OUT  = sqrt(params.fv_r2_max_cm2) + 1.0                  # outside FV by r

# Cluster ctor: ec=Q to make :SS_in_ROI vs :SS_outside_ROI selectable via es.
mk(; x=0.0, y=0.0, z=Z_MID, ec=Q_MeV, es=Q_MeV) = Cluster(x, y, z, ec, es)

# ---------------------------------------------------------------------------
# 1. status = :vetoed_skin → :skin_vetoed (clusters ignored)
# ---------------------------------------------------------------------------

@testset "1. :vetoed_skin → :skin_vetoed (clusters ignored)" begin
    cs_in_fv  = [mk()]
    cs_empty  = Cluster[]
    cs_two    = [mk(z=30.0), mk(z=60.0)]
    @test classify_event(:vetoed_skin, empty_stack(), cs_in_fv, params) === :skin_vetoed
    @test classify_event(:vetoed_skin, empty_stack(), cs_empty, params) === :skin_vetoed
    @test classify_event(:vetoed_skin, empty_stack(), cs_two,   params) === :skin_vetoed
end

# ---------------------------------------------------------------------------
# 2. status = :rejected_fv → :outside_FV (clusters ignored)
# ---------------------------------------------------------------------------

@testset "2. :rejected_fv → :outside_FV (clusters ignored)" begin
    @test classify_event(:rejected_fv, empty_stack(), [mk()],     params) === :outside_FV
    @test classify_event(:rejected_fv, empty_stack(), Cluster[],  params) === :outside_FV
end

# ---------------------------------------------------------------------------
# 3. status = :escaped → :escaped (clusters ignored)
# ---------------------------------------------------------------------------

@testset "3. :escaped → :escaped" begin
    @test classify_event(:escaped, empty_stack(), Cluster[], params) === :escaped
    @test classify_event(:escaped, empty_stack(), [mk()],    params) === :escaped
end

# ---------------------------------------------------------------------------
# 4. status = :completed, empty clusters → :escaped
# ---------------------------------------------------------------------------

@testset "4. :completed + empty clusters → :escaped" begin
    @test classify_event(:completed, empty_stack(), Cluster[], params) === :escaped
end

# ---------------------------------------------------------------------------
# 5. status = :completed, ≥ 2 clusters → :MS_rejected
# ---------------------------------------------------------------------------

@testset "5. :completed + 2 clusters → :MS_rejected" begin
    @test classify_event(:completed, empty_stack(), [mk(z=40.0), mk(z=60.0)], params) === :MS_rejected
    @test classify_event(:completed, empty_stack(), [mk(z=40.0), mk(z=60.0), mk(z=80.0)], params) === :MS_rejected
end

# ---------------------------------------------------------------------------
# 6. SS, cluster z below FV → :outside_FV
# ---------------------------------------------------------------------------

@testset "6. SS cluster z < fv_z_min → :outside_FV" begin
    @test classify_event(:completed, empty_stack(), [mk(z=Z_LOW)], params) === :outside_FV
end

# ---------------------------------------------------------------------------
# 7. SS, cluster z above FV → :outside_FV
# ---------------------------------------------------------------------------

@testset "7. SS cluster z > fv_z_max → :outside_FV" begin
    @test classify_event(:completed, empty_stack(), [mk(z=Z_HIGH)], params) === :outside_FV
end

# ---------------------------------------------------------------------------
# 8. SS, cluster outside r² → :outside_FV
# ---------------------------------------------------------------------------

@testset "8. SS cluster r² > fv_r2_max → :outside_FV" begin
    @test classify_event(:completed, empty_stack(), [mk(x=R_OUT, y=0.0, z=Z_MID)], params) === :outside_FV
end

# ---------------------------------------------------------------------------
# 9. SS in FV with es at Q → :SS_in_ROI
# ---------------------------------------------------------------------------

@testset "9. SS in FV, es == Q → :SS_in_ROI" begin
    @test classify_event(:completed, empty_stack(), [mk(es=Q_MeV)], params) === :SS_in_ROI
end

# ---------------------------------------------------------------------------
# 10. SS in FV with es far from Q → :SS_outside_ROI
# ---------------------------------------------------------------------------

@testset "10. SS in FV, es far from Q → :SS_outside_ROI" begin
    @test classify_event(:completed, empty_stack(), [mk(es=Q_MeV + 5*HW_MeV)], params) === :SS_outside_ROI
    @test classify_event(:completed, empty_stack(), [mk(es=Q_MeV - 5*HW_MeV)], params) === :SS_outside_ROI
    @test classify_event(:completed, empty_stack(), [mk(es=1.0)],              params) === :SS_outside_ROI
end

# ---------------------------------------------------------------------------
# 11. Coverage: every outcome (except :companion_vetoed) is reachable
# ---------------------------------------------------------------------------

@testset "11. coverage of CLASSIFY_EVENT_OUTCOMES (sans :companion_vetoed)" begin
    reached = Set{Symbol}()
    push!(reached, classify_event(:vetoed_skin, empty_stack(), Cluster[], params))           # :skin_vetoed
    push!(reached, classify_event(:rejected_fv, empty_stack(),   Cluster[], params))           # :outside_FV
    push!(reached, classify_event(:completed, empty_stack(),           Cluster[], params))           # :escaped
    push!(reached, classify_event(:completed, empty_stack(), [mk(z=40.0), mk(z=60.0)], params))     # :MS_rejected
    push!(reached, classify_event(:completed, empty_stack(), [mk(es=Q_MeV)], params))                # :SS_in_ROI
    push!(reached, classify_event(:completed, empty_stack(), [mk(es=1.0)],   params))                # :SS_outside_ROI
    expected = setdiff(Set(CLASSIFY_EVENT_OUTCOMES), Set([:companion_vetoed]))
    # All six output symbols are reachable from the inputs above.
    @test reached == expected
end

# ---------------------------------------------------------------------------
# 12. Edge case: vetoed_skin takes priority over clusters
# ---------------------------------------------------------------------------

@testset "12. vetoed_skin priority over even SS-in-ROI clusters" begin
    # Even a perfect SS-in-ROI cluster cannot rescue a skin-vetoed event.
    @test classify_event(:vetoed_skin, empty_stack(), [mk(es=Q_MeV)], params) === :skin_vetoed
end

# ---------------------------------------------------------------------------
# Slow-check fixtures: build stacks with controlled :Skin contributions.
# ---------------------------------------------------------------------------

const E_SKIN_VETO_MEV = params.E_skin_veto_keV / 1000.0

function _stack_with_skin_total(skin_E_MeV::Float64)
    s = PhotonStack()
    push_row!(s; nm=0, parent_region=:CB, region=:Skin,
                  interaction=INT_PHOTO,
                  x=0.0, y=0.0, z=50.0,
                  epre=skin_E_MeV, edep=skin_E_MeV)
    s
end

# ---------------------------------------------------------------------------
# 13. Slow skin: stack skin sum > threshold overrides clusters
# ---------------------------------------------------------------------------

@testset "13. slow skin: cumulative > threshold → :skin_vetoed" begin
    s = _stack_with_skin_total(1.5 * E_SKIN_VETO_MEV)
    cs = [mk(es=Q_MeV)]   # would otherwise be :SS_in_ROI
    @test classify_event(:completed, s, cs, params) === :skin_vetoed
end

# ---------------------------------------------------------------------------
# 14. Slow FV on SS: cluster outside FV (above visible) → :outside_FV
# ---------------------------------------------------------------------------

@testset "14. slow FV on SS: cluster outside FV → :outside_FV" begin
    cs = [mk(z=Z_HIGH, ec=Q_MeV, es=Q_MeV)]
    @test classify_event(:completed, empty_stack(), cs, params) === :outside_FV
end

# ---------------------------------------------------------------------------
# 15. Slow FV on MS: one cluster outside FV → :outside_FV (overrides MS)
# ---------------------------------------------------------------------------

@testset "15. slow FV on MS: any cluster outside FV → :outside_FV" begin
    cs = [mk(z=Z_MID,  ec=Q_MeV, es=Q_MeV),
          mk(z=Z_HIGH, ec=0.5,    es=0.5)]
    @test classify_event(:completed, empty_stack(), cs, params) === :outside_FV
end

# ---------------------------------------------------------------------------
# 16. SS cluster outside FV with sub-visible ec → :SS_outside_ROI
#     (documents the new spec: only ec > E_visible_keV triggers FV reject)
# ---------------------------------------------------------------------------

@testset "16. lone sub-visible cluster anywhere → :escaped" begin
    sub_visible = params.E_visible_keV / 1000.0 * 0.5
    cs = [mk(z=Z_HIGH, ec=sub_visible, es=sub_visible)]
    # Cluster is sub-visible-threshold → select_FV ignores it AND it doesn't
    # count towards select_SC. With 0 visible clusters there is no signal
    # to classify → :escaped.
    @test classify_event(:completed, empty_stack(), cs, params) === :escaped
    # Same outcome regardless of position (sub-visible everywhere is invisible).
    cs_in = [mk(z=Z_MID, ec=sub_visible, es=sub_visible)]
    @test classify_event(:completed, empty_stack(), cs_in, params) === :escaped
end

# ---------------------------------------------------------------------------
# 17. SS bug regression: 1 visible cluster in FV/ROI + N sub-visible
#     z-spread clusters → must be classified :SS_in_ROI, NOT :MS_rejected.
#
#     This is the bug that produced 0 SS_in_ROI events from CTH_Tl208 at
#     1e9 samples while CB_Tl208 produced ~92. Vertical CTH trajectories
#     create sub-visible clusters at z-spaced positions; the old
#     length-based select_SC counted them and rejected the event as MS.
# ---------------------------------------------------------------------------

@testset "17. SS in FV/ROI + N sub-visible z-spread → :SS_in_ROI (bug regression)" begin
    sub = params.E_visible_keV / 1000.0 * 0.5     # sub-visible
    cs = [mk(z=Z_MID, ec=Q_MeV, es=Q_MeV),         # visible, in FV, in ROI
          mk(z=Z_LOW,  ec=sub,    es=sub),          # sub-vis, outside FV
          mk(z=Z_HIGH, ec=sub,    es=sub)]          # sub-vis, outside FV
    @test classify_event(:completed, empty_stack(), cs, params) === :SS_in_ROI
end

@testset "18. SS in FV but es outside ROI + sub-vis clusters → :SS_outside_ROI" begin
    sub = params.E_visible_keV / 1000.0 * 0.5
    cs = [mk(z=Z_MID, ec=Q_MeV, es=Q_MeV + 5*HW_MeV),
          mk(z=Z_HIGH, ec=sub,   es=sub)]
    @test classify_event(:completed, empty_stack(), cs, params) === :SS_outside_ROI
end

@testset "19. 0 visible clusters (all sub-vis) → :escaped" begin
    sub = params.E_visible_keV / 1000.0 * 0.5
    cs = [mk(z=Z_MID, ec=sub, es=sub),
          mk(z=Z_LOW, ec=sub, es=sub)]
    # No visible signal at all — treat as escaped (no detectable activity).
    @test classify_event(:completed, empty_stack(), cs, params) === :escaped
end

@testset "20. 2 visible clusters + many sub-vis → :MS_rejected" begin
    sub = params.E_visible_keV / 1000.0 * 0.5
    cs = [mk(z=30.0, ec=Q_MeV, es=Q_MeV),
          mk(z=40.0, ec=sub,    es=sub),
          mk(z=60.0, ec=0.5,    es=0.5),             # second visible
          mk(z=80.0, ec=sub,    es=sub)]
    @test classify_event(:completed, empty_stack(), cs, params) === :MS_rejected
end

@testset "21. visible cluster picked for ROI is the visible one, not clusters[1]" begin
    # Sub-visible cluster sorted first; visible cluster is later in the list.
    # Must use the VISIBLE cluster's es for the ROI decision.
    sub  = params.E_visible_keV / 1000.0 * 0.5
    # Sub-visible has es == sub (~5 keV), nowhere near ROI.
    # Visible has es == Q_MeV (in ROI).
    cs = [mk(z=Z_LOW,  ec=sub,    es=sub),       # NOT visible — would give :SS_outside_ROI if used
          mk(z=Z_MID,  ec=Q_MeV, es=Q_MeV)]      # visible + in ROI
    @test classify_event(:completed, empty_stack(), cs, params) === :SS_in_ROI
end

println("\n  ── test_classify_event3.jl: classify_event (with slow checks) OK ──\n")
