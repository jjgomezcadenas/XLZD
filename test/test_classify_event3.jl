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
    @test classify_event(:vetoed_skin, cs_in_fv, params) === :skin_vetoed
    @test classify_event(:vetoed_skin, cs_empty, params) === :skin_vetoed
    @test classify_event(:vetoed_skin, cs_two,   params) === :skin_vetoed
end

# ---------------------------------------------------------------------------
# 2. status = :rejected_fv → :SS_outside_FV (clusters ignored)
# ---------------------------------------------------------------------------

@testset "2. :rejected_fv → :SS_outside_FV (clusters ignored)" begin
    @test classify_event(:rejected_fv, [mk()],     params) === :SS_outside_FV
    @test classify_event(:rejected_fv, Cluster[],  params) === :SS_outside_FV
end

# ---------------------------------------------------------------------------
# 3. status = :escaped → :escaped (clusters ignored)
# ---------------------------------------------------------------------------

@testset "3. :escaped → :escaped" begin
    @test classify_event(:escaped, Cluster[], params) === :escaped
    @test classify_event(:escaped, [mk()],    params) === :escaped
end

# ---------------------------------------------------------------------------
# 4. status = :completed, empty clusters → :escaped
# ---------------------------------------------------------------------------

@testset "4. :completed + empty clusters → :escaped" begin
    @test classify_event(:completed, Cluster[], params) === :escaped
end

# ---------------------------------------------------------------------------
# 5. status = :completed, ≥ 2 clusters → :MS_rejected
# ---------------------------------------------------------------------------

@testset "5. :completed + 2 clusters → :MS_rejected" begin
    @test classify_event(:completed, [mk(z=40.0), mk(z=60.0)], params) === :MS_rejected
    @test classify_event(:completed, [mk(z=40.0), mk(z=60.0), mk(z=80.0)], params) === :MS_rejected
end

# ---------------------------------------------------------------------------
# 6. SS, cluster z below FV → :SS_outside_FV
# ---------------------------------------------------------------------------

@testset "6. SS cluster z < fv_z_min → :SS_outside_FV" begin
    @test classify_event(:completed, [mk(z=Z_LOW)], params) === :SS_outside_FV
end

# ---------------------------------------------------------------------------
# 7. SS, cluster z above FV → :SS_outside_FV
# ---------------------------------------------------------------------------

@testset "7. SS cluster z > fv_z_max → :SS_outside_FV" begin
    @test classify_event(:completed, [mk(z=Z_HIGH)], params) === :SS_outside_FV
end

# ---------------------------------------------------------------------------
# 8. SS, cluster outside r² → :SS_outside_FV
# ---------------------------------------------------------------------------

@testset "8. SS cluster r² > fv_r2_max → :SS_outside_FV" begin
    @test classify_event(:completed, [mk(x=R_OUT, y=0.0, z=Z_MID)], params) === :SS_outside_FV
end

# ---------------------------------------------------------------------------
# 9. SS in FV with es at Q → :SS_in_ROI
# ---------------------------------------------------------------------------

@testset "9. SS in FV, es == Q → :SS_in_ROI" begin
    @test classify_event(:completed, [mk(es=Q_MeV)], params) === :SS_in_ROI
end

# ---------------------------------------------------------------------------
# 10. SS in FV with es far from Q → :SS_outside_ROI
# ---------------------------------------------------------------------------

@testset "10. SS in FV, es far from Q → :SS_outside_ROI" begin
    @test classify_event(:completed, [mk(es=Q_MeV + 5*HW_MeV)], params) === :SS_outside_ROI
    @test classify_event(:completed, [mk(es=Q_MeV - 5*HW_MeV)], params) === :SS_outside_ROI
    @test classify_event(:completed, [mk(es=1.0)],              params) === :SS_outside_ROI
end

# ---------------------------------------------------------------------------
# 11. Coverage: every outcome (except :companion_vetoed) is reachable
# ---------------------------------------------------------------------------

@testset "11. coverage of CLASSIFY_EVENT_OUTCOMES (sans :companion_vetoed)" begin
    reached = Set{Symbol}()
    push!(reached, classify_event(:vetoed_skin, Cluster[], params))           # :skin_vetoed
    push!(reached, classify_event(:rejected_fv,   Cluster[], params))           # :SS_outside_FV
    push!(reached, classify_event(:completed,           Cluster[], params))           # :escaped
    push!(reached, classify_event(:completed, [mk(z=40.0), mk(z=60.0)], params))     # :MS_rejected
    push!(reached, classify_event(:completed, [mk(es=Q_MeV)], params))                # :SS_in_ROI
    push!(reached, classify_event(:completed, [mk(es=1.0)],   params))                # :SS_outside_ROI
    expected = setdiff(Set(CLASSIFY_EVENT_OUTCOMES), Set([:companion_vetoed]))
    # All six output symbols are reachable from the inputs above.
    @test reached == expected
end

# ---------------------------------------------------------------------------
# 12. Edge case: vetoed_skin takes priority over clusters
# ---------------------------------------------------------------------------

@testset "12. vetoed_skin priority over even SS-in-ROI clusters" begin
    # Even a perfect SS-in-ROI cluster cannot rescue a skin-vetoed event.
    @test classify_event(:vetoed_skin, [mk(es=Q_MeV)], params) === :skin_vetoed
end

println("\n  ── test_classify_event3.jl: classify_event OK ──\n")
