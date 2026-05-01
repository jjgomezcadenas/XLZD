# test/test_stack_histograms3.jl — Unit tests for StackHistogramSet.
#
# Hand-built PhotonStack inputs, no MC. Each testset constructs a
# stack via push_row! and asserts the resulting bin counts in the
# StackHistogramSet.

using Test

isdefined(Main, :XLZD3) || include("../src3/XLZD3.jl")
using .XLZD3

const params = MCParams()

# Helper: build a stack from a list of (region, interaction, edep) triples
# at the same (x, y, z) position by default.
function _stack_from(rows::Vector{Tuple{Symbol,Symbol,Float64}};
                      x=0.0, y=0.0, z=50.0)
    s = PhotonStack()
    for (reg, inter, ed) in rows
        push_row!(s; nm=0, parent_region=:CB, region=reg,
                      interaction=inter, x=x, y=y, z=z,
                      epre=ed, edep=ed)
    end
    s
end

# ---------------------------------------------------------------------------
# 1. Empty stack: ng_max bin 1 (=0 rows) increments; per-event totals = 0
# ---------------------------------------------------------------------------

@testset "1. empty stack → ng_max[1] += 1, all per-event counts in bin 1" begin
    sh = StackHistogramSet()
    s  = PhotonStack()
    update_stack_histograms!(sh, s, params)
    @test sh.ng_max_counts[1] == 1
    @test sh.n_photo_counts[1]        == 1
    @test sh.n_compton_counts[1]      == 1
    @test sh.n_pair_counts[1]         == 1
    @test sh.n_below_thresh_counts[1] == 1
    @test sum(sh.first_interaction_counts) == 0
    @test sum(sh.inclusive_edep_counts)    == 0
    @test sum(sh.E_first_counts)           == 0
    @test sum(sh.Δz_counts)                == 0
end

# ---------------------------------------------------------------------------
# 2. Single :TPC PHOTO row: first_interaction[1] = 1; n_photo[2] = 1
# ---------------------------------------------------------------------------

@testset "2. single :TPC PHOTO" begin
    sh = StackHistogramSet()
    s  = _stack_from([(:TPC, INT_PHOTO, 2.448)])
    update_stack_histograms!(sh, s, params)
    @test sh.first_interaction_counts[1]  == 1   # PHOTO
    @test sh.n_photo_counts[2]            == 1   # 1 photo this event
    @test sh.n_compton_counts[1]          == 1   # 0 compton
    @test sh.E_first_counts[ floor(Int, 2.448 / 2.7 * 270) + 1 ] == 1
    @test sum(sh.Δz_counts) == 0   # only one :TPC row → no Δz pairs
    @test sh.region_interaction_counts[1, 1] == 1   # (:TPC, PHOTO)
end

# ---------------------------------------------------------------------------
# 3. PAIR + COMPTON sequence: first_interaction[3] (PAIR) = 1
# ---------------------------------------------------------------------------

@testset "3. PAIR-then-COMPTON" begin
    sh = StackHistogramSet()
    s  = _stack_from([(:TPC, INT_PAIR,    1.593),
                      (:TPC, INT_COMPTON, 0.300)])
    update_stack_histograms!(sh, s, params)
    @test sh.first_interaction_counts[3]  == 1   # PAIR
    @test sh.n_pair_counts[2]             == 1
    @test sh.n_compton_counts[2]          == 1
    @test sh.region_interaction_counts[1, 3] == 1   # (:TPC, PAIR)
    @test sh.region_interaction_counts[1, 2] == 1   # (:TPC, COMPTON)
end

# ---------------------------------------------------------------------------
# 4. Δz: two :TPC rows at z = 50.0 and z = 50.10 → Δz = 0.1 cm in bin 1
# ---------------------------------------------------------------------------

@testset "4. Δz from first :TPC deposit" begin
    sh = StackHistogramSet()
    s  = PhotonStack()
    push_row!(s; nm=0, parent_region=:CB, region=:TPC, interaction=INT_COMPTON,
                  x=0.0, y=0.0, z=50.00, epre=2.0, edep=0.5)
    push_row!(s; nm=1, parent_region=:TPC, region=:TPC, interaction=INT_PHOTO,
                  x=0.0, y=0.0, z=50.10, epre=1.5, edep=1.5)
    update_stack_histograms!(sh, s, params)
    @test sh.Δz_counts[1] == 1   # 0.1 cm; bin 1 = [0, 0.5)
end

# ---------------------------------------------------------------------------
# 5. Per-region × per-interaction counts: skin Compton + inert photo
# ---------------------------------------------------------------------------

@testset "5. region × interaction matrix" begin
    sh = StackHistogramSet()
    s  = _stack_from([(:TPC,  INT_COMPTON, 0.5),
                      (:Skin, INT_COMPTON, 0.04),
                      (:Inert, INT_PHOTO,  0.6)])
    update_stack_histograms!(sh, s, params)
    @test sh.region_interaction_counts[1, 2] == 1   # TPC, COMPTON
    @test sh.region_interaction_counts[2, 2] == 1   # Skin, COMPTON
    @test sh.region_interaction_counts[3, 1] == 1   # Inert, PHOTO
    @test sum(sh.region_interaction_counts)  == 3
end

# ---------------------------------------------------------------------------
# 6. inclusive_edep: sum across all rows
# ---------------------------------------------------------------------------

@testset "6. inclusive_edep totals" begin
    sh = StackHistogramSet()
    s  = _stack_from([(:TPC, INT_COMPTON, 0.5),
                      (:TPC, INT_PHOTO,   0.3),
                      (:Skin, INT_COMPTON, 0.04)])
    update_stack_histograms!(sh, s, params)
    expected = 0.5 + 0.3 + 0.04
    bin = floor(Int, expected / 2.7 * 270) + 1
    @test sh.inclusive_edep_counts[bin] == 1
    @test sum(sh.inclusive_edep_counts) == 1
end

# ---------------------------------------------------------------------------
# 7. ng_max counts the number of rows in the event (in bin = n+1)
# ---------------------------------------------------------------------------

@testset "7. ng_max bins" begin
    sh = StackHistogramSet()
    s  = _stack_from([(:TPC, INT_COMPTON, 0.1),
                      (:TPC, INT_COMPTON, 0.1),
                      (:TPC, INT_PHOTO,   0.1)])
    update_stack_histograms!(sh, s, params)
    @test sh.ng_max_counts[4] == 1     # 3 rows → bin 3+1 = 4
end

# ---------------------------------------------------------------------------
# 8. merge_stack_histograms! sums counts elementwise
# ---------------------------------------------------------------------------

@testset "8. merge_stack_histograms!" begin
    sh1 = StackHistogramSet()
    sh2 = StackHistogramSet()
    update_stack_histograms!(sh1, _stack_from([(:TPC, INT_PHOTO, 1.0)]), params)
    update_stack_histograms!(sh2, _stack_from([(:TPC, INT_PAIR,  1.5)]), params)
    merge_stack_histograms!(sh1, sh2)
    @test sh1.first_interaction_counts[1] == 1   # PHOTO from sh1
    @test sh1.first_interaction_counts[3] == 1   # PAIR from sh2
    @test sum(sh1.first_interaction_counts) == 2
end

println("\n  ── test_stack_histograms3.jl: StackHistogramSet OK ──\n")
