# test/test_stack_histograms3.jl — Unit tests for StackHistogramSet.
#
# Stripped to the three diagnostic histograms the cut-flow refactor
# keeps: interaction_type_freq, path_length_LXe, region_interaction.

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
# 1. Empty stack: only path_length entry (= 0) lands.
# ---------------------------------------------------------------------------

@testset "1. empty stack → path_length[bin 1] += 1; nothing else" begin
    sh = StackHistogramSet()
    s  = PhotonStack()                        # path_length_LXe = 0
    update_stack_histograms!(sh, s, params)
    @test sum(sh.interaction_type_freq)      == 0
    @test sum(sh.region_interaction_counts)  == 0
    @test sh.path_length_LXe_counts[1]        == 1
end

# ---------------------------------------------------------------------------
# 2. interaction_type_freq counts per row across all 4 types.
# ---------------------------------------------------------------------------

@testset "2. interaction_type_freq counts each row" begin
    sh = StackHistogramSet()
    s  = _stack_from([(:TPC, INT_PHOTO,   1.0),
                      (:TPC, INT_COMPTON, 0.5),
                      (:TPC, INT_COMPTON, 0.3),
                      (:TPC, INT_PAIR,    1.6),
                      (:Skin, INT_BELOW_THRESH, 0.02)])
    update_stack_histograms!(sh, s, params)
    @test sh.interaction_type_freq[1] == 1     # PHOTO
    @test sh.interaction_type_freq[2] == 2     # COMPTON
    @test sh.interaction_type_freq[3] == 1     # PAIR
    @test sh.interaction_type_freq[4] == 1     # BELOW_THRESH
    @test sum(sh.interaction_type_freq) == 5
end

# ---------------------------------------------------------------------------
# 3. region_interaction matrix: rows × cols = (TPC/Skin/Inert) × (P/C/Pair/BT).
# ---------------------------------------------------------------------------

@testset "3. region × interaction matrix" begin
    sh = StackHistogramSet()
    s  = _stack_from([(:TPC,   INT_COMPTON, 0.5),
                      (:Skin,  INT_COMPTON, 0.04),
                      (:Inert, INT_PHOTO,   0.6)])
    update_stack_histograms!(sh, s, params)
    @test sh.region_interaction_counts[1, 2] == 1   # TPC, COMPTON
    @test sh.region_interaction_counts[2, 2] == 1   # Skin, COMPTON
    @test sh.region_interaction_counts[3, 1] == 1   # Inert, PHOTO
    @test sum(sh.region_interaction_counts)  == 3
end

# ---------------------------------------------------------------------------
# 4. path_length_LXe is read off `stack.path_length_LXe`, not from rows.
# ---------------------------------------------------------------------------

@testset "4. path_length_LXe binned from stack accumulator" begin
    sh = StackHistogramSet(path_length_n_bins=10, path_length_max_cm=100.0)
    s  = PhotonStack()
    s.path_length_LXe = 35.0                 # bin (35/100*10)+1 = 4
    update_stack_histograms!(sh, s, params)
    @test sh.path_length_LXe_counts[4] == 1
    @test sum(sh.path_length_LXe_counts) == 1
end

@testset "5. path_length_LXe out-of-range silently skipped" begin
    sh = StackHistogramSet(path_length_n_bins=10, path_length_max_cm=10.0)
    s  = PhotonStack()
    s.path_length_LXe = 100.0                 # past max
    update_stack_histograms!(sh, s, params)
    @test sum(sh.path_length_LXe_counts) == 0
end

# ---------------------------------------------------------------------------
# 6. merge_stack_histograms! sums every field.
# ---------------------------------------------------------------------------

@testset "6. merge_stack_histograms! sums elementwise" begin
    sh1 = StackHistogramSet()
    sh2 = StackHistogramSet()
    update_stack_histograms!(sh1, _stack_from([(:TPC, INT_PHOTO,   1.0)]), params)
    update_stack_histograms!(sh2, _stack_from([(:TPC, INT_COMPTON, 0.5)]), params)
    merge_stack_histograms!(sh1, sh2)
    @test sh1.interaction_type_freq[1] == 1                   # PHOTO from sh1
    @test sh1.interaction_type_freq[2] == 1                   # COMPTON from sh2
    @test sum(sh1.region_interaction_counts) == 2
    @test sum(sh1.path_length_LXe_counts)    == 2             # both contributed bin 1
end

@testset "7. merge_stack_histograms! asserts compatible binning" begin
    sh1 = StackHistogramSet(path_length_n_bins=100)
    sh2 = StackHistogramSet(path_length_n_bins=50)            # mismatched
    @test_throws AssertionError merge_stack_histograms!(sh1, sh2)
end

println("\n  ── test_stack_histograms3.jl: StackHistogramSet OK ──\n")
