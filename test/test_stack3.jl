# test/test_stack3.jl — StackRow / PhotonStack unit tests.

using Test

isdefined(Main, :XLZD3) || include("../src3/XLZD3.jl")
using .XLZD3

# ---------------------------------------------------------------------------
# StackRow basic construction
# ---------------------------------------------------------------------------

@testset "StackRow construction and field access" begin
    r = StackRow(1, 0, :CB, :TPC, INT_PAIR,
                 1.0, 2.0, 3.0, 2.615, 1.593)
    @test r.ng           == 1
    @test r.nm           == 0
    @test r.parent_region === :CB
    @test r.region       === :TPC
    @test r.interaction  === INT_PAIR
    @test r.x            == 1.0
    @test r.y            == 2.0
    @test r.z            == 3.0
    @test r.epre         == 2.615
    @test r.edep         == 1.593
end

# ---------------------------------------------------------------------------
# Interaction symbol constants
# ---------------------------------------------------------------------------

@testset "Interaction constants" begin
    @test INT_PHOTO        === :PHOTO
    @test INT_COMPTON      === :COMPTON
    @test INT_PAIR         === :PAIR
    @test INT_BELOW_THRESH === :BELOW_THRESH
    @test length(unique([INT_PHOTO, INT_COMPTON, INT_PAIR, INT_BELOW_THRESH])) == 4
end

# ---------------------------------------------------------------------------
# PhotonStack defaults
# ---------------------------------------------------------------------------

@testset "PhotonStack defaults" begin
    s = PhotonStack()
    @test length(s)   == 0
    @test s.next_ng   == 1
    @test isempty(s.rows)
end

# ---------------------------------------------------------------------------
# push_row! — sequential ng assignment
# ---------------------------------------------------------------------------

@testset "push_row! assigns sequential ng" begin
    s = PhotonStack()
    ng1 = push_row!(s; nm=0, parent_region=:CB,  region=:TPC, interaction=INT_PAIR,
                       x=0.0, y=0.0, z=0.0, epre=2.615, edep=1.593)
    ng2 = push_row!(s; nm=1, parent_region=:TPC, region=:TPC, interaction=INT_COMPTON,
                       x=1.0, y=0.0, z=0.0, epre=0.511, edep=0.230)
    ng3 = push_row!(s; nm=2, parent_region=:TPC, region=:TPC, interaction=INT_PHOTO,
                       x=2.0, y=0.0, z=0.0, epre=0.281, edep=0.281)
    @test ng1 == 1
    @test ng2 == 2
    @test ng3 == 3
    @test length(s)  == 3
    @test s.next_ng  == 4
    @test s.rows[1].ng == 1
    @test s.rows[2].ng == 2
    @test s.rows[3].ng == 3
end

# ---------------------------------------------------------------------------
# push_row! preserves all fields
# ---------------------------------------------------------------------------

@testset "push_row! preserves all fields" begin
    s = PhotonStack()
    push_row!(s; nm=7, parent_region=:Inert, region=:Skin, interaction=INT_BELOW_THRESH,
                x=-3.5, y=12.25, z=42.0, epre=0.025, edep=0.025)
    r = s.rows[1]
    @test r.ng            == 1
    @test r.nm            == 7
    @test r.parent_region === :Inert
    @test r.region        === :Skin
    @test r.interaction   === INT_BELOW_THRESH
    @test r.x             == -3.5
    @test r.y             == 12.25
    @test r.z             == 42.0
    @test r.epre          == 0.025
    @test r.edep          == 0.025
end

# ---------------------------------------------------------------------------
# empty! resets ng counter and clears rows
# ---------------------------------------------------------------------------

@testset "empty! resets state" begin
    s = PhotonStack()
    push_row!(s; nm=0, parent_region=:CB, region=:TPC, interaction=INT_PHOTO,
                x=0.0, y=0.0, z=0.0, epre=2.448, edep=2.448)
    push_row!(s; nm=0, parent_region=:CB, region=:TPC, interaction=INT_PHOTO,
                x=0.0, y=0.0, z=0.0, epre=2.448, edep=2.448)
    push_row!(s; nm=0, parent_region=:CB, region=:TPC, interaction=INT_PHOTO,
                x=0.0, y=0.0, z=0.0, epre=2.448, edep=2.448)
    @test length(s) == 3
    @test s.next_ng == 4

    empty!(s)
    @test length(s) == 0
    @test s.next_ng == 1
    @test isempty(s.rows)

    # Subsequent push starts ng numbering over.
    ng = push_row!(s; nm=0, parent_region=:CB, region=:TPC, interaction=INT_PHOTO,
                      x=0.0, y=0.0, z=0.0, epre=1.0, edep=1.0)
    @test ng == 1
    @test s.next_ng == 2
end

# ---------------------------------------------------------------------------
# Mother-chain consistency on a hand-built stack
# ---------------------------------------------------------------------------

@testset "Mother-chain consistency: pair → 511 → Compton → photo" begin
    # Source photon (2.615 MeV from Tl-208 in :CTH) makes a pair in :TPC.
    # One of the 511 keV annihilation γ Compton-scatters in :TPC, then
    # the scattered photon photoelectric-absorbs in :TPC.
    s = PhotonStack()
    pair_ng = push_row!(s; nm=0, parent_region=:CTH, region=:TPC,
                          interaction=INT_PAIR,
                          x=0.0, y=0.0, z=50.0, epre=2.615, edep=1.593)
    comp_ng = push_row!(s; nm=pair_ng, parent_region=:TPC, region=:TPC,
                          interaction=INT_COMPTON,
                          x=0.5, y=0.0, z=50.4, epre=0.511, edep=0.230)
    phot_ng = push_row!(s; nm=comp_ng, parent_region=:TPC, region=:TPC,
                          interaction=INT_PHOTO,
                          x=1.0, y=0.0, z=50.7, epre=0.281, edep=0.281)

    @test pair_ng == 1
    @test comp_ng == 2
    @test phot_ng == 3
    @test s.rows[1].nm == 0
    @test s.rows[2].nm == 1
    @test s.rows[3].nm == 2

    # Walk the chain from row 3 up to source (nm = 0).
    chain = Int[]
    cur = phot_ng
    while cur != 0
        push!(chain, cur)
        cur = s.rows[cur].nm
    end
    @test chain == [3, 2, 1]   # reverse order from leaf to source

    # Energy bookkeeping for the 511 keV branch:
    # Compton epre - edep = 281 keV, which equals next photo's epre.
    @test isapprox(s.rows[2].epre - s.rows[2].edep, s.rows[3].epre; atol=1e-9)
end

println("\n  ── test_stack3.jl: StackRow / PhotonStack OK ──\n")
