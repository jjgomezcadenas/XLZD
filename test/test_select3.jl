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
# select_SC — pure cluster-count predicate
# ---------------------------------------------------------------------------

@testset "1. select_SC: empty list → false" begin
    @test select_SC(Cluster[]) === false
end

@testset "2. select_SC: one cluster → true" begin
    @test select_SC([mk(Q_MeV, Q_MeV)]) === true
end

@testset "3. select_SC: two clusters → false" begin
    @test select_SC([mk(1.0, 1.0), mk(0.5, 0.5)]) === false
end

@testset "4. select_SC: three clusters → false" begin
    @test select_SC([mk(0.3, 0.3), mk(0.7, 0.7), mk(1.2, 1.2)]) === false
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

println("\n  ── test_select3.jl: select_SC + select_ROI OK ──\n")
