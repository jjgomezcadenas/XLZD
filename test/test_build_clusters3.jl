# test/test_build_clusters3.jl — build_clusters(rng, stack, params) → Vector{Cluster}.
#                                  Builds clusters; does NOT classify the event.
#
# Hand-built PhotonStack inputs, no MC, no tracking. Each testset
# constructs a stack with push_row! and asserts the cluster count,
# cluster energies, and energy-weighted centroids. Smearing is exercised
# in the final testset (statistical check on `cluster.es`).

using Test
using Random
using Printf
using Statistics

isdefined(Main, :XLZD3) || include("../src3/XLZD3.jl")
using .XLZD3

const params = MCParams()  # Δz_threshold_mm = 3.0 → 0.30 cm
const Δz     = Δz_threshold_cm(params)

# Helper: push a generic row with sensible defaults.
function pushrow!(s::PhotonStack; nm=0, parent_region=:CB, region=:TPC,
                   interaction=INT_PHOTO, x=0.0, y=0.0, z=0.0, epre=1.0, edep=1.0)
    push_row!(s; nm=nm, parent_region=parent_region, region=region,
                  interaction=interaction, x=x, y=y, z=z, epre=epre, edep=edep)
end

# Fixed-seed RNG used in all deterministic-shape testsets. Smearing
# values are not asserted in cases 1–8; testset 9 covers smearing.
fresh_rng() = MersenneTwister(0xABCD)

# ---------------------------------------------------------------------------
# 1. Empty stack → 0 clusters
# ---------------------------------------------------------------------------

@testset "1. empty stack → 0 clusters" begin
    s = PhotonStack()
    cs = build_clusters(fresh_rng(), s, params)
    @test cs isa Vector{Cluster}
    @test isempty(cs)
end

# ---------------------------------------------------------------------------
# 2. Single :TPC :PHOTO row → 1 cluster, ec = edep, centroid = (x,y,z)
# ---------------------------------------------------------------------------

@testset "2. single :TPC :PHOTO row → 1 cluster" begin
    s = PhotonStack()
    pushrow!(s; region=:TPC, interaction=INT_PHOTO,
              x=1.0, y=2.0, z=50.0, epre=2.448, edep=2.448)
    cs = build_clusters(fresh_rng(), s, params)
    @test length(cs) == 1
    @test cs[1].ec ≈ 2.448
    @test cs[1].xc ≈ 1.0
    @test cs[1].yc ≈ 2.0
    @test cs[1].zc ≈ 50.0
    @test cs[1].es != 0.0           # smearing happened
end

# ---------------------------------------------------------------------------
# 3. Two :TPC rows Δz = 1 mm → 1 cluster, ec summed, energy-weighted centroid
# ---------------------------------------------------------------------------

@testset "3. two :TPC rows Δz < 3mm → 1 cluster, ec summed" begin
    s = PhotonStack()
    # row A at z=50.00, edep=0.5; row B at z=50.10 (1 mm), edep=0.3
    pushrow!(s; region=:TPC, interaction=INT_COMPTON,
              x=0.0, y=0.0, z=50.00, epre=2.0, edep=0.5)
    pushrow!(s; region=:TPC, interaction=INT_PHOTO,
              x=2.0, y=4.0, z=50.10, epre=1.5, edep=0.3)
    cs = build_clusters(fresh_rng(), s, params)
    @test length(cs) == 1
    @test cs[1].ec ≈ 0.8
    # Energy-weighted: x = (0·0.5 + 2·0.3) / 0.8 = 0.75
    #                  y = (0·0.5 + 4·0.3) / 0.8 = 1.50
    #                  z = (50·0.5 + 50.1·0.3) / 0.8 = 50.0375
    @test cs[1].xc ≈ 0.75
    @test cs[1].yc ≈ 1.5
    @test cs[1].zc ≈ 50.0375
end

# ---------------------------------------------------------------------------
# 4. Two :TPC rows Δz = 5 cm → 2 clusters in z order
# ---------------------------------------------------------------------------

@testset "4. two :TPC rows Δz > 3mm → 2 clusters" begin
    s = PhotonStack()
    pushrow!(s; region=:TPC, x=0.0, y=0.0, z=50.0, edep=0.5)
    pushrow!(s; region=:TPC, x=0.0, y=0.0, z=55.0, edep=0.7)   # Δz = 5 cm
    cs = build_clusters(fresh_rng(), s, params)
    @test length(cs) == 2
    @test cs[1].zc ≈ 50.0
    @test cs[1].ec ≈ 0.5
    @test cs[2].zc ≈ 55.0
    @test cs[2].ec ≈ 0.7
end

# ---------------------------------------------------------------------------
# 5. Mix of :TPC + :Skin + :Inert → only :TPC contribute
# ---------------------------------------------------------------------------

@testset "5. mixed regions: only :TPC clustered" begin
    s = PhotonStack()
    pushrow!(s; region=:Skin,  x=0.0, y=0.0, z=10.0, edep=0.4)
    pushrow!(s; region=:Inert, x=0.0, y=0.0, z=20.0, edep=0.3)
    pushrow!(s; region=:TPC,   x=0.0, y=0.0, z=50.0, edep=0.5)
    pushrow!(s; region=:Gas,   x=0.0, y=0.0, z=80.0, edep=0.1)
    cs = build_clusters(fresh_rng(), s, params)
    @test length(cs) == 1
    @test cs[1].ec ≈ 0.5
    @test cs[1].zc ≈ 50.0
end

# ---------------------------------------------------------------------------
# 6. :BELOW_THRESH in :TPC contributes; :BELOW_THRESH outside :TPC excluded
# ---------------------------------------------------------------------------

@testset "6. :BELOW_THRESH respects region filter" begin
    s = PhotonStack()
    # In :TPC: contributes
    pushrow!(s; region=:TPC, interaction=INT_BELOW_THRESH,
              x=0.0, y=0.0, z=50.0, epre=0.020, edep=0.020)
    # In :Inert: filtered out
    pushrow!(s; region=:Inert, interaction=INT_BELOW_THRESH,
              x=0.0, y=0.0, z=51.0, epre=0.015, edep=0.015)
    cs = build_clusters(fresh_rng(), s, params)
    @test length(cs) == 1
    @test cs[1].ec ≈ 0.020
    @test cs[1].zc ≈ 50.0
end

# ---------------------------------------------------------------------------
# 7. Pair vertex + 2× 511 keV photo, all within Δz < 3 mm → 1 cluster
# ---------------------------------------------------------------------------

@testset "7. pair vertex + two 511 keV all close → 1 cluster" begin
    s = PhotonStack()
    # Pair vertex: 2.615 MeV → kinetic = 2.615 - 1.022 = 1.593 MeV
    pair_ng = push_row!(s; nm=0, parent_region=:CTH, region=:TPC,
                          interaction=INT_PAIR,
                          x=0.0, y=0.0, z=50.00, epre=2.615, edep=1.593)
    # Both 511 keV annihilation γ photoelectric-absorb within Δz < 3 mm
    push_row!(s; nm=pair_ng, parent_region=:TPC, region=:TPC,
                interaction=INT_PHOTO,
                x=0.05, y=0.0, z=50.05, epre=0.511, edep=0.511)
    push_row!(s; nm=pair_ng, parent_region=:TPC, region=:TPC,
                interaction=INT_PHOTO,
                x=-0.05, y=0.0, z=50.10, epre=0.511, edep=0.511)
    cs = build_clusters(fresh_rng(), s, params)
    @test length(cs) == 1
    @test cs[1].ec ≈ 1.593 + 2*0.511   # = 2.615 MeV total visible
end

# ---------------------------------------------------------------------------
# 8. Three rows in z order with mixed Δz → 2 clusters with right energy split
# ---------------------------------------------------------------------------

@testset "8. mixed Δz pattern (1mm, 5cm) → 2 clusters" begin
    s = PhotonStack()
    pushrow!(s; region=:TPC, x=0.0, y=0.0, z=50.00, edep=0.4)
    pushrow!(s; region=:TPC, x=0.0, y=0.0, z=50.10, edep=0.6)   # +1 mm → joins cluster 1
    pushrow!(s; region=:TPC, x=0.0, y=0.0, z=55.10, edep=0.9)   # +5 cm → starts cluster 2
    cs = build_clusters(fresh_rng(), s, params)
    @test length(cs) == 2
    @test cs[1].ec ≈ 1.0     # 0.4 + 0.6
    @test cs[2].ec ≈ 0.9
    # Energy-weighted z for cluster 1: (50·0.4 + 50.1·0.6)/1.0 = 50.06
    @test cs[1].zc ≈ 50.06
    @test cs[2].zc ≈ 55.10
end

# ---------------------------------------------------------------------------
# 9. Smearing distribution: ec at Q_ββ → mean(es) ≈ Q, std(es) ≈ σ_E·Q
# ---------------------------------------------------------------------------

@testset "9. smearing: mean(es) ≈ ec, std(es) ≈ σ_E·ec" begin
    Q_MeV    = params.Q_betabeta_keV / 1000.0
    σ_E_MeV  = params.σ_E_over_E * Q_MeV
    rng = MersenneTwister(0x1234)
    N = 5000
    es_samples = Float64[]
    for _ in 1:N
        s = PhotonStack()
        pushrow!(s; region=:TPC, x=0.0, y=0.0, z=50.0, epre=Q_MeV, edep=Q_MeV)
        cs = build_clusters(rng, s, params)
        @assert length(cs) == 1
        push!(es_samples, cs[1].es)
    end
    μ = mean(es_samples)
    σ = std(es_samples)
    @printf("\n     smearing at Q_ββ: μ=%.4f MeV, σ=%.4f MeV  (expected μ=%.4f, σ=%.4f)\n",
            μ, σ, Q_MeV, σ_E_MeV)
    # Allow ±3·SE/√N on the mean and ±10% on σ
    @test abs(μ - Q_MeV)  < 3 * σ_E_MeV / sqrt(N)
    @test 0.90 * σ_E_MeV < σ < 1.10 * σ_E_MeV
end

println("\n  ── test_build_clusters3.jl: build_clusters OK ──\n")
