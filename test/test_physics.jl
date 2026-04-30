using Random

@testset "Cross-section interpolation at 2.448 MeV" begin
    xcom = load_xcom("data/nist.csv")

    # Expected values from spec (Section 4.1), 2% tolerance
    @test isapprox(σ_photo(xcom, 2.448),   8.75e-4; rtol=0.02)
    @test isapprox(σ_Compton(xcom, 2.448),  3.27e-2; rtol=0.02)
    @test isapprox(σ_pair(xcom, 2.448),     4.5e-3;  rtol=0.10)  # less precise spec value

    # Total (photo + compton + pair)
    σ_tot = σ_photo(xcom, 2.448) + σ_Compton(xcom, 2.448) + σ_pair(xcom, 2.448)
    @test isapprox(σ_tot, 3.81e-2; rtol=0.05)

    # μ_total_lin at LXe density
    μ = μ_total_lin(xcom, 2.448, 2.953)
    @test isapprox(μ, 0.113; rtol=0.05)

    # Pair production must be zero below threshold
    @test σ_pair(xcom, 0.5) == 0.0
    @test σ_pair(xcom, 1.0) == 0.0
end

@testset "Branching fractions at 2.448 MeV" begin
    xcom = load_xcom("data/nist.csv")
    pe = σ_photo(xcom, 2.448)
    co = σ_Compton(xcom, 2.448)
    pp = σ_pair(xcom, 2.448)
    tot = pe + co + pp

    f_photo   = pe / tot
    f_compton = co / tot
    f_pair    = pp / tot

    # Photoelectric ~2–3%
    @test 0.01 < f_photo < 0.04
    # Compton ~83–89%
    @test 0.83 < f_compton < 0.89
    # Pair ~9–15%
    @test 0.09 < f_pair < 0.15
    # Sum to 1
    @test isapprox(f_photo + f_compton + f_pair, 1.0; atol=1e-12)
end

@testset "Klein-Nishina: scattered energy < incident" begin
    rng = MersenneTwister(99)
    for E in [0.1, 0.5, 1.0, 2.448]
        for _ in 1:20
            E_scatt, cos_θ = sample_klein_nishina(rng, E)
            @test 0.0 < E_scatt < E
            @test -1.0 ≤ cos_θ ≤ 1.0
        end
    end
end

@testset "rotate_direction: preserves norm and angle" begin
    rng = MersenneTwister(77)
    for _ in 1:20
        # Random unit direction
        u = 2.0 * rand(rng) - 1.0
        φ = 2π * rand(rng)
        s = sqrt(1.0 - u * u)
        dx, dy, dz = s * cos(φ), s * sin(φ), u

        cos_θ = 2.0 * rand(rng) - 1.0
        ϕ = 2π * rand(rng)
        dx2, dy2, dz2 = rotate_direction(dx, dy, dz, cos_θ, ϕ)

        # Unit norm preserved
        norm2 = dx2^2 + dy2^2 + dz2^2
        @test isapprox(norm2, 1.0; atol=1e-12)

        # Dot product ≈ cos(θ)
        dot = dx * dx2 + dy * dy2 + dz * dz2
        @test isapprox(dot, cos_θ; atol=1e-10)
    end
end
