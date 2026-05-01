# test/test_physics3.jl — Verify the ported physics primitives.

using Test
using Random
using Printf
using Statistics

isdefined(Main, :XLZD3) || include("../src3/XLZD3.jl")
using .XLZD3

const xcom_path = joinpath(@__DIR__, "..", "data", "nist.csv")
@assert isfile(xcom_path)

xcom = load_xcom(xcom_path)
const ρ_LXE = 2.953  # g/cm³

# ---------------------------------------------------------------------------
# XCOM loader & cross sections
# ---------------------------------------------------------------------------

@testset "XCOMTable load" begin
    @test length(xcom.energy_MeV) > 5
    @test issorted(xcom.energy_MeV)
    @test all(xcom.σ_photo .>= 0)
    @test all(xcom.σ_Compton .>= 0)
    @test all(xcom.σ_pair .>= 0)
    # Log grids are consistent with the linear ones
    @test xcom.log_E ≈ log.(xcom.energy_MeV)
end

@testset "Cross-section accessors at γ-relevant energies" begin
    for E in (0.583, 2.448, 2.615)
        sP = σ_photo(xcom, E)
        sC = σ_Compton(xcom, E)
        sX = σ_pair(xcom, E)
        @test sP > 0
        @test sC > 0
        @test sX >= 0  # may be 0 below 1.022 MeV
        # Compton dominant in this regime for LXe
        @test sC > sP

        μ = μ_total_lin(xcom, E, ρ_LXE)
        @test μ ≈ (sP + sC + sX) * ρ_LXE
        # Sanity: μ at 2.6 MeV in LXe ~ 0.09–0.12 cm⁻¹
        if E > 2.0
            @test 0.07 < μ < 0.15
        end
    end

    # σ_pair below threshold returns 0 exactly
    @test σ_pair(xcom, 0.5) == 0.0
    @test σ_pair(xcom, 1.0) == 0.0   # 1.0 < 1.022
    # Above threshold
    @test σ_pair(xcom, 2.0) > 0
end

# ---------------------------------------------------------------------------
# Klein-Nishina kinematics + distribution
# ---------------------------------------------------------------------------

@testset "Klein-Nishina kinematics" begin
    rng = MersenneTwister(0xCAFE)
    E0 = 2.448
    # Sample many; verify each (E', cos θ) satisfies 1/E' − 1/E = (1−cosθ)/m_e c²
    for _ in 1:1000
        E_prime, cos_θ = sample_klein_nishina(rng, E0)
        @test 0 < E_prime ≤ E0
        @test -1 ≤ cos_θ ≤ 1
        # Kinematic constraint
        lhs = 1.0 / E_prime - 1.0 / E0
        rhs = (1.0 - cos_θ) / ME_C2_MEV
        @test isapprox(lhs, rhs; atol=1e-10)
    end
end

@testset "Klein-Nishina mean energy and forward bias at 2.6 MeV" begin
    rng = MersenneTwister(0xBEEF)
    E0 = 2.615
    N = 100_000
    cos_θs = Float64[]
    E_primes = Float64[]
    for _ in 1:N
        E_p, cos_θ = sample_klein_nishina(rng, E0)
        push!(cos_θs, cos_θ)
        push!(E_primes, E_p)
    end
    # At 2.6 MeV the KN cross section is forward-biased: ⟨cos θ⟩ ≈ 0.45.
    # Just check it's clearly positive (forward) and below the low-energy
    # value (~0.7 at 0.1 MeV).
    @test 0.3 < mean(cos_θs) < 0.6
    # Mean scattered energy is well below the input
    @test mean(E_primes) < E0
    @test mean(E_primes) > E0 / 4   # not collapsed to 0
end

# ---------------------------------------------------------------------------
# rotate_direction
# ---------------------------------------------------------------------------

@testset "rotate_direction — identity, normalisation, no rotation" begin
    rng = MersenneTwister(0xAA)
    # Identity: θ = 0 (cos θ = 1, any φ) → output = input
    for _ in 1:50
        # random unit vector
        u = 2*rand(rng) - 1
        s = sqrt(1 - u*u)
        φ = 2π * rand(rng)
        d = (s*cos(φ), s*sin(φ), u)
        for φ_rot in (0.0, 1.5, 4.7)
            d2 = rotate_direction(d..., 1.0, φ_rot)
            @test all(isapprox.(d2, d; atol=1e-12))
        end
    end
end

@testset "rotate_direction — output is a unit vector" begin
    rng = MersenneTwister(0xBB)
    for _ in 1:200
        u = 2*rand(rng) - 1
        s = sqrt(1 - u*u)
        φ_pos = 2π * rand(rng)
        d = (s*cos(φ_pos), s*sin(φ_pos), u)
        cos_θ = 2*rand(rng) - 1
        φ = 2π * rand(rng)
        dx, dy, dz = rotate_direction(d..., cos_θ, φ)
        @test isapprox(dx*dx + dy*dy + dz*dz, 1.0; atol=1e-10)
    end
end

@testset "rotate_direction — angle preserved with respect to input" begin
    rng = MersenneTwister(0xCC)
    for _ in 1:100
        u = 2*rand(rng) - 1
        s = sqrt(1 - u*u)
        φ_pos = 2π * rand(rng)
        d = (s*cos(φ_pos), s*sin(φ_pos), u)
        cos_θ_target = 2*rand(rng) - 1
        φ = 2π * rand(rng)
        d2 = rotate_direction(d..., cos_θ_target, φ)
        # Dot product of input with output equals cos θ
        dot = d[1]*d2[1] + d[2]*d2[2] + d[3]*d2[3]
        @test isapprox(dot, cos_θ_target; atol=1e-10)
    end
end

@testset "rotate_direction — pole stability (input near +ẑ)" begin
    # Inputs nearly aligned with ẑ used to take the second basis branch
    d = (1e-5, 0.0, 1.0)
    norm_d = sqrt(d[1]^2 + d[2]^2 + d[3]^2)
    d_unit = d ./ norm_d
    d2 = rotate_direction(d_unit..., 0.5, 1.0)
    @test isapprox(d2[1]^2 + d2[2]^2 + d2[3]^2, 1.0; atol=1e-10)

    d = (0.0, 0.0, 1.0)
    d2 = rotate_direction(d..., 0.0, 0.0)
    @test isapprox(d2[1]^2 + d2[2]^2 + d2[3]^2, 1.0; atol=1e-10)
end

# ---------------------------------------------------------------------------
# Documentation
# ---------------------------------------------------------------------------

println()
println("── XCOM cross sections (LXe) ──")
@printf("  E (MeV)   σ_photo    σ_Compton  σ_pair    μ_total (cm⁻¹)\n")
for E in (0.583, 1.0, 2.448, 2.615)
    @printf("  %5.3f    %8.3e   %8.3e  %8.3e   %.4f\n",
            E, σ_photo(xcom, E), σ_Compton(xcom, E),
            σ_pair(xcom, E), μ_total_lin(xcom, E, ρ_LXE))
end
println()
