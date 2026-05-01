# test/test_sampling3.jl — verify the source-sampling primitives.

using Test
using Random
using Printf
using Statistics

isdefined(Main, :XLZD3) || include("../src3/XLZD3.jl")
using .XLZD3

# Helpers used in expectations
_trapz(y, x) = sum(0.5 * (y[i]+y[i+1]) * (x[i+1]-x[i]) for i in 1:length(x)-1)

# Build a detector and a real EffectiveSource set for tests
const lxe_csv      = joinpath(@__DIR__, "..", "data", "lxe_detector.csv")
const lxe_nist     = joinpath(@__DIR__, "..", "data", "nist_lxe.csv")
const ti_path      = joinpath(@__DIR__, "..", "data", "nist_ti.csv")
const geom_csv     = joinpath(@__DIR__, "..", "data", "lz_cryo_geometry.csv")
const extras_csv   = joinpath(@__DIR__, "..", "data", "lz_cryo_extras.csv")
const surfaces_csv = joinpath(@__DIR__, "..", "data", "lz_cryo_surface_sources.csv")

mat_LXe = load_material("LXe", 2.953, lxe_nist)
mat_Ti  = load_material("Ti",  4.510, ti_path)
det     = build_lxe_detector(lxe_csv, mat_LXe)
cryo    = build_cryostat(geom_csv, extras_csv, surfaces_csv)
indiv   = build_individual_sources(cryo, mat_Ti)
effs    = build_effective_sources(indiv, cryo, mat_Ti)
by_name = Dict(e.name => e for e in effs)

const RNG_SEED = 0x5AFE5A
const N_SAMPLE = 200_000

# ---------------------------------------------------------------------------
# build_cdf + sample_u
# ---------------------------------------------------------------------------

@testset "build_cdf — uniform dN/du" begin
    u = collect(range(0.0, 1.0, length=21))
    f = ones(length(u))
    cdf = build_cdf(u, f)
    @test cdf[1]   == 0.0
    @test cdf[end] == 1.0
    @test all(diff(cdf) .>= -1e-15)   # non-decreasing
    # CDF of a uniform PDF on [0, 1] is just u itself
    @test cdf ≈ u rtol=1e-10
end

@testset "sample_u reproduces a uniform input dN/du" begin
    rng = MersenneTwister(RNG_SEED)
    u_bins = collect(range(0.0, 1.0, length=51))
    f = ones(length(u_bins))
    cdf = build_cdf(u_bins, f)
    samples = [sample_u(rng, u_bins, cdf) for _ in 1:N_SAMPLE]
    @test 0.0 ≤ minimum(samples) ≤ 0.05
    @test 0.95 ≤ maximum(samples) ≤ 1.0
    # Uniform on [0,1] → mean = 0.5
    @test isapprox(mean(samples), 0.5, atol=0.005)
end

@testset "sample_u reproduces a linear input dN/du = (a + b·u)" begin
    rng = MersenneTwister(RNG_SEED + 1)
    u_bins = collect(range(0.0, 1.0, length=101))
    a, b = 0.4, 1.2
    f = a .+ b .* u_bins              # > 0 on [0, 1]
    cdf = build_cdf(u_bins, f)
    samples = [sample_u(rng, u_bins, cdf) for _ in 1:N_SAMPLE]
    # ⟨u⟩ for a linear PDF p(u) ∝ a + b·u on [0, 1] is
    #     (a/2 + b/3) / (a + b/2)
    expected_mean = (a/2 + b/3) / (a + b/2)
    @test isapprox(mean(samples), expected_mean, atol=0.005)
end

@testset "sample_u reproduces ⟨u⟩ from a real EffectiveSource" begin
    rng = MersenneTwister(RNG_SEED + 2)
    eff = by_name["CB_Bi214"]
    cdf = build_cdf(eff.u_bins, eff.dNdu)
    expected_mean = _trapz(eff.u_bins .* eff.dNdu, eff.u_bins) /
                    _trapz(eff.dNdu, eff.u_bins)
    samples = [sample_u(rng, eff.u_bins, cdf) for _ in 1:N_SAMPLE]
    @test isapprox(mean(samples), expected_mean, rtol=0.01)
    @test all(eff.u_bins[1] .<= samples .<= eff.u_bins[end])
end

# ---------------------------------------------------------------------------
# sample_barrel_entry
# ---------------------------------------------------------------------------

@testset "sample_barrel_entry — geometry and direction" begin
    rng = MersenneTwister(RNG_SEED + 10)
    eff = by_name["CB_Bi214"]
    cdf = build_cdf(eff.u_bins, eff.dNdu)
    R = det.R_ICV_inner

    # Sample many and check spatial constraints
    rs = Float64[]; zs = Float64[]; phis = Float64[]
    cosθs = Float64[]
    for _ in 1:N_SAMPLE
        x, y, z, dx, dy, dz = sample_barrel_entry(rng, det, cdf, eff.u_bins)
        push!(rs, sqrt(x*x + y*y))
        push!(zs, z)
        push!(phis, atan(y, x))
        # Direction is unit
        @assert isapprox(dx*dx+dy*dy+dz*dz, 1.0; atol=1e-10)
        # Direction points inward: dot with (-x/R, -y/R, 0) > 0
        push!(cosθs, -(dx*x + dy*y) / R)
    end
    @test all(rs .≈ R)
    @test all(det.z_RFR_bottom .≤ zs .≤ det.z_gate)
    # φ uniform → mean cos(φ) ≈ 0
    @test abs(mean(cos.(phis))) < 0.02
    @test abs(mean(sin.(phis))) < 0.02
    # Direction inward: mean cosθ matches the integrated ⟨u⟩ of the source
    expected_u = _trapz(eff.u_bins .* eff.dNdu, eff.u_bins) /
                 _trapz(eff.dNdu, eff.u_bins)
    @test isapprox(mean(cosθs), expected_u, rtol=0.02)
    @test all(cosθs .> 0)            # all inward (θ < π/2 from the inward normal)
end

# ---------------------------------------------------------------------------
# sample_endcap_entry
# ---------------------------------------------------------------------------

@testset "sample_endcap_entry — top head" begin
    rng = MersenneTwister(RNG_SEED + 20)
    eff = by_name["CTH_Bi214"]
    cdf = build_cdf(eff.u_bins, eff.dNdu)
    disk = icv_top_inner_disk(det)
    R = det.R_ICV_inner
    d = R / ICV_TOP_ASPECT
    z_eq = disk.z_equator

    # Sample many; verify ellipsoidal surface eq, inward normal, direction.
    cosθs = Float64[]
    for _ in 1:N_SAMPLE
        x, y, z, dx, dy, dz = sample_endcap_entry(rng, det, :endcap_top,
                                                   cdf, eff.u_bins)
        # Ellipsoid: (r/R)^2 + (z_off/d)^2 = 1, with z_off = z - z_eq ≥ 0 (head bulges up)
        r2 = x*x + y*y
        z_off = z - z_eq
        @test z_off ≥ -1e-9
        @test isapprox(r2/R^2 + (z_off/d)^2, 1.0; atol=1e-8)

        # Direction is unit
        @test isapprox(dx*dx+dy*dy+dz*dz, 1.0; atol=1e-10)

        # Direction inward
        nx, ny, nz = inward_normal(disk, x, y, z)
        push!(cosθs, dx*nx + dy*ny + dz*nz)
    end
    expected_u = _trapz(eff.u_bins .* eff.dNdu, eff.u_bins) /
                 _trapz(eff.dNdu, eff.u_bins)
    @test isapprox(mean(cosθs), expected_u, rtol=0.02)
    @test all(cosθs .> 0)
end

@testset "sample_endcap_entry — bottom head (3:1)" begin
    rng = MersenneTwister(RNG_SEED + 30)
    eff = by_name["CBH_Bi214"]
    cdf = build_cdf(eff.u_bins, eff.dNdu)
    disk = icv_bot_inner_disk(det)
    R = det.R_ICV_inner
    d = R / ICV_BOT_ASPECT
    z_eq = disk.z_equator

    cosθs = Float64[]
    for _ in 1:N_SAMPLE
        x, y, z, dx, dy, dz = sample_endcap_entry(rng, det, :endcap_bottom,
                                                   cdf, eff.u_bins)
        # Bottom head bulges down: z_off = z - z_eq ≤ 0
        r2 = x*x + y*y
        z_off = z - z_eq
        @test z_off ≤ 1e-9
        @test isapprox(r2/R^2 + (z_off/d)^2, 1.0; atol=1e-8)

        @test isapprox(dx*dx+dy*dy+dz*dz, 1.0; atol=1e-10)
        nx, ny, nz = inward_normal(disk, x, y, z)
        push!(cosθs, dx*nx + dy*ny + dz*nz)
    end
    expected_u = _trapz(eff.u_bins .* eff.dNdu, eff.u_bins) /
                 _trapz(eff.dNdu, eff.u_bins)
    @test isapprox(mean(cosθs), expected_u, rtol=0.02)
    @test all(cosθs .> 0)
end

# ---------------------------------------------------------------------------
# Dispatching sample_entry
# ---------------------------------------------------------------------------

@testset "sample_entry dispatches on EffectiveSource.region" begin
    rng = MersenneTwister(RNG_SEED + 40)
    for name in ("CB_Tl208", "CTH_Tl208", "CBH_Tl208")
        eff = by_name[name]
        x, y, z, dx, dy, dz = sample_entry(rng, det, eff)
        # Sanity: position inside ICV inner radius and direction unit
        @test sqrt(x*x + y*y) ≤ det.R_ICV_inner + 1e-6
        @test isapprox(dx*dx + dy*dy + dz*dz, 1.0; atol=1e-10)
    end
end

println()
@printf("  Drew %d samples from each of CB / CTH / CBH (Bi-214), ⟨cosθ⟩ matched\n",
        N_SAMPLE)
println()
