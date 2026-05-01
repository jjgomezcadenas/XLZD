# test/test_transport3.jl — Slab and transmission factor tests.

using Test

isdefined(Main, :XLZD3) || include("../src3/XLZD3.jl")
using .XLZD3

@testset "Slab basic" begin
    s = Slab(0.2, 1.0, "test")
    @test s.μ == 0.2
    @test s.t_cm == 1.0
    @test s.label == "test"
    @test optical_depth([s]) ≈ 0.2
    @test optical_depth([s, Slab(0.1, 2.0)]) ≈ 0.4
    @test optical_depth(Slab[]) == 0.0
end

@testset "Empty slab list → transmission = 1" begin
    u_bins = collect(range(0.005, 0.995, length=10))
    T = transmission_factor(Slab[], u_bins)
    @test all(T .== 1.0)
    @test length(T) == length(u_bins)
end

@testset "Single-slab transmission" begin
    μ_test = 0.2
    t_test = 1.5
    u_bins = [0.1, 0.5, 1.0]
    T = transmission_factor([Slab(μ_test, t_test)], u_bins)
    expected = [exp(-μ_test * t_test / u) for u in u_bins]
    @test T ≈ expected
end

@testset "Multi-slab transmission sums optical depths" begin
    s1 = Slab(0.2, 1.0)
    s2 = Slab(0.1, 2.0)
    u_bins = [0.5, 1.0]
    τ = 0.2 * 1.0 + 0.1 * 2.0
    T = transmission_factor([s1, s2], u_bins)
    expected = [exp(-τ / u) for u in u_bins]
    @test T ≈ expected
end

@testset "Transmission monotonic in u (more grazing → more attenuation)" begin
    s = Slab(0.2, 1.0)
    u_bins = collect(range(0.1, 1.0, length=10))
    T = transmission_factor([s], u_bins)
    @test all(diff(T) .> 0.0)   # T increases with u (less grazing)
end
