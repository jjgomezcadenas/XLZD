using Random, Statistics

@testset "build_geometry — LZ defaults" begin
    params = Params()
    geom = build_geometry(params)

    # TPC: R=72.8, H=145.6 → ~7.16 tonnes
    @test isapprox(geom.R_lxe, 72.8; atol=0.1)
    @test isapprox(geom.L_lxe, 145.6; atol=0.1)

    # Full LXe (TPC + skin): R=80.3 → ~8.7 tonnes
    @test isapprox(geom.R_lxe_full, 80.3; atol=0.1)
    @test isapprox(geom.LXe_mass_kg / 1000.0, 8.7; rtol=0.05)

    # FV from bbonu paper: z ∈ [26, 96], r ≤ 39 → ~967 kg
    @test isapprox(geom.FV_LXe_mass_kg, 967.0; rtol=0.05)
    @test geom.FV_LXe_mass_kg < geom.LXe_mass_kg
end

@testset "BFV — buffer is FV + 1 mm" begin
    params = Params()
    bfv = build_bfv(params)

    @test bfv.z_min == params.fv_z_min_cm - 0.1
    @test bfv.z_max == params.fv_z_max_cm + 0.1

    r_fv = sqrt(params.fv_r2_max_cm2)
    r_bfv = sqrt(bfv.r2_max)
    @test isapprox(r_bfv, r_fv + 0.1; atol=1e-10)
end

@testset "in_bfv / in_fv — boundary tests" begin
    params = Params()
    bfv = build_bfv(params)

    # Center of detector: inside both FV and BFV
    @test in_fv(0.0, 0.0, 60.0, params)
    @test in_bfv(0.0, 0.0, 60.0, bfv)

    # Just outside FV z_max but inside BFV (within 1 mm)
    z_between = params.fv_z_max_cm + 0.05  # 96.05
    @test !in_fv(0.0, 0.0, z_between, params)
    @test in_bfv(0.0, 0.0, z_between, bfv)

    # Just outside BFV z_max
    z_outside = params.fv_z_max_cm + 0.2  # 96.2
    @test !in_bfv(0.0, 0.0, z_outside, bfv)

    # Just outside FV radially but inside BFV (within 1 mm)
    r_fv = sqrt(params.fv_r2_max_cm2)
    r_between = r_fv + 0.05
    @test !in_fv(r_between, 0.0, 60.0, params)
    @test in_bfv(r_between, 0.0, 60.0, bfv)

    # Far outside both
    @test !in_fv(100.0, 0.0, 60.0, params)
    @test !in_bfv(100.0, 0.0, 60.0, bfv)
end

@testset "in_skin / in_active_lxe" begin
    params = Params()
    L = params.geom_L_lxe

    # Point inside field cage (active LXe)
    @test in_active_lxe(0.0, 0.0, 100.0, params, L)
    @test !in_skin(0.0, 0.0, 100.0, params, L)

    # Point in skin region (between FC outer wall and ICV inner wall)
    r_skin = (params.geom_R_skin_inner + params.geom_R_skin_outer) / 2.0
    @test in_skin(r_skin, 0.0, 100.0, params, L)
    @test !in_active_lxe(r_skin, 0.0, 100.0, params, L)

    # Point at FC boundary
    r_fc = params.geom_R_skin_inner - 0.1  # just inside FC
    @test in_active_lxe(r_fc, 0.0, 100.0, params, L)

    # Point outside ICV inner wall
    r_outside = params.geom_R_skin_outer + 1.0
    @test !in_skin(r_outside, 0.0, 100.0, params, L)
    @test !in_active_lxe(r_outside, 0.0, 100.0, params, L)
end

@testset "sample_entry_point — barrel" begin
    source = SourceConfig(
        label="test_barrel", E_MeV=2.448, entry=:barrel,
        R_entry=74.3, H_entry=159.35, z_min_entry=70.0
    )
    rng = MersenneTwister(42)
    L = 297.0

    for _ in 1:100
        x, y, z = sample_entry_point(rng, source, L)
        r = sqrt(x^2 + y^2)
        @test isapprox(r, 74.3; atol=1e-10)
        @test z >= 70.0
        @test z <= 70.0 + 159.35
    end
end

@testset "sample_entry_point — endcap" begin
    source = SourceConfig(
        label="test_endcap", E_MeV=2.448, entry=:endcap,
        R_entry=74.3
    )
    rng = MersenneTwister(42)
    L = 297.0

    n_top = 0; n_bot = 0
    for _ in 1:1000
        x, y, z = sample_entry_point(rng, source, L)
        r = sqrt(x^2 + y^2)
        @test r <= 74.3 + 1e-10
        @test z == 0.0 || z == L
        if z == 0.0; n_bot += 1; else; n_top += 1; end
    end
    # Roughly 50/50 split
    @test 400 < n_top < 600
    @test 400 < n_bot < 600
end

@testset "sample_u — flat" begin
    source = SourceConfig(
        label="test", E_MeV=2.448, entry=:barrel,
        R_entry=74.3, angular=:flat
    )
    rng = MersenneTwister(42)

    N = 10000
    us = [sample_u(rng, source) for _ in 1:N]
    @test all(0.0 .≤ us .≤ 1.0)
    # Mean should be ~0.5 for uniform
    @test isapprox(mean(us), 0.5; atol=0.02)
end

@testset "sample_u — shaped" begin
    source = SourceConfig(
        label="test", E_MeV=2.448, entry=:barrel,
        R_entry=74.3, angular=:shaped,
        a=-0.633, b=3.16, u_min=0.3, norm=0.994
    )
    rng = MersenneTwister(42)

    N = 10000
    us = [sample_u(rng, source) for _ in 1:N]
    @test all(us .>= 0.3 - 1e-10)
    @test all(us .<= 1.0 + 1e-10)
    # Mean should be > 0.5 (distribution peaked toward u=1)
    @test mean(us) > 0.6
end

@testset "sample_entry_direction — barrel inward" begin
    source = SourceConfig(
        label="test", E_MeV=2.448, entry=:barrel,
        R_entry=74.3, H_entry=159.35, z_min_entry=70.0,
        angular=:flat
    )
    rng = MersenneTwister(42)
    L = 297.0

    for _ in 1:100
        x, y, z = sample_entry_point(rng, source, L)
        dx, dy, dz = sample_entry_direction(rng, source, x, y, z, L)

        # Direction is a unit vector
        norm = sqrt(dx^2 + dy^2 + dz^2)
        @test isapprox(norm, 1.0; atol=1e-10)

        # Dot product with outward radial should be negative (pointing inward)
        r = sqrt(x^2 + y^2)
        dot_radial = (x * dx + y * dy) / r
        @test dot_radial < 0.0
    end
end

@testset "sample_entry_direction — endcap inward" begin
    source = SourceConfig(
        label="test", E_MeV=2.448, entry=:endcap,
        R_entry=74.3, angular=:flat
    )
    rng = MersenneTwister(42)
    L = 297.0

    for _ in 1:100
        x, y, z = sample_entry_point(rng, source, L)
        dx, dy, dz = sample_entry_direction(rng, source, x, y, z, L)

        norm = sqrt(dx^2 + dy^2 + dz^2)
        @test isapprox(norm, 1.0; atol=1e-10)

        # Bottom endcap (z=0): dz should be positive (inward = up)
        # Top endcap (z=L): dz should be negative (inward = down)
        if z < 1.0
            @test dz > 0.0
        else
            @test dz < 0.0
        end
    end
end

@testset "sample_companion_energy" begin
    rng = MersenneTwister(42)
    N = 10000
    energies = [sample_companion_energy(rng) for _ in 1:N]

    # All returned energies are valid
    valid = Set([0.583, 0.860, 0.763, 0.0])
    @test all(e ∈ valid for e in energies)

    # Check approximate branching ratios (10% tolerance)
    n_583 = count(e -> e ≈ 0.583, energies)
    n_860 = count(e -> e ≈ 0.860, energies)
    n_763 = count(e -> e ≈ 0.763, energies)
    n_none = count(e -> e ≈ 0.0, energies)

    @test isapprox(n_583 / N, 0.85; atol=0.02)
    @test isapprox(n_860 / N, 0.12; atol=0.02)
    @test isapprox(n_763 / N, 0.02; atol=0.01)
    @test isapprox(n_none / N, 0.01; atol=0.01)
end

@testset "path_to_cylinder_exit" begin
    R = 100.0

    # Photon at center heading along +x → exits at R
    d = path_to_cylinder_exit(0.0, 0.0, 100.0, 1.0, 0.0, 0.0, R, 0.0, 200.0)
    @test isapprox(d, R; atol=1e-6)

    # Photon at center heading along +z → exits at z_max = 200
    d = path_to_cylinder_exit(0.0, 0.0, 50.0, 0.0, 0.0, 1.0, R, 0.0, 200.0)
    @test isapprox(d, 150.0; atol=1e-6)

    # Photon at center heading along -z → exits at z_min = 0
    d = path_to_cylinder_exit(0.0, 0.0, 50.0, 0.0, 0.0, -1.0, R, 0.0, 200.0)
    @test isapprox(d, 50.0; atol=1e-6)

    # Photon near wall heading outward → short exit distance
    d = path_to_cylinder_exit(99.0, 0.0, 100.0, 1.0, 0.0, 0.0, R, 0.0, 200.0)
    @test isapprox(d, 1.0; atol=0.1)

    # Axial exit dominates
    d = path_to_cylinder_exit(0.0, 0.0, 199.0, 0.0, 0.0, 1.0, R, 0.0, 200.0)
    @test isapprox(d, 1.0; atol=1e-6)

    # Extended volume: z_min = -20, photon heading down exits at z = -20
    d = path_to_cylinder_exit(0.0, 0.0, 0.0, 0.0, 0.0, -1.0, R, -20.0, 200.0)
    @test isapprox(d, 20.0; atol=1e-6)
end
