# src3/geometry.jl — Pure geometric primitives.
#
# Two types: `GCyl` (cylindrical shell) and `GDisk` (head/disc, flat or
# ellipsoidal). Both expose a uniform interface — area, volume, mass,
# uniform-on-surface sampling, inward normal — used by both the analytic
# transport integrals and the per-photon Monte Carlo.
#
# No materials, no activities, no detector parameters. Coordinate frame:
# axis along ẑ; z = 0 at the cathode by project convention.

# ---------------------------------------------------------------------------
# GCyl — cylindrical shell
# ---------------------------------------------------------------------------

"""
    GCyl(R_inner, wall_thickness, z_min, z_max)

Cylindrical shell with axis along ẑ. Inner surface at radius `R_inner`,
outer surface at `R_inner + wall_thickness`. Axial extent `[z_min, z_max]`.
All lengths in cm.
"""
struct GCyl
    R_inner::Float64
    wall_thickness::Float64
    z_min::Float64
    z_max::Float64

    function GCyl(R_inner::Real, wall_thickness::Real,
                  z_min::Real, z_max::Real)
        R_inner > 0 ||
            error("GCyl: R_inner must be positive (got $R_inner)")
        wall_thickness > 0 ||
            error("GCyl: wall_thickness must be positive (got $wall_thickness)")
        z_max > z_min ||
            error("GCyl: z_max ($z_max) must exceed z_min ($z_min)")
        new(Float64(R_inner), Float64(wall_thickness),
            Float64(z_min), Float64(z_max))
    end
end

"Outer radius (cm)."
R_outer(g::GCyl) = g.R_inner + g.wall_thickness

"Axial height z_max − z_min (cm)."
height(g::GCyl) = g.z_max - g.z_min

"Lateral surface area of the inner cylinder, 2π R_inner H (cm²)."
area_inner(g::GCyl) = 2π * g.R_inner * height(g)

"Lateral surface area of the outer cylinder, 2π R_outer H (cm²)."
area_outer(g::GCyl) = 2π * R_outer(g) * height(g)

"Exact shell volume π(R_o² − R_i²) H (cm³)."
volume_shell(g::GCyl) = π * (R_outer(g)^2 - g.R_inner^2) * height(g)

"Volume enclosed by the inner cylinder (cm³)."
volume_inner(g::GCyl) = π * g.R_inner^2 * height(g)

"Mass (kg) of the shell given material density ρ (g/cm³)."
mass(g::GCyl, ρ::Real) = volume_shell(g) * ρ / 1000.0

"""
    sample_inner_surface(rng, g::GCyl) -> (x, y, z)

Sample a point uniformly on the inner cylindrical surface
(area = 2π R_inner H).
"""
function sample_inner_surface(rng::AbstractRNG, g::GCyl)
    φ = 2π * rand(rng)
    z = g.z_min + height(g) * rand(rng)
    (g.R_inner * cos(φ), g.R_inner * sin(φ), z)
end

"""
    inward_normal(g::GCyl, x, y, z) -> (nx, ny, nz)

Unit vector pointing toward the cylinder axis at the point (x, y, z).
The z-coordinate is unused (cylindrical symmetry).
"""
function inward_normal(::GCyl, x::Real, y::Real, ::Real)
    r = sqrt(x*x + y*y)
    r > 0 || error("inward_normal(GCyl): point on axis is degenerate")
    (-x/r, -y/r, 0.0)
end


# ---------------------------------------------------------------------------
# GDisk — head/disc (flat or ellipsoidal)
# ---------------------------------------------------------------------------

"""
    GDisk(R, wall_thickness, z_equator, aspect_ratio, points_outward)

Head attached to a cylinder of equatorial radius `R` at axial position
`z_equator`. The shape parameter `aspect_ratio = n` defines the head depth
along the axis as `R / n`:

  * `aspect_ratio = Inf` → flat disc (depth = 0)
  * `aspect_ratio = 1`   → hemisphere
  * `aspect_ratio = 2`   → 2:1 ellipsoidal head
  * `aspect_ratio = 3`   → 3:1 ellipsoidal head

`points_outward` is `:up` (head bulges in +ẑ) or `:down` (bulges in −ẑ).
All lengths in cm.
"""
struct GDisk
    R::Float64
    wall_thickness::Float64
    z_equator::Float64
    aspect_ratio::Float64
    points_outward::Symbol

    function GDisk(R::Real, wall_thickness::Real, z_equator::Real,
                   aspect_ratio::Real, points_outward::Symbol)
        R > 0 ||
            error("GDisk: R must be positive (got $R)")
        wall_thickness >= 0 ||
            error("GDisk: wall_thickness must be ≥ 0 (got $wall_thickness)")
        aspect_ratio > 0 ||
            error("GDisk: aspect_ratio must be positive (got $aspect_ratio)")
        points_outward in (:up, :down) ||
            error("GDisk: points_outward must be :up or :down (got $points_outward)")
        new(Float64(R), Float64(wall_thickness), Float64(z_equator),
            Float64(aspect_ratio), points_outward)
    end
end

"Head depth along the axis (cm). 0 for a flat disc."
depth(g::GDisk) = isinf(g.aspect_ratio) ? 0.0 : g.R / g.aspect_ratio

"True iff the head is flat (zero depth)."
is_flat(g::GDisk) = isinf(g.aspect_ratio) || g.aspect_ratio == 0.0

"Axial coordinate of the head apex (the point furthest from z_equator)."
function z_apex(g::GDisk)
    sgn = g.points_outward === :up ? +1.0 : -1.0
    g.z_equator + sgn * depth(g)
end

"""
    area_inner(g::GDisk)

Inside surface area of the head (cm²):

  * flat disc: π R²
  * hemisphere (n = 1): 2π R²
  * oblate ellipsoidal head (n > 1): half of the full oblate-spheroid
    surface, with semi-axes a = R (equatorial) and c = R/n (polar)
  * prolate (n < 1, depth > R): half of the prolate spheroid surface
"""
function area_inner(g::GDisk)
    is_flat(g) && return π * g.R^2
    n = g.aspect_ratio
    n ≈ 1.0 && return 2π * g.R^2

    a = g.R
    c = g.R / n
    if c < a
        # Oblate: depth < equatorial radius
        e = sqrt(1.0 - (c/a)^2)
        S_full = 2π * a^2 + π * (c^2 / e) * log((1 + e) / (1 - e))
    else
        # Prolate: depth > equatorial radius (not used for LZ heads)
        e = sqrt(1.0 - (a/c)^2)
        S_full = 2π * a^2 + 2π * a * c * asin(e) / e
    end
    return S_full / 2.0
end

"""
    volume_shell(g::GDisk)

Thin-shell approximation `area_inner × wall_thickness`. Exact for flat
discs; for ellipsoidal heads with t/R ≪ 1 the relative error is O(t/R)
(~1% for LZ heads).
"""
volume_shell(g::GDisk) = area_inner(g) * g.wall_thickness

"Mass (kg) of the head shell given material density ρ (g/cm³)."
mass(g::GDisk, ρ::Real) = volume_shell(g) * ρ / 1000.0

"""
    sample_inner_surface(rng, g::GDisk) -> (x, y, z)

Sample a point uniformly by surface area on the inside of the head.
For a flat disc this is uniform on the disc. For a curved head it uses
rejection sampling against the meridional area element.
"""
function sample_inner_surface(rng::AbstractRNG, g::GDisk)
    sgn = g.points_outward === :up ? +1.0 : -1.0
    if is_flat(g)
        φ = 2π * rand(rng)
        r = g.R * sqrt(rand(rng))
        return (r * cos(φ), r * sin(φ), g.z_equator)
    end
    R = g.R
    d = depth(g)
    # Parametrize meridian by u ∈ [0, π/2]:
    #   r(u) = R sin(u),  z_off(u) = sgn · d · cos(u)
    # Surface element: dA/du = 2π · r · sqrt((dr/du)² + (dz/du)²)
    #                       = 2π · R sin(u) · sqrt(R² cos² u + d² sin² u)
    # Maximum of sqrt(...) is max(R, d); peak of the full element at u = π/2
    M = 2π * R * max(R, d)
    while true
        u = (π/2) * rand(rng)
        elem = 2π * R * sin(u) * sqrt(R^2 * cos(u)^2 + d^2 * sin(u)^2)
        if rand(rng) * M ≤ elem
            φ = 2π * rand(rng)
            r = R * sin(u)
            z_off = sgn * d * cos(u)
            return (r * cos(φ), r * sin(φ), g.z_equator + z_off)
        end
    end
end

"""
    inward_normal(g::GDisk, x, y, z) -> (nx, ny, nz)

Unit vector pointing into the cylinder (away from the head's outside)
at the surface point (x, y, z).
"""
function inward_normal(g::GDisk, x::Real, y::Real, z::Real)
    sgn = g.points_outward === :up ? +1.0 : -1.0
    if is_flat(g)
        return (0.0, 0.0, -sgn)
    end
    R = g.R
    d = depth(g)
    z_off = z - g.z_equator
    # Outward normal of the level set (r/R)² + (z_off/d)² = 1 is
    # ∇ = (x/R², y/R², z_off/d²); inward normal is its negative.
    nx = -x / R^2
    ny = -y / R^2
    nz = -z_off / d^2
    nrm = sqrt(nx^2 + ny^2 + nz^2)
    nrm > 0 || error("inward_normal(GDisk): degenerate point")
    (nx/nrm, ny/nrm, nz/nrm)
end

"""
    path_through_shell(x, y, z, dx, dy, dz, g::GDisk) -> Float64

Path length (cm) of a ray through the head shell, evaluated as
`wall_thickness / |cos α|` where α is the angle between the ray and the
local inward normal at (x, y, z). Exact for flat discs; for ellipsoidal
heads the relative error is O(t/R).
Returns 0 for grazing rays where |cos α| < 1e-10.
"""
function path_through_shell(x::Real, y::Real, z::Real,
                            dx::Real, dy::Real, dz::Real, g::GDisk)
    nx, ny, nz = inward_normal(g, x, y, z)
    cos_α = -(dx*nx + dy*ny + dz*nz)  # angle vs outward normal
    abs(cos_α) < 1e-10 && return 0.0
    g.wall_thickness / abs(cos_α)
end
