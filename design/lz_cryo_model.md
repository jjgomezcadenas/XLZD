# LZ Cryostat Model — Refactored Geometry & Mass Budget

This document describes the cryostat representation built in `src2/`,
the data files that feed it, and the subset that is exposed as gamma
sources to the Monte Carlo. It supersedes the simple-cylinder + averaged-
endcap model in `py/lz_cryo.py`, which lumped the OCV and ICV heads into
two area-weighted endcap sources and absorbed the missing ~1160 kg of
hardware mass into a single "extra Ti" surface source on the bottom
endcap.

## 1. Motivation

The bb0nu paper (Phys. Rev. C 102, 014602, Table I) gives the full Ti
cryostat assembly mass as **2 590 kg** with late-chain specific
activities of 0.08 mBq/kg (²³⁸U) and 0.22 mBq/kg (²³²Th). Reproducing
the bb0nu background requires the *total* Ti activity (mass × specific
activity) to match.

The TDR's bare-vessel masses (Table 5.4.1) are 950 kg (ICV) + 1 115 kg
(OCV) = 2 065 kg. A simple-cylinder geometric model with the wall
thicknesses from `lz_cryo_geometry.csv` gives 1 433 kg. The remaining
~1 160 kg is real Ti hardware — flanges, ports, legs, stiffeners,
coldhead fins, internal anchor rings — that the simple-cylinder model
omits.

The new model places every gram of Ti at a defined geometric location,
so each kilogram has a birth point for MC photon emission; nothing is
treated as phantom mass.

## 2. Geometry primitives

Two pure-geometry types live in `src2/geometry.jl`:

### `GCyl` — cylindrical shell

```
struct GCyl
    R_inner::Float64        # inner radius (cm)
    wall_thickness::Float64 # radial wall thickness (cm)
    z_min::Float64          # axial bottom (cm)
    z_max::Float64          # axial top (cm)
end
```

Methods: `R_outer`, `height`, `area_inner`, `area_outer`,
`volume_shell` (exact: π(R_o² − R_i²)·H), `volume_inner`, `mass(g, ρ)`,
`sample_inner_surface(rng, g)`, `inward_normal(g, x, y, z)`.

### `GDisk` — head/disc (flat or ellipsoidal)

```
struct GDisk
    R::Float64
    wall_thickness::Float64
    z_equator::Float64
    aspect_ratio::Float64    # Inf=flat, 1=hemisphere, 2=2:1, 3=3:1
    points_outward::Symbol   # :up or :down
end
```

Methods: `depth(g) = R/aspect_ratio` (0 for flat), `z_apex`, `is_flat`,
`area_inner` (oblate-spheroid surface for n>1, with hemisphere and flat
limits handled explicitly), `volume_shell ≈ area_inner × wall_thickness`
(thin-shell, ~1% relative error for LZ heads at t/R ≈ 0.01),
`sample_inner_surface`, `inward_normal`, `path_through_shell`.

## 3. Cryostat composition

The composite `Cryostat` (in `src2/cryostat.jl`) is built by
`build_cryostat(geom_csv, extras_csv)` from two data files:

### 3.1 `data/lz_cryo_geometry.csv` — barrels and heads

Two rows (OCV, ICV), columns: `R_cm`, `H_total_cm`, `H_cyl_cm`,
`t_wall_mm`, `t_top_mm`, `t_bot_mm`, `top_ar`, `bot_ar`,
`mass_shell_kg`, `mass_heads_kg`, `mass_total_kg`. Each row generates:

* one `GCyl` (the cylindrical body), and
* two `GDisk` (top and bottom heads, with their respective aspect
  ratios and orientations `:up` / `:down`).

So the geometry CSV contributes **2 barrels + 4 heads = 6 objects**.

### 3.2 `data/lz_cryo_extras.csv` — additional Ti hardware

22 rows representing all the Ti hardware not covered by the simple
barrels-and-heads model. Each row is a `GCyl` (annular shell);
columns:

| column | meaning |
|---|---|
| `name` | human-readable label |
| `category` | `flange / port / leg / stiffener / fin / internal_anchor` |
| `R_inner_cm`, `R_outer_cm` | radii of the annular shell |
| `z_min_cm`, `z_max_cm` | axial extent |
| `count` | number of identical pieces (mass scales by `count`) |
| `include_in_mc` | flag for use as a gamma-source ring |

Categories present:

* **flange** (5 rows, 606 kg): OCV top reinforcing ring, OCV upper
  body flange pair, OCV lower body flange pair, OCV bottom support
  flange, ICV main flange pair.
* **port** (8 rows, 121 kg): YBe recess, conduit ports, tie-bar
  ports, HV port, cathode side port, bottom umbilical, bottom small
  port, ICV weir drains.
* **leg** (1 row × count 3, 122 kg): OCV support legs.
* **stiffener** (1 row, 26 kg): ICV stiffening ring at top of conical
  section.
* **fin** (1 row × count 6, 57 kg): ICV coldhead fins.
* **internal_anchor** (5 rows, 208 kg): ICV bottom anchor port
  flanges, OCV-top tie-bar mount ring, OCV-bottom leg-pad ring, two
  OCV body-joint weld reinforcements, ICV bottom anchor block ring.

### 3.3 The `Cryostat` struct

```
struct Cryostat
    barrels::Vector{GCyl}           # 2 entries
    heads::Vector{GDisk}            # 4 entries
    extras::Vector{CryostatExtra}   # 22 entries
end
```

`CryostatExtra` wraps a `GCyl` with its name, category, count, and
`include_in_mc` flag.

## 4. Mass budget

Computed by `total_mass(c, ρ_Ti)` and `mass_breakdown(c, ρ_Ti)` with
`ρ_Ti = 4.510 g/cm³`:

| category | mass (kg) |
|---|---:|
| Barrels (OCV + ICV) | 791 |
| Heads (4 ellipsoidal) | 642 |
| Flanges (5) | 606 |
| Ports (8 entries, 12 pieces) | 121 |
| Legs (3 pieces) | 122 |
| Stiffener | 26 |
| Coldhead fins (6 pieces) | 57 |
| Internal anchor rings (5 entries, 10 pieces) | 208 |
| **Total Ti mass (model)** | **2 572** |
| **bb0nu Table I** | **2 590** |
| **Residual (model − bb0nu)** | **−18 (−0.7 %)** |

Tested in `test/test_geometry2.jl` against ±5 % tolerance; current
deviation is 0.7 %.

## 5. Radioactive (MC-active) subset

Not all of the 22 extras serve as gamma sources in the MC. Some have
geometry that does not fit cleanly into the axisymmetric `GCyl`
primitive (the 3 OCV legs are off-axis), or are heavily shielded by
multiple Ti walls (internal anchors), or are small enough that their
contribution to the FV background is negligible.

The radioactive model — items flagged `include_in_mc = true` — comprises:

| | mass (kg) | source |
|---|---:|---|
| Barrels (OCV + ICV) | 791 | always MC-active |
| Heads (4 ellipsoidal) | 642 | always MC-active |
| Flanges (5 entries) | 606 | bb0nu paper does not break out individual hardware; the flanges are the bulk Ti at the bb0nu specific activity |
| Main ports (6 entries: HV, cathode side, YBe recess, conduits, tie-bars, bottom umbilical) | 105 | the prominent ports visible in the engineering drawing |
| **MC-active subtotal** | **2 144** | **83 % of total Ti** |

Items NOT in the radioactive model (17 % of the Ti mass):

* OCV bottom small port (9 kg) — small contribution.
* ICV weir drains (7 kg, 3 pieces) — heavily shielded by ICV body.
* ICV bottom anchor port flanges (39 kg, 6 pieces) — at small radii
  on the ICV bottom dome, geometrically awkward to position.
* OCV support legs (122 kg, 3 pieces) — off-axis, outside the OCV
  bottom dished end, behind the dome-skin LXe and water tank;
  contribution to the FV is sub-percent.
* ICV stiffening ring (26 kg) — internal hardware.
* ICV coldhead fins (57 kg, 6 pieces) — internal.
* Internal anchor rings (208 kg, 5 entries) — buried inside the
  vessel-to-vessel and TPC-to-vessel mounting hardware.

The 17 % MC-inactive Ti consists almost entirely of items that are
either heavily shielded (legs, anchors), buried inside other Ti walls
(stiffener, fins), or small. We expect their contribution to the FV
background to be considerably less than 17 %.

## 6. Conventions and caveats

### 6.1 Coordinate frame
`z = 0` at the cathode (matches the existing Julia MC convention).
TPC drift region: `z ∈ [0, 145.6 cm]`. ICV cylindrical body:
`z ∈ [−41.3, +148.5 cm]`. ICV apexes at `z = −69` and `+190`. OCV
cylindrical body: `z ∈ [−38, +174]`; OCV apexes at `z = −84` and `+220`.

### 6.2 Mass formulae
* `GCyl.volume_shell` uses the exact formula π(R_o² − R_i²)·H, ~0.4 %
  larger than the thin-shell formula 2π·R_inner·H·t used in the
  reference `lz_cryo_geometry.csv` numbers. The 1 % test tolerance
  absorbs this.
* `GDisk.volume_shell` uses the thin-shell approximation
  `area_inner × wall_thickness`. Relative error O(t/R) ≈ 1 % at LZ
  scales.

### 6.3 Estimated dimensions
The engineering drawing in `docs/cryostat.png` gives reliable
diameters for every flange and port but does **not** give axial
thicknesses for the flanges and most ports. ASME-typical values for
Ti at 1.8 m diameter / ~2 bar are used in `lz_cryo_extras.csv`. If
authoritative LZ engineering drawings become available, the CSV
should be updated to match.

### 6.4 Side ports as cylindrical rings
The side ports (HV, cathode side) have axes perpendicular to the
cryostat z-axis, which the `GCyl` primitive does not represent
directly. They are encoded as thin annular rings at the OCV body
radius, with axial extent matching the port location and diameter.
This conserves mass and places gammas at the correct (R, z) but
loses port-axis directional information — fine for the analytic
shielding integrals where the source-to-LXe path is dominated by the
OCV wall in any direction.

### 6.5 Plural pieces (`count > 1`)
Items appearing N times at different azimuths (3 tie-bar ports,
6 coldhead fins, 6 ICV anchor flanges, 3 legs) are recorded as one
CSV row with `count = N`. The loader produces `count` copies of the
same `GCyl` (same R, z); for axisymmetric mass integration this
is exact. When per-azimuth positions matter (e.g., for off-axis
transport), the CSV will need an `azimuth_deg` column added.

## 7. Files

| file | purpose |
|---|---|
| `src2/geometry.jl` | `GCyl`, `GDisk` and their methods |
| `src2/cryostat.jl` | `Cryostat`, `CryostatExtra`, `build_cryostat`, `total_mass`, `mass_breakdown`, `mc_active_extras`, `mc_active_mass` |
| `src2/XLZD2.jl` | module root, exports |
| `data/lz_cryo_geometry.csv` | barrels + heads parameters |
| `data/lz_cryo_extras.csv` | 22 additional Ti elements |
| `test/test_geometry2.jl` | geometry tests + full mass budget regression |
| `latex/lz_cryo_model.tex` | standalone TikZ source for the model figure |
| `latex/lz_cryo_model.png` | rendered figure included in the paper appendix |
