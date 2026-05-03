# LZ Cryostat Model

<!-- Numbers hand-copied from data/lz_cryo_*.csv as of 2026-05-03.
     If you change a CSV, ask Claude to resync this file (and its
     companion design/tikz/cryostat.tex). -->

Two concentric Ti vessels (Outer Cryostat Vessel + Inner Cryostat
Vessel) with ellipsoidal heads, a 6 cm LXe skin in the gap between the
ICV and the field cage, and ~1138 kg of additional Ti structures
(flanges, ports, support legs, weld reinforcements, internal anchors)
distributed along both vessels. The skin LXe is **tracked** in the
Julia MC (with a 100 keV deposit veto), not pre-attenuated in Python.

Activities are from the LZ bb0nu paper Table I. Geometry comes from
TDR Table 5.4.1 (wall thicknesses) and Table 3.1.1 (overall dimensions).

## 1. Vessel geometry

**Source: `data/lz_cryo_geometry.csv`** (1 row per vessel)

| element | R (cm) | H_total (cm) | H_cyl (cm) | t_wall (mm) | t_top (mm) | t_bot (mm) | top AR | bot AR |
|---------|--------|--------------|------------|-------------|------------|------------|--------|--------|
| OCV     | 91.5   | 304.0        | 212.5      | 7.0         | 9.0        | 15.0       | 2:1    | 2:1    |
| ICV     | 83.0   | 259.0        | 189.83     | 9.0         | 8.0        | 12.0       | 2:1    | 3:1    |

`H_total = H_cyl + R/top_ar + R/bot_ar` (cylindrical body + two
ellipsoidal head depths). Aspect ratio `n` means head depth = `R / n`
along the axis (n=2 → semi-axis ratio 2:1 oblate; n=3 → 3:1 oblate).
ICV inner radius = 82.1 cm (= `R - t_wall = 83.0 - 0.9`), used as the
γ-entry surface for the cryostat sources.

## 2. Vessel mass per piece

Derived from CSV row × Ti density (4.510 g/cm³). Barrels use exact
shell volume `π(R_o² − R_i²)·H`; heads use the thin-shell approximation
`area_inner × t_wall` (relative error O(t/R) ≈ 1% for LZ heads).

| element | shell (kg) | heads (kg) | total (kg) |
|---------|-----------|------------|-----------|
| OCV     | 385.7     | 392.9      | 778.6     |
| ICV     | 401.9     | 249.2      | 651.1     |
| **vessels total** | — | — | **1429.7** |

## 3. Extras (flanges, ports, legs, anchors)

**Source: `data/lz_cryo_extras.csv`** (one row per Ti element). Each
row is an annular shell with `count` multiplicity and an `include_in_mc`
flag controlling whether it is exposed as an MC source.

### 3.1 MC-active extras (driven by `include_in_mc = true`)

| name                              | category        | R_in (cm) | R_out (cm) | z_min (cm) | z_max (cm) | count | mass (kg) |
|-----------------------------------|-----------------|-----------|------------|------------|------------|-------|-----------|
| OCV_top_reinforcing_ring          | flange          | 92.2      | 96.2       | 174.0      | 179.0      | 1     | 53.4      |
| OCV_upper_body_flange_pair        | flange          | 92.2      | 105.0      | 97.0       | 103.0      | 1     | 214.6     |
| OCV_lower_body_flange_pair        | flange          | 92.2      | 105.0      | -13.0      | -7.0       | 1     | 214.6     |
| OCV_bottom_support_flange         | flange          | 92.2      | 92.8       | -42.0      | -34.0      | 1     | 12.6      |
| ICV_main_flange_pair              | flange          | 83.0      | 90.5       | 127.0      | 133.0      | 1     | 110.6     |
| OCV_top_YBe_recess                | port            | 15.0      | 18.0       | 180.0      | 190.0      | 1     | 14.0      |
| OCV_top_conduit_ports             | port            | 10.0      | 12.0       | 180.0      | 198.0      | 2     | 22.4      |
| OCV_top_tie_bar_ports             | port            | 7.0       | 8.5        | 180.0      | 198.0      | 3     | 17.8      |
| OCV_HV_port                       | port            | 92.2      | 93.0       | 44.0       | 56.0       | 1     | 25.2      |
| OCV_cathode_side_port             | port            | 92.2      | 92.7       | 85.0       | 95.0       | 1     | 13.1      |
| OCV_bottom_umbilical_port         | port            | 12.7      | 15.7       | -90.0      | -80.0      | 1     | 12.1      |
| **MC-active extras total**        |                 |           |            |            |            |       | **710.4** |

### 3.2 Inactive extras (`include_in_mc = false`, mass kept for accounting)

`OCV_bottom_small_port`, `ICV_weir_drains`, `ICV_bottom_anchor_port_flanges`,
`OCV_support_legs`, `ICV_stiffening_ring`, `ICV_coldhead_fins`,
`OCV_top_tie_bar_mount_ring`, `OCV_bottom_leg_pad_ring`,
`OCV_upper_joint_weld_reinforcement`, `OCV_lower_joint_weld_reinforcement`,
`ICV_bottom_anchor_block_ring`. Combined mass ≈ 432.0 kg; either
heavily shielded by surrounding structures or contributing negligibly
to the FV background (kept on the books, not exposed as sources).

### 3.3 Mass totals by category (all 22 extras, MC and non-MC)

| category         | mass (kg) |
|------------------|-----------|
| flange           | 605.8     |
| internal_anchor  | 208.0     |
| leg              | 122.4     |
| port             | 120.6     |
| fin              | 56.6      |
| stiffener        | 25.5      |
| **extras total** | **1138.9** |

**Grand total Ti**: 1429.7 (vessels) + 1138.9 (extras) = **2568.6 kg**, vs
the bb0nu paper's 2590 kg (1% under-count from CSV-vs-paper roundoff).

## 4. Surface sources (non-Ti)

**Source: `data/lz_cryo_surface_sources.csv`** (one row per non-Ti
surface). Currently a single MLI insulation layer attached to the ICV
body.

| name                      | attached_to | material | total_mass (kg) | act_U238_late (mBq/kg) | act_Th232_late (mBq/kg) |
|---------------------------|-------------|----------|-----------------|------------------------|--------------------------|
| Cryostat_insulation_MLI   | ICV_body    | MLI      | 13.8            | 11.1                   | 7.79                     |

The MLI is modeled as a `PSurface` (zero-thickness, no self-shielding).
Its activity is added to whichever effective source matches its
attachment face — here `:ICV_body` → effective source CB.

## 5. Specific activities

Late-chain Ti specific activities from the bb0nu paper Table I, applied
uniformly to every Ti element (vessels + extras):

```
TI_BB0NU_U238_LATE_MBQKG  = 0.08    # mBq/kg, ²³⁸U-late chain (Bi-214 line)
TI_BB0NU_TH232_LATE_MBQKG = 0.22    # mBq/kg, ²³²Th-late chain (Tl-208 line)
```

MLI carries its own measured specific activities (much higher than Ti;
see Table 4).

## 6. γ rates

Total γ-production rate per chain at full 4π, summed over all MC-active
PObjects. Computed as `mass_kg × spec_act_mBqkg × 1e-3 × BR × seconds_per_year`.

Branching ratios:
- Bi-214 → 2.448 MeV γ: BR = 1.55%
- Tl-208 → 2.615 MeV γ: BR_chain = 35.9% (²³²Th-late → ²⁰⁸Tl), BR_γ ≈ 100%
- Tl-208 cascade companion: ~99% of decays yield a coincident γ (lumped at 583 keV)

| family | mass (kg) | Bi-214 γ/yr | Tl-208 γ/yr |
|--------|-----------|-------------|--------------|
| Ti vessels + extras (full 2568 kg) | 2568.6 | 1.005e+05 | 6.402e+06 |
| MLI insulation                     | 13.8   | 7.495e+04 | 1.218e+06 |

(Inward γ/yr per effective source — after each source's slab
self-shielding — is computed by `make_gamma_source` and aggregated by
`build_effective_sources`. See per-source numbers in
`test/test_pobjects3.jl` console output.)

## 7. Effective sources

Three regions × 2 isotopes (Bi-214, Tl-208) + 3 cascade-companion
sources (Tl-208c at 583 keV) = **9 effective sources**, defined as the
inward γ flux at the ICV inner surface after each source's downstream
Ti slab attenuation:

| name        | region        | isotope | exit surface              | downstream Ti slab (per contribution) |
|-------------|---------------|---------|---------------------------|----------------------------------------|
| CB_Bi214    | barrel        | Bi-214  | ICV inner cylinder, R=82.1 cm | none for ICV-attached, ICV body wall (9 mm) for OCV-attached |
| CTH_Bi214   | endcap_top    | Bi-214  | ICV top head (2:1, t=8 mm) | same pattern, ICV top wall for OCV-attached |
| CBH_Bi214   | endcap_bottom | Bi-214  | ICV bottom head (3:1, t=12 mm) | ICV bottom wall for OCV-attached |
| CB_Tl208    | barrel        | Tl-208  | as above                  | as above                               |
| CTH_Tl208   | endcap_top    | Tl-208  | as above                  | as above                               |
| CBH_Tl208   | endcap_bottom | Tl-208  | as above                  | as above                               |
| CB_Tl208c   | barrel        | Tl-208c | as above                  | as above                               |
| CTH_Tl208c  | endcap_top    | Tl-208c | as above                  | as above                               |
| CBH_Tl208c  | endcap_bottom | Tl-208c | as above                  | as above                               |

Each effective source's contribution list (which physical objects feed
which effective source) is tabulated in `src3/effective_sources.jl` —
search for `cb_contribs`, `cth_contribs`, `cbh_contribs`. The MLI
surface is appended to the effective source matching its `attached_to`
face (here CB).

The Tl-208c sources are not user-selectable; the driver attaches them
automatically to the corresponding `_Tl208` main source for the
companion-veto step.

## 8. Sampling primitives

The cryostat sources sample γ entry from the **ICV inner exit surface**
(γ have already traversed any downstream Ti at this point). Three
geometric primitives in `src3/sampling.jl`:

- `sample_barrel_entry`        — uniform on cylinder R = 82.1 cm,
                                  z ∈ [z_RFR_bottom, z_gate]
- `sample_endcap_entry(:endcap_top)`
                                — uniform-by-area on the 2:1 oblate
                                  ellipsoidal inner head; equator at
                                  z = z_ICV_top − 82.1/2
- `sample_endcap_entry(:endcap_bottom)`
                                — uniform-by-area on the 3:1 oblate
                                  inner head; equator at
                                  z = z_LXe_bottom + 82.1/3

The angular distribution `dN/du(u = cos θ_inward)` is sampled via
inverse-CDF from each effective source's pre-computed `dNdu` array.

## 9. References

1. LZ bb0nu paper — Table I (specific activities, total Ti mass)
2. LZ TDR — Table 3.1.1 (vessel dimensions), Table 5.4.1 (wall thicknesses)
3. LZ Detector Summary — `docs/LZ_detector_summary.md`
4. Cryostat MC specification — `design/xlzd_code.md`
