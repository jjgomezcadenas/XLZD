# LZ Field Cage Model

<!-- Numbers hand-copied from data/lz_fc_barrels.csv and
     data/lz_fc_grids.csv as of 2026-05-03. If you change a CSV,
     ask Claude to resync this file (and design/tikz/fieldcage.tex). -->

Six effective sources sitting **inside** the LXe (no Ti vessel between
them and the active volume), built directly from two CSVs. Four are
smeared cylindrical shells over the field-cage axial extent; two are
annular slabs at the gate / cathode planes representing the grid
holders. The bottom-shield grid is dropped from the model (γ paths to
the FV blocked by cathode + RFR).

Specific activities and masses come from the LZ bb0nu paper Table I.
Geometry is summarized in `docs/LZ_detector_summary.md` §2-§5.

## 1. Component inventory

### 1.1 Barrel components

**Source: `data/lz_fc_barrels.csv`** (4 rows, all share z ∈ [-13.75, 145.6] cm)

| name | material | R_inner (cm) | mass (kg) | Bi-214 (mBq/kg) | Tl-208 (mBq/kg) | description                                     |
|------|----------|--------------|-----------|-----------------|-----------------|-------------------------------------------------|
| FCRN | Ti       | 74.3         | 93.0      | 0.35            | 0.24            | Field-cage Ti rings (57 drift + 8 RFR)          |
| FCRS | SS       | 74.3         | 0.06      | 1350            | 2010            | Field-shaping resistors (alumina ceramic; SS proxy) |
| FCSE | SS       | 74.3         | 5.02      | 5.82            | 1.88            | TPC sensors (loop antennas + sensors)           |
| FCPT | PTFE     | 72.8         | 184.0     | 0.04            | 0.01            | PTFE reflector walls (5 mm panels)              |

### 1.2 Grid holders

**Source: `data/lz_fc_grids.csv`** (2 rows, both SS, same annular footprint)

| name | R_in (cm) | R_out (cm) | z_face (cm) | normal_sign | mass (kg) | Bi-214 (mBq/kg) | Tl-208 (mBq/kg) | description                                              |
|------|-----------|------------|-------------|-------------|-----------|-----------------|-----------------|----------------------------------------------------------|
| FCTG | 72.8      | 80.3       | 145.6       | −1          | 44.55     | 2.63            | 1.46            | Anode+gate collapsed at gate plane; opens upward         |
| FCBG | 72.8      | 80.3       |   0.0       | +1          | 22.275    | 2.63            | 1.46            | Cathode only; bottom-shield dropped; opens downward      |

`normal_sign` encodes which way the inward emission normal points (+1 = +ẑ,
−1 = −ẑ). The slab itself sits *behind* the LXe-facing face (FCTG above
the gate plane in the gas/anode region; FCBG below the cathode plane in
the RFR), so γ leaving the LXe-facing face go into the active LXe drift.

### 1.3 Mass and source-mode rationale

- **FCBG = ¼ of bb0nu's 89.1 kg lump** (4 grids of comparable mass; cathode-only).
- **FCTG = ½ of 89.1 kg** (anode + gate grids collapsed; their 1.3 cm gap is small vs. the distance to the FV top at z = 96 cm).
- **Bottom-shield grid dropped** — its z position (below RFR) puts it behind the cathode in optical path; γ from there cannot reach the FV without traversing the cathode + RFR LXe, and the contribution is absorbed into model uncertainty.

Combined FC mass: 93 + 0.06 + 5.02 + 184 + 44.55 + 22.275 = **348.9 kg**.

## 2. Back-derived geometry

The CSV stores `mass_kg` directly (the engineering input); the loader
back-derives the wall thickness (PCyl) or slab height (PAnnularDisk)
at startup so the slab self-shielding integral has the correct optical
depth without an explicit thickness column.

For PCyl barrels: `t = mass·1000 / (ρ · 2π · R · H)`, with H = 159.35 cm
(z_max − z_min for all four barrel rows).

For PAnnularDisk grids: `H = mass·1000 / (ρ_SS · π · (R_out² − R_in²))`.

| name | ρ (g/cm³) | back-derived wall t (PCyl) or H (slab) | unit |
|------|-----------|----------------------------------------|------|
| FCRN | 4.510 (Ti)   | 0.277  | cm |
| FCRS | 7.930 (SS)   | 1.0×10⁻⁴ | cm |
| FCSE | 7.930 (SS)   | 0.00866 | cm |
| FCPT | 2.200 (PTFE) | 1.147  | cm |
| FCTG | 7.930 (SS)   | 1.558  | cm (axial slab) |
| FCBG | 7.930 (SS)   | 0.779  | cm (axial slab) |

## 3. Self-shielding (optical depth at 2.448 MeV)

μ·t per source. Materials use the NIST XCOM tables in `data/`:
`nist_ti.csv` (Ti), `nist_fe.csv` (Fe as proxy for SS), `nist_teflon.csv` (PTFE).

| name | μ_lin at 2.448 MeV (cm⁻¹) | μ·t  | inward fraction (vs. ½) |
|------|---------------------------|------|--------------------------|
| FCRN | 0.196 (Ti)               | 0.048 | 0.90 |
| FCRS | 0.312 (SS)               | 3.0×10⁻⁵ | 0.99 |
| FCSE | 0.312 (SS)               | 0.0027 | 0.98 |
| FCPT | 0.0840 (PTFE)            | 0.0964 | 0.84 |
| FCTG | 0.312 (SS)               | 0.485 | 0.56 |
| FCBG | 0.312 (SS)               | 0.243 | 0.70 |

The "inward fraction" is the ratio of the source's `total_per_yr` (at
its LXe-facing exit surface, after self-shielding) to the no-shielding
limit of `½ × total γ rate`. FCRS, FCSE are essentially un-shielded
(thin shells); FCTG and FCBG are the only sources with non-trivial
self-shielding because their slab cross-section is real solid SS.

The SS uses **Fe as a chemical proxy** (Fe is 70% of SS-304 by mass;
Cr/Ni Z values are nearly identical, so μ/ρ matches to <2%). Density
is set to ρ_SS = 7.93 g/cm³ (SS-304), not ρ_Fe.

## 4. γ rates per source

4π γ-production rates, computed as `mass_kg × spec_act × 1e-3 × BR ×
seconds_per_year`. For Bi-214: BR_γ = 1.55%. For Tl-208: BR_chain =
35.9% (²³²Th-late → ²⁰⁸Tl), BR_γ ≈ 100%. For the Tl-208 cascade
companion: ~99% (lumped at 583 keV).

| name | mass (kg) | Bi-214 γ/yr | Tl-208 γ/yr | Tl-208c γ/yr |
|------|-----------|-------------|--------------|---------------|
| FCRN | 93.0   | 1.592e+04 | 2.529e+05 | 2.504e+05 |
| FCRS | 0.06   | 3.962e+04 | 1.366e+06 | 1.353e+06 |
| FCSE | 5.02   | 1.429e+04 | 1.069e+05 | 1.058e+05 |
| FCPT | 184.0  | 3.600e+03 | 2.085e+04 | 2.064e+04 |
| FCTG | 44.55  | 5.731e+04 | 7.369e+05 | 7.295e+05 |
| FCBG | 22.275 | 2.866e+04 | 3.684e+05 | 3.647e+05 |

Inward γ/yr per effective source (after self-shielding) is computed
by `make_gamma_source` and stored in `EffectiveSource.total_per_yr`;
see the FC test output in `test/test_fc_effective_sources3.jl` for
the actual numbers.

## 5. Effective sources

6 components × 3 isotopes (Bi-214, Tl-208, Tl-208c) = **18 effective
sources**. Each has a single contribution (the source PObject itself)
with empty downstream slab list (no Ti to traverse — sources sit in
the LXe).

Region symbols (used by `src3/sampling.jl` to dispatch entry-point
geometry):

- `:fc_barrel`  for FCRN, FCRS, FCSE, FCPT (PCyl emission on the inner cylinder)
- `:fc_annular` for FCTG, FCBG (PAnnularDisk emission on the LXe-facing face)

Effective-source names follow `<COMPONENT>_<ISOTOPE>` exactly: `FCRN_Bi214`,
`FCBG_Tl208c`, etc. The Tl-208c sources are not user-selectable; the
driver attaches them automatically to the corresponding `_Tl208` main
source for the companion-veto step.

## 6. Sampling primitives

Two new entry-point samplers in `src3/sampling.jl`, each reading
the source geometry directly from `eff.contributions[1].source.producer`:

- `sample_fc_barrel_entry(rng, g::GCyl, cdf, u_bins)` — uniform on the
  inner cylindrical surface at radius `g.R_inner`, z ∈ [g.z_min, g.z_max];
  inward normal = −r̂. Used for FCRN/FCRS/FCSE/FCPT.
- `sample_fc_annular_entry(rng, p::PAnnularDisk, cdf, u_bins)` — uniform
  by area on the annular face at `z = p.z_face` (sample r² uniformly
  in [R_in², R_out²], φ uniform in [0, 2π]); inward normal =
  `p.normal_sign · ẑ`. Used for FCTG/FCBG.

The angular distribution `dN/du(u = cos θ_inward)` is sampled via
inverse-CDF from each effective source's pre-computed `dNdu` array,
just as for the cryostat sources.

## 7. Cross-check vs. LZ bb0nu paper

The original spec (with bottom-grid mass split as ½, not ¼) reproduced
the LZ paper to ~1.3× (`design/lz_fieldcage_model.md` previous version).
With the cathode-only FCBG model, the smoke-test totals at 1e6 samples
(default cuts: σ_rms/E = 0.7 %, ROI = ±1 FWHM, FV [26,96] × r ≤ 39 cm):

| component  | smoke γ/yr (1e6 samples, 8 threads) |
|------------|--------------------------------------|
| FCRN_Bi214 | 0.90 |
| FCRS_Bi214 | 2.12 |
| FCSE_Bi214 | 0.79 |
| FCPT_Bi214 | 0.090 |
| FCTG_Bi214 | 0.032 |
| FCBG_Bi214 | 0.181 |
| **Bi-214 total** | **4.11** events/yr |

vs. spec (full-mass FCBG) **2.75** and bb0nu paper **2.06** events/yr.
Tl-208 contributions need production-scale statistics (≥1e9 samples) to
resolve the Compton-tail fraction inside the ROI window.

## 8. References

1. LZ bb0nu paper — Table I (specific activities, masses)
2. LZ TDR — Table 3.6.3 (grid wire parameters), Table 9.2.7 (grid masses)
3. LZ Detector Summary — `docs/LZ_detector_summary.md` §2-§5
4. LZ Instrument Paper — Table 1 (grid voltages, wire ø, pitch, counts)
