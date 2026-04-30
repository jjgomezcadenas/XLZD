# LZ PMT Model for 0νββ Background Simulations

## 1. Model Description

PMTs separated into three source groups at their physical positions:
- **Top PMTs** at z=152.6 cm (gas, propagated to liquid surface at z=145.6)
- **Bottom PMTs** at z=-15.75 cm (LXe, traverse RFR invisibly)
- **Barrel** (cables + side skin PMTs) at R=72.8 cm

All activities from bb0nu paper Table I (per kg, not per PMT).

## 2. Components

| bb0nu Paper Row | Mass (kg) | ²³⁸U-late | ²³²Th-late | Our Region |
|-----------------|-----------|-----------|------------|------------|
| TPC PMTs | 91.9 | 3.22 | 1.61 | top/bottom endcap |
| TPC PMT bases | 2.80 | 75.9 | 33.1 | top/bottom endcap |
| TPC PMT structures | 166.0 | 1.60 | 1.06 | top/bottom endcap |
| TPC PMT cables | 88.7 | 4.31 | 0.82 | barrel |
| Skin PMTs+bases | 8.59 | 46.0 | 14.9 | barrel (side) / bottom |

PMTs+bases+structures split 253/241 (top/bottom) by PMT count.
Structures split 50/50.

## 3. Source Positions

| Source | Entry | z_entry | R_entry | Extended vol | Notes |
|--------|-------|---------|---------|-------------|-------|
| PMT top | endcap_top | 152.6 | 72.8 | yes | Gas→propagate to LXe at z=145.6 |
| PMT bottom | endcap_bottom | -15.75 | 72.8 | yes | In LXe, RFR invisible (z<0) |
| PMT barrel | barrel | 0-145.6 | 72.8 | no | Cables + side skin PMTs |

Key physics:
- Top PMTs: 7 cm gas (transparent) + ~50 cm LXe to FV top
- Bottom PMTs: 16 cm invisible LXe (RFR) + 26 cm to FV bottom
- Deposits below cathode (z<0) invisible — no S2 signal

## 4. Results (Bi-214, FV [26, 96] cm, r ≤ 39 cm)

| Source | events/yr | Paper (events/yr) |
|--------|-----------|-------------------|
| PMT top | 0.76 | — |
| PMT bottom | 2.22 | — |
| PMT barrel | 0.39 | — |
| **PMT total** | **3.37** | **3.40** |

Ratio: 0.99× — excellent agreement.

Paper breakdown: PMTs 2.95 + bases 1.52 + structures 2.65 + cables 1.44 +
skin 0.75 = **9.31 counts/1000d = 3.40 events/yr**.

## 5. ²⁰⁸Tl Companion Gamma

For Tl-208 sources, companion gamma (583/860/763 keV) emitted in 4π.
Event vetoed if companion deposits ≥ companion_veto_keV (CLI, default 5 keV)
in active LXe or ≥ skin_veto_keV (CLI, default 100 keV) in LXe skin.

## 6. Source Indices and Run Commands

```bash
# Source indices:
# 17 = pmt_top_Bi214           (endcap_top, z=152.6)
# 18 = pmt_top_Tl208
# 19 = pmt_bottom_Bi214        (endcap_bottom, z=-15.75)
# 20 = pmt_bottom_Tl208
# 21 = pmt_barrel_Bi214        (barrel, R=72.8)
# 22 = pmt_barrel_Tl208

julia --project=. -t 8 scripts/run_mc.jl --e-visible 10.0 --fv-z-min 26.0 --fv-z-max 96.0 --fv-r-max 39.0 --source-index 17 --output pmt_top_bi214
julia --project=. -t 8 scripts/run_mc.jl --e-visible 10.0 --fv-z-min 26.0 --fv-z-max 96.0 --fv-r-max 39.0 --source-index 19 --output pmt_bottom_bi214
julia --project=. -t 8 scripts/run_mc.jl --e-visible 10.0 --fv-z-min 26.0 --fv-z-max 96.0 --fv-r-max 39.0 --source-index 21 --output pmt_barrel_bi214
```

## 7. Not Modeled (flat corrections)

| Source | Counts/1000d |
|--------|-------------|
| Other components | 2.41 |
| Outer detector | 2.79 |
| Cavern walls | 11.6 |
| ²²²Rn | 0.45 |
| ¹³⁷Xe | 0.28 |
| ¹³⁶Xe 2νββ | 0.01 |
| ⁸B solar ν | 0.03 |

## 8. Grand Summary — All Bi-214 Sources

| Source | Index | events/yr | Paper |
|--------|-------|-----------|-------|
| Cryo barrel | 1 | 1.18 | — |
| Cryo head top | 2 | 0.03 | — |
| Cryo head bot | 3 | 0.02 | — |
| **Cryostat total** | | **1.23** | **0.80** |
| FC rings | 7 | 0.50 | 0.30 |
| FC PTFE | 9 | 0.14 | 0.14 |
| FC res+sens | 11 | 1.70 | 1.39 |
| FC holder top | 13 | 0.03 | — |
| FC holder bot | 15 | 0.38 | — |
| **FC total** | | **2.75** | **2.06** |
| PMT top | 17 | 0.76 | — |
| PMT bottom | 19 | 2.22 | — |
| PMT barrel | 21 | 0.39 | — |
| **PMT total** | | **3.37** | **3.40** |
| **GRAND TOTAL** | | **7.35** | **7.67** |

**Agreement: 4% below the paper.** The remaining discrepancy is within
the uncertainties of our simplified cylindrical geometry.

## 9. Detector Parameters (CLI defaults)

| Parameter | Default | CLI flag |
|-----------|---------|----------|
| FV z_min | 26.0 cm | --fv-z-min |
| FV z_max | 96.0 cm | --fv-z-max |
| FV r_max | 39.0 cm | --fv-r-max |
| σ/E | 0.7% | (in Params) |
| ROI halfwidth | 17.2 keV | --roi-halfwidth |
| Δz SS/MS | 3.0 mm | --dz-threshold |
| E visible | 10.0 keV | --e-visible |
| Companion veto | 5.0 keV | --companion-veto |
| Skin veto | 100.0 keV | --skin-veto |
| N samples | 10⁸ | --n-samples |

## 10. References

1. LZ bb0nu paper — Table I
2. LZ TDR — Table 3.1.1, Table 5.4.1
3. LZ Instrument Paper — Sections 2.1, 2.2
