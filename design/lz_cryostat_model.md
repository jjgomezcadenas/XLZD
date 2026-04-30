# LZ Cryostat Model for 0νββ Background Simulations

## 1. Model Description

Two concentric Ti cylinders (ICV, OCV) with real wall thicknesses from
TDR Table 5.4.1. LXe skin (6 cm) is tracked in the Julia MC, not
pre-attenuated in Python. Extra Ti mass (flanges, support) assigned to
bottom endcap (heavily shielded by dome LXe).

All activities from LZ bb0nu paper Table I.

## 2. Geometry

| Parameter | ICV | OCV |
|-----------|-----|-----|
| R (cm) | 83.0 | 91.5 |
| H (cm) | 259.0 | 304.0 |
| Wall (mm) | 9.0 | 7.0 |
| Top head | 2:1, 8 mm | 2:1, 9 mm |
| Bottom head | 3:1, 12 mm | 2:1, 15 mm |
| Shell mass (kg) | 401.9 | 385.7 |
| Head mass (kg) | 249.2 | 392.9 |
| Total (kg) | 651.1 | 778.6 |

Total vessel geometric mass: 1429.7 kg
Extra Ti (flanges, ports, support): 1160.3 kg → bottom endcap
**Grand total: 2590 kg** (bb0nu paper)

ICV inner radius: 82.1 cm (entry point for cryostat gammas)

## 3. Radioactivity

| Source | Mass (kg) | ²³⁸U-late (mBq/kg) | ²³²Th-late (mBq/kg) |
|--------|-----------|---------------------|----------------------|
| Ti vessel | 2590 | 0.08 | 0.22 |
| MLI insulation | 13.8 | 11.1 | 7.79 |

## 4. Source Positions

| Source | Type | Position | Notes |
|--------|------|----------|-------|
| Barrel (shell) | barrel at R=82.1 | z ∈ [0, 145.6] | Skin tracked in MC |
| Head top | endcap_top at z=190 | In gas above PMTs | Propagates through gas to LXe surface |
| Head bottom | endcap_bottom at z=-44 | In LXe below dome | Traverses ~44 cm invisible LXe |

## 5. Python Transport

Gammas attenuated through Ti walls only (no LXe skin). LXe skin is
tracked with full Compton physics in the Julia MC. Skin deposits ≥ 100 keV
trigger skin veto.

## 6. Results (Bi-214, FV [26, 96] cm, r ≤ 39 cm)

| Source | events/yr | Paper (events/yr) |
|--------|-----------|-------------------|
| Barrel (Ti + MLI) | 1.18 | — |
| Head top | 0.03 | — |
| Head bottom | 0.02 | — |
| **Cryostat total** | **1.23** | **0.80** |

Ratio: 1.5×. MLI (upper limit activity) accounts for most of the excess.

## 7. Run Commands

```bash
# Source indices (from run_mc.jl):
#  1 = cryo_barrel_Bi214
#  2 = cryo_head_top_Bi214
#  3 = cryo_head_bot_Bi214
#  4 = cryo_barrel_Tl208
#  5 = cryo_head_top_Tl208
#  6 = cryo_head_bot_Tl208

julia --project=. -t 8 scripts/run_mc.jl --e-visible 10.0 --fv-z-min 26.0 --fv-z-max 96.0 --fv-r-max 39.0 --source-index 1 --output cryo_barrel_bi214
julia --project=. -t 8 scripts/run_mc.jl --e-visible 10.0 --fv-z-min 26.0 --fv-z-max 96.0 --fv-r-max 39.0 --source-index 2 --output cryo_head_top_bi214
julia --project=. -t 8 scripts/run_mc.jl --e-visible 10.0 --fv-z-min 26.0 --fv-z-max 96.0 --fv-r-max 39.0 --source-index 3 --output cryo_head_bot_bi214
```

## 8. References

1. LZ bb0nu paper — Table I
2. LZ TDR — Table 5.4.1 (wall thicknesses), Table 3.1.1 (geometry)
