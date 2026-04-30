# LZ Field Cage + Grids Model for 0νββ Background Simulations

## 1. Model Description

Field cage barrel split into three sub-sources (rings, PTFE, resistors+sensors)
at their correct radial positions. Grid holders modeled as short barrel sources
at the TPC perimeter. All values from bb0nu paper Table I.

## 2. Components

### Barrel (3 sub-sources)

| Component | Material | Mass (kg) | ²³⁸U-late | ²³²Th-late | R_entry (cm) |
|-----------|----------|-----------|-----------|------------|-------------|
| Rings | Ti | 93.0 | 0.35 | 0.24 | 74.3 (FC outer) |
| PTFE walls | PTFE | 184.0 | 0.04 | 0.01 | 72.8 (TPC inner) |
| Resistors+sensors | mixed | 5.08 | ~1350/5.82 | ~2010/1.88 | 74.3 (FC outer) |

### Endcap (grid holders as short barrels)

| Component | Mass (kg) | ²³⁸U-late | ²³²Th-late | Position |
|-----------|-----------|-----------|------------|----------|
| Top holders | 44.15 | 2.63 | 1.46 | R=72.8, H=2.5cm, z=143.1 |
| Bottom holders | 44.15 | 2.63 | 1.46 | R=72.8, H=2.5cm, z=-1.25 |

Wire mesh mass (0.8 kg) dropped as negligible.

## 3. Results (Bi-214, FV [26, 96] cm, r ≤ 39 cm)

| Source | events/yr | Paper (events/yr) | Ratio |
|--------|-----------|-------------------|-------|
| FC rings | 0.50 | 0.30 | 1.7× |
| FC PTFE | 0.14 | 0.14 | 1.0× |
| FC res+sens | 1.70 | 1.39 | 1.2× |
| FC holder top | 0.03 | — | — |
| FC holder bot | 0.38 | — | — |
| **FC holders total** | **0.41** | **0.23** | **1.8×** |
| **FC total** | **2.75** | **2.06** | **1.3×** |

## 4. Source Indices and Run Commands

```bash
# Source indices:
#  7 = fc_rings_Bi214         (barrel, R=74.3)
#  8 = fc_rings_Tl208
#  9 = fc_ptfe_Bi214          (barrel, R=72.8)
# 10 = fc_ptfe_Tl208
# 11 = fc_ressens_Bi214       (barrel, R=74.3)
# 12 = fc_ressens_Tl208
# 13 = fc_holder_top_Bi214    (barrel, R=72.8, H=2.5cm, z=143.1)
# 14 = fc_holder_top_Tl208
# 15 = fc_holder_bot_Bi214    (barrel, R=72.8, H=2.5cm, z=-1.25)
# 16 = fc_holder_bot_Tl208

julia --project=. -t 8 scripts/run_mc.jl --e-visible 10.0 --fv-z-min 26.0 --fv-z-max 96.0 --fv-r-max 39.0 --source-index 7 --output fc_rings_bi214
julia --project=. -t 8 scripts/run_mc.jl --e-visible 10.0 --fv-z-min 26.0 --fv-z-max 96.0 --fv-r-max 39.0 --source-index 9 --output fc_ptfe_bi214
julia --project=. -t 8 scripts/run_mc.jl --e-visible 10.0 --fv-z-min 26.0 --fv-z-max 96.0 --fv-r-max 39.0 --source-index 11 --output fc_ressens_bi214
julia --project=. -t 8 scripts/run_mc.jl --e-visible 10.0 --fv-z-min 26.0 --fv-z-max 96.0 --fv-r-max 39.0 --source-index 13 --output fc_holder_top_bi214
julia --project=. -t 8 scripts/run_mc.jl --e-visible 10.0 --fv-z-min 26.0 --fv-z-max 96.0 --fv-r-max 39.0 --source-index 15 --output fc_holder_bot_bi214
```

## 5. References

1. LZ bb0nu paper — Table I
2. LZ TDR — Table 3.1.1 (geometry)
3. LZ Instrument Paper — Table 1 (grid parameters)
